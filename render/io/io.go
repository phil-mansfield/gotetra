package io

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"log"
	"os"
	"reflect"

	"unsafe"

	"github.com/phil-mansfield/gotetra/render/geom"
)

// TODO: swtich from logging error statements to returning error codes.

const (
	// Endianness used by default when writing catalogs. Catalogs of any
	// endianness can be read.
	DefaultEndiannessFlag int32 = -1
)

// CatalogHeader describes meta-information about the current catalog. This
// serves as a standardized version of the different header formats.
type CatalogHeader struct {
	Cosmo CosmologyHeader

	Mass       float64 // Mass of one particle
	Count      int64   // Number of particles in catalog
	TotalCount int64   // Number of particles in all catalogs
	CountWidth int64   // Number of particles "on one side": TotalCount^(1/3)

	Idx        int64   // Index of catalog: x-major ordering is used
	GridWidth  int64   // Number of gird cells "on one side"
	Width      float64 // Width of the catalog's bounding box
	TotalWidth float64 // Width of the sim's bounding box
}

/* 
The binary format used for phase sheets is as follows:
    |-- 1 --||-- 2 --||-- ... 3 ... --||-- ... 4 ... --||-- ... 5 ... --|

    1 - (int32) Flag indicating the endianness of the file. 0 indicates a big
        endian byte ordering and -1 indicates a little endian byte order.
    2 - (int32) Size of a Header struct. Should be checked for consistency.
    3 - (sheet.Header) Header file containing meta-information about the
        sheet fragment.
    4 - ([][3]float32) Contiguous block of x, y, z coordinates. Given in Mpc.
    5 - ([][3]float32) Contiguous block of v_x, v_y, v_z coordinates.
 */
type SheetHeader struct {
	Cosmo CosmologyHeader
	Count, CountWidth int64
	SegmentWidth, GridWidth, GridCount int64 // GridWidth = SegmentWidth + 1
	Idx, Cells int64

	Mass float64
	TotalWidth float64

	Origin, Width geom.Vec
	VelocityOrigin, VelocityWidth geom.Vec
}

// CosmologyHeader contains information describing the cosmological
// context in which the simulation was run.
type CosmologyHeader struct {
	Z      float64
	OmegaM float64
	OmegaL float64
	H100   float64
}


// gadgetHeader is the formatting for meta-information used by Gadget 2.
type gadgetHeader struct {
	NPart                                     [6]uint32
	Mass                                      [6]float64
	Time, Redshift                            float64
	FlagSfr, FlagFeedback                     int32
	NPartTotal                                [6]uint32
	FlagCooling, NumFiles                     int32
	BoxSize, Omega0, OmegaLambda, HubbleParam float64
	FlagStellarAge, HashTabSize               int32

	Padding [88]byte
}

// readInt32 returns single 32-bit interger from the given file using the
// given endianness.
func readInt32(r io.Reader, order binary.ByteOrder) int32 {
	var n int32
	if err := binary.Read(r, order, &n); err != nil {
		panic(err)
	}
	return n
}

// Standardize returns a Header that corresponds to the source
// Gadget 2 header.
func (gh *gadgetHeader) Standardize() *CatalogHeader {
	h := &CatalogHeader{}

	h.Count = int64(gh.NPart[1] + gh.NPart[0]<<32)
	h.TotalCount = int64(gh.NPartTotal[1] + gh.NPartTotal[0]<<32)
	h.Mass = float64(gh.Mass[1])
	h.TotalWidth = float64(gh.BoxSize)
	h.Width = -1.0

	h.Cosmo.Z = gh.Redshift
	h.Cosmo.OmegaM = gh.Omega0
	h.Cosmo.OmegaL = gh.OmegaLambda
	h.Cosmo.H100 = gh.HubbleParam

	return h
}

// WrapDistance takes a value and interprets it as a position defined within
// a periodic domain of width h.BoxSize.
func (h *gadgetHeader) WrapDistance(x float64) float64 {
	if x < 0 {
		return x + h.BoxSize
	} else if x >= h.BoxSize {
		return x - h.BoxSize
	}
	return x
}

// ReadGadgetHeader reads a Gadget catalog and returns a standardized
// gotetra containing its information.
func ReadGadgetHeader(path string, order binary.ByteOrder) *CatalogHeader {
	f, err := os.Open(path)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	gh := &gadgetHeader{}

	_ = readInt32(f, order)
	binary.Read(f, binary.LittleEndian, gh)
	h := gh.Standardize()

	return h
}

// ReadGadgetParticlesAt reads a Gadget file and writes all the particles within
// it to the given particle buffer, ps. floatBuf and intBuf are used internally.
// The length of all three buffers must be equal to  the number of particles in
// the catalog.
func ReadGadgetParticlesAt(
	path string,
	order binary.ByteOrder,
	xs, vs []geom.Vec,
	ids []int64,
) {
	f, err := os.Open(path)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	gh := &gadgetHeader{}

	_ = readInt32(f, order)
	binary.Read(f, binary.LittleEndian, gh)
	_ = readInt32(f, order)

	h := gh.Standardize()

	if int64(len(xs)) != h.Count {
		panic(fmt.Sprintf(
			"Incorrect length for xs buffer. Found %d, expected %d",
			len(xs), h.Count,
		))
	} else if int64(len(vs)) != h.Count {
		panic(fmt.Sprintf(
			"Incorrect length for vs buffer. Found %d, expected %d",
			len(vs), h.Count,
		))
	} else if int64(len(ids)) != h.Count {
		panic(fmt.Sprintf(
			"Incorrect length for int buffer. Found %d, expected %d",
			len(ids), h.Count,
		))
	}

	_ = readInt32(f, order)
	readVecAsByte(f, order, xs)
	_ = readInt32(f, order)
	_ = readInt32(f, order)
	readVecAsByte(f, order, vs)
	_ = readInt32(f, order)
	_ = readInt32(f, order)
	readInt64AsByte(f, order, ids)
	//_ = readInt32(f, order)
	
	// Convert gadget velocities out of code units.
	rootA := float32(math.Sqrt(float64(gh.Time)))
	for i := range xs {
		for j := 0; j < 3; j++ {
			xs[i][j] = float32(gh.WrapDistance(float64(xs[i][j])))
			vs[i][j] = vs[i][j] * rootA
		}
	}
}

// ReadGadget reads the gadget particle catalog located at the given location
// and written with the given endianness. Its header and particle sequence
// are returned in a standardized format.
func ReadGadget(
	path string, order binary.ByteOrder,
) (hd *CatalogHeader, xs, vs []geom.Vec, ids []int64) {

	hd = ReadGadgetHeader(path, order)
	xs = make([]geom.Vec,  hd.Count)
	vs = make([]geom.Vec,  hd.Count)
	ids = make([]int64, hd.Count)

	ReadGadgetParticlesAt(path, order, xs, vs, ids)
	return hd, xs, vs, ids
}

// endianness is a utility function converting an endianness flag to a
// byte order.
func endianness(flag int32) binary.ByteOrder {
	if flag == 0 {
		return binary.LittleEndian
	} else if flag == -1 {
		return binary.BigEndian
	} else {
		panic("Unrecognized endianness flag.")
	}
}

// readSheetHeaderAt reads the header of a sheet file into the buffer hdBuf.
func readSheetHeaderAt(
	file string, hdBuf *SheetHeader,
) (*os.File, binary.ByteOrder, error) {
	f, err := os.OpenFile(file, os.O_RDONLY, os.ModePerm)
	if err != nil { return nil, binary.LittleEndian, err }

	// order doesn't matter for this read, since flags are symmetric.
	order := endianness(readInt32(f, binary.LittleEndian))

	headerSize := readInt32(f, order)
	if headerSize != int32(unsafe.Sizeof(SheetHeader{})) {
		return nil, binary.LittleEndian, 
		fmt.Errorf("Expected catalog.SheetHeader size of %d, found %d.",
			unsafe.Sizeof(SheetHeader{}), headerSize,
		)
	}

	_, err = f.Seek(4 + 4, 0)
	if err != nil { return nil, binary.LittleEndian, err }

	err = binary.Read(f, order, hdBuf)
	if err != nil { return nil, binary.LittleEndian, err }

	hdBuf.Count = hdBuf.CountWidth*hdBuf.CountWidth*hdBuf.CountWidth
	return f, order, nil
}

// ReadHeaderAt reads the header in the given file into the target Header.
func ReadSheetHeaderAt(file string, hdBuf *SheetHeader) error {
	f, _, err := readSheetHeaderAt(file, hdBuf)
	if err != nil { return err }
	if err = f.Close(); err != nil { return err }
	return nil
}

// ReadPositionsAt reads the velocities in the given file into a buffer.
func ReadSheetPositionsAt(file string, xsBuf []geom.Vec) error {
	h := &SheetHeader{}
	f, order, err := readSheetHeaderAt(file, h)
	if err != nil { return nil }

	if h.GridCount != int64(len(xsBuf)) {
		return fmt.Errorf("Position buffer has length %d, but file %s has %d "+ 
			"vectors.", len(xsBuf), file, h.GridCount)
	}

	// Go to block 4 in the file.
	// The file pointer should already be here, but let's just be safe, okay?
	f.Seek(int64(4 + 4 + int(unsafe.Sizeof(SheetHeader{}))), 0)
	if err := readVecAsByte(f, order, xsBuf); err != nil { return err }

	if err := f.Close(); err != nil { return err }
	return nil
}

// ReadVelocitiesAt reads the velocities in the given file into a buffer.
func ReadSheetVelocitiesAt(file string, vsBuf []geom.Vec) error {
	h := &SheetHeader{}
	f, order, err := readSheetHeaderAt(file, h)
	if err != nil { return err }
	if h.GridCount != int64(len(vsBuf)) {
		return fmt.Errorf("Velocity buffer has length %d, but file %s has %d " + 
			"vectors.", len(vsBuf), file, h.GridCount)
	}

	f.Seek(int64(4 + 4 + int(unsafe.Sizeof(SheetHeader{})) +
		3 * 4 * len(vsBuf)), 0)
	if err := readVecAsByte(f, order, vsBuf); err != nil { return err }
	if err := f.Close(); err != nil { return err }
	return nil
}

// Write writes a grid of position and velocity vectors to a file, defined
// by the given header.
func WriteSheet(file string, h *SheetHeader, xs, vs []geom.Vec) {
	if int(h.GridCount) != len(xs) {
		log.Fatalf("Header count %d for file %s does not match xs length, %d",
			h.GridCount, file, len(xs))
	} else if int(h.GridCount) != len(vs) {
		log.Fatalf("Header count %d for file %s does not match vs length, %d",
			h.GridCount, file, len(xs))
	} else if h.GridWidth*h.GridWidth*h.GridWidth != h.GridCount {
		log.Fatalf("Header CountWidth %d doesn't match Count %d",
			h.GridWidth, h.GridCount)
	}

	f, err := os.Create(file)
	endiannessFlag := int32(0)
	order := endianness(endiannessFlag)

	if err = binary.Write(f, order, endiannessFlag); err != nil {
		log.Fatalf(err.Error())
	}
	if err = binary.Write(
		f, order, int32(unsafe.Sizeof(SheetHeader{})),
		); err != nil {
		log.Fatalf(err.Error())
	}

	if err = binary.Write(f, order, h); err != nil {
		log.Fatalf(err.Error())
	}
	if err = writeVecAsByte(f, order, xs); err != nil {
		log.Fatalf(err.Error())
	}
	if err = writeVecAsByte(f, order, vs); err != nil {
		log.Fatalf(err.Error())
	}
}

// CellBounds returns the bounds  of a given sheet aligned to the boundaries of
// the voxels with a given granularity.
func (hd *SheetHeader) CellBounds(cells int) *geom.CellBounds {
	cb := &geom.CellBounds{}
	cellWidth := hd.TotalWidth / float64(cells)

	for j := 0; j < 3; j++ {
		cb.Origin[j] = int(
			math.Floor(float64(hd.Origin[j]) / cellWidth),
		)
		cb.Width[j] = 1 + int(
			math.Floor(float64(hd.Origin[j] + hd.Width[j]) / cellWidth),
		)

		cb.Width[j] -= cb.Origin[j]
	}

	return cb
}

// readVecAsByte reads vectors from rd without needing intermediary interface
// allocations.
func readVecAsByte(rd io.Reader, end binary.ByteOrder, buf []geom.Vec) error {
	bufLen := len(buf)

	hd := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
	hd.Len *= 12
	hd.Cap *= 12
	
	byteBuf := *(*[]byte)(unsafe.Pointer(&hd))
	_, err := rd.Read(byteBuf)
	if err != nil { return err }

	if !isSysOrder(end) {
		for i := 0; i < bufLen * 3; i++ {
			for j := 0; j < 2; j++ {
				idx1, idx2 := i*4 + j, i*4 + 3 - j
				byteBuf[idx1], byteBuf[idx2] = byteBuf[idx2], byteBuf[idx1]
			}
		}
	}

	hd.Len /= 12
	hd.Cap /= 12

	return nil
}

// writeVecAsByte writes an array of vectors to wr without needing to make
// intermediate interface allocations.
func writeVecAsByte(wr io.Writer, end binary.ByteOrder, buf []geom.Vec) error {
	bufLen := len(buf)

	hd := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
	hd.Len *= 12
	hd.Cap *= 12
	
	byteBuf := *(*[]byte)(unsafe.Pointer(&hd))

	if !isSysOrder(end) {
		for i := 0; i < bufLen * 3; i++ {
			for j := 0; j < 2; j++ {
				idx1, idx2 := i*4 + j, i*4 + 3 - j
				byteBuf[idx1], byteBuf[idx2] = byteBuf[idx2], byteBuf[idx1]
			}
		}
	}

	_, err := wr.Write(byteBuf)
	if err != nil { return err }

	hd.Len /= 12
	hd.Cap /= 12

	return nil
}

// readInt64AsByte reads int64 values from rd.
func readInt64AsByte(rd io.Reader, end binary.ByteOrder, buf []int64) error {
	bufLen := len(buf)

	hd := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
	hd.Len *= 8
	hd.Cap *= 8
	
	byteBuf := *(*[]byte)(unsafe.Pointer(&hd))
	_, err := rd.Read(byteBuf)
	if err != nil { return err }

	if !isSysOrder(end) {
		for i := 0; i < bufLen; i++ {
			for j := 0; j < 4; j++ {
				idx1, idx2 := i*8 + j, i*8 + 7 - j
				byteBuf[idx1], byteBuf[idx2] = byteBuf[idx2], byteBuf[idx1]
			}
		}
	}

	hd.Len /= 8
	hd.Cap /= 8

	return nil
}

// isSysOrder returns true if end is the system endianness order.
func isSysOrder(end binary.ByteOrder) bool {
	buf32 := []int32{1}

	hd := *(*reflect.SliceHeader)(unsafe.Pointer(&buf32))
	hd.Len *= 4
	hd.Cap *= 4

	buf8 := *(*[]int8)(unsafe.Pointer(&hd))
	if buf8[0] == 1 {
		return binary.LittleEndian == end
	}
	return binary.BigEndian == end
}
