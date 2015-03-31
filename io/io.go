package io

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"log"
	"os"

	"unsafe"

	"github.com/phil-mansfield/gotetra/geom"
)

// TODO: swtich from logging error statements to returning error codes.

const (
	// Endianness used by default when writing catalogs. Catalogs of any
	// endianness can be read.
	DefaultEndiannessFlag int32 = -1
)

// This is a terrible idea and shouldn't exist.
type Particle struct {
	Xs, Vs geom.Vec
	Id int64
}

// CatalogHeader describes meta-information about the current catalog.
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
    3 - (sheet.Header) Header dile containing meta-information about the
        sheet fragment.
    4 - ([][3]float32) Contiguous block of x, y, z coordinates. Given in Mpc.
    5 - ([][3]float32) Contiguous block of v_x, v_y, v_z coordinates.
 */
type SheetHeader struct {
	Cosmo CosmologyHeader
	Count, CountWidth int64
	SegmentWidth, GridWidth, GridCount int64
	Idx, Cells int64

	Mass float64
	TotalWidth float64

	Origin geom.Vec
	Width geom.Vec
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
//
// This call signature, and espeically the Particle type are all a consequence
// of soem shockingly poor early design decisions.
func ReadGadgetParticlesAt(
	path string,
	order binary.ByteOrder,
	floatBuf []float32,
	intBuf []int64,
	ps []Particle,
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

	if int64(len(floatBuf)) != 3*h.Count {
		panic(fmt.Sprintf(
			"Incorrect length for float buffer. Found %d, expected %d",
			len(floatBuf), 3*h.Count,
		))
	} else if int64(len(intBuf)) != h.Count {
		panic(fmt.Sprintf(
			"Incorrect length for int buffer. Found %d, expected %d",
			len(intBuf), h.Count,
		))
	} else if int64(len(ps)) != h.Count {
		panic(fmt.Sprintf(
			"Incorrect length for Particle buffer. Found %d, expected %d",
			len(ps), h.Count,
		))
	}

	_ = readInt32(f, order)
	binary.Read(f, order, floatBuf)
	_ = readInt32(f, order)

	for i := range ps {
		ps[i].Xs[0] = float32(gh.WrapDistance(float64(floatBuf[3*i+0])))
		ps[i].Xs[1] = float32(gh.WrapDistance(float64(floatBuf[3*i+1])))
		ps[i].Xs[2] = float32(gh.WrapDistance(float64(floatBuf[3*i+2])))
	}

	_ = readInt32(f, order)
	binary.Read(f, order, floatBuf)
	_ = readInt32(f, order)

	rootA := float32(math.Sqrt(float64(gh.Time)))
	for i := range ps {
		ps[i].Vs[0] = floatBuf[3*i+0] * rootA
		ps[i].Vs[1] = floatBuf[3*i+1] * rootA
		ps[i].Vs[2] = floatBuf[3*i+2] * rootA
	}

	_ = readInt32(f, order)
	binary.Read(f, order, intBuf)
	_ = readInt32(f, order)

	for i := range ps {
		ps[i].Id = intBuf[i]
	}
}

// ReadGadget reads the gadget particle catalog located at the given location
// and written with the given endianness. Its header and particle sequence
// are returned in a standardized format.
func ReadGadget(
	path string, order binary.ByteOrder,
) (*CatalogHeader, []Particle) {

	h := ReadGadgetHeader(path, order)
	floatBuf := make([]float32, 3 * h.Count)
	intBuf := make([]int64, h.Count)
	ps := make([]Particle, h.Count)
	ReadGadgetParticlesAt(path, order, floatBuf, intBuf, ps)
	return h, ps
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

	binary.Read(f, order, hdBuf)
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
		return fmt.Errorf("Position buffer has length %d, but file %s has %d " + 
			"vectors.", len(xsBuf), file, h.GridCount)
	}

	// Go to block 4 in the file.
	// The file pointer should already be here, but let's just be safe, okay?
	f.Seek(int64(4 + 4 + int(unsafe.Sizeof(SheetHeader{}))), 0)
	if err := binary.Read(f, order, xsBuf); err != nil { return err }

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

	// Go to block 5 in the file.
	f.Seek(int64(4 + 4 + int(unsafe.Sizeof(SheetHeader{})) +
		3 * 4 * len(vsBuf)), 0)
	if err := binary.Read(f, order, vsBuf); err != nil { return err }
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
	if err = binary.Write(f, order, int32(unsafe.Sizeof(SheetHeader{}))); err != nil {
		log.Fatalf(err.Error())
	}
	if err = binary.Write(f, order, h); err != nil {
		log.Fatalf(err.Error())
	}
	if err = binary.Write(f, order, xs); err != nil {
		log.Fatalf(err.Error())
	}
	if err = binary.Write(f, order, vs); err != nil {
		log.Fatalf(err.Error())
	}
}

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
