package sheet

import (
	"encoding/binary"
	"io"
	"log"
	"os"
	"unsafe"

	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/catalog"
)

/* 
The binary format used is as follows:
    |-- 1 --||-- 2 --||-- ... 3 ... --||-- ... 4 ... --||-- ... 5 ... --|

    1 - (int32) Flag indicating the endianness of the file. 0 indicates a big
        endian byte ordering and -1 indicates a little endian byte order.
    2 - (int32) Size of a Header struct. Should be checked for consistency.
    3 - (sheet.Header) Header dile containing meta-information about the
        sheet fragment.
    4 - ([][3]float32) Contiguous block of x, y, z coordinates. Given in Mpc.
    5 - ([][3]float32) Contiguous block of v_x, v_y, v_z coordinates.
 */

type Header struct {
	Cosmo catalog.CosmologyHeader
	Count, CountWidth int64
	SegmentWidth, GridWidth, GridCount int64
	Idx, Cells int64

	Mass float64
	TotalWidth float64
	Mins geom.Vec
	Widths geom.Vec
}

func readHeaderAt(file string, hdBuf *Header) (*os.File, binary.ByteOrder) {
	f, err := os.OpenFile(file, os.O_RDONLY, os.ModePerm)
	if err != nil {
		log.Fatal(err.Error())
	}

	// order doesn't matter for this read, since flags are symmetric.
	order := endianness(readInt32(f, binary.LittleEndian))
	headerSize := readInt32(f, order)
	if headerSize != int32(unsafe.Sizeof(Header{})) {
		log.Fatalf("Expected sheet.Header size of %d, found %d.",
			unsafe.Sizeof(Header{}), headerSize)
	}

	_, err = f.Seek(4 + 4, 0)
	if err != nil {
		log.Fatalf(err.Error())
	}

	binary.Read(f, order, hdBuf)

	return f, order
}

// ReadHeaderAt reads the header in the given file into the target Header.
func ReadHeaderAt(file string, hdBuf *Header) {
	f, _ := readHeaderAt(file, hdBuf)
	if err := f.Close(); err != nil {
		log.Fatalf(err.Error())
	}
}

// ReadPositionsAt reads the velocities in the given file into a buffer.
func ReadPositionsAt(file string, xsBuf []geom.Vec) {
	h := &Header{}
	f, order := readHeaderAt(file, h)
	if h.Count != int64(len(xsBuf)) {
		log.Fatalf("Position buffer has length %d, but file %s has %d vectors.",
			len(xsBuf), file, h.Count)
	}

	// Go to block 4 in the file.
	// The file pointer should already be here, but let's just be safe, okay?
	f.Seek(int64(4 + 4 + int(unsafe.Sizeof(Header{}))), 0)
	if err := binary.Read(f, order, xsBuf); err != nil {
		log.Fatalf(err.Error())
	}

	if err := f.Close(); err != nil {
		log.Fatalf(err.Error())
	}
}

// ReadVelocitiesAt reads the velocities in the given file into a buffer.
func ReadVelocitiesAt(file string, vsBuf []geom.Vec) {
	h := &Header{}
	f, order := readHeaderAt(file, h)
	if h.Count != int64(len(vsBuf)) {
		log.Fatalf("Velocity buffer has length %d, but file %s has %d vectors.",
			len(vsBuf), file, h.Count)
	}

	// Go to block 5 in the file.
	f.Seek(int64(4 + 4 + int(unsafe.Sizeof(Header{})) + 3 * 4 * len(vsBuf)), 0)
	if err := binary.Read(f, order, vsBuf); err != nil {
		log.Fatalf(err.Error())
	}

	if err := f.Close(); err != nil {
		log.Fatalf(err.Error())
	}
}

// Write writes a grid of position and velocity vectors to a file, defined
// by the given header.
func Write(file string, h *Header, xs, vs []geom.Vec) {
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
	if err = binary.Write(f, order, int32(unsafe.Sizeof(Header{}))); err != nil {
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

// readInt32 returns single 32-bit interger from the given file using the
// given endianness.
func readInt32(r io.Reader, order binary.ByteOrder) int32 {
	var n int32
	if err := binary.Read(r, order, &n); err != nil {
		log.Fatalf(err.Error())
	}
	return n
}

// endianness is a utility function converting an endianness flag to a
// byte order.
func endianness(flag int32) binary.ByteOrder {
	if flag == 0 {
		return binary.LittleEndian
	} else if flag == -1 {
		return binary.BigEndian
	} else {
		log.Fatalf("Unrecognized endianness flag, %d (must be 0 or ~0)", flag)
		panic("Impossible")
	}
}
