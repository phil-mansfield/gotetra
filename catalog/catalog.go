/*package catalog contains funcitons for reading from and writing to
files containing catalogs of particle locations and velocities.

Catalogs written by external programs are standardized to a single type when
read. Currently the only such extrenal catalog type that is supported the
Gadget 2 particle catalog.

The binary format used is as follows:
    |-- 1 --||-- 2 --||-- 3 --||-- ... 4 ... --||-- ... 5 ... --|
    
    1 - (int32) Flag indicating the endianness of the file. 0 indicates a big
        endian byte ordering and -1 indicates a little endian byte order.
    2 - (int32) Size of a Header struct. Should be checked for consistency.
    3 - (int32) Size of a Particle struct. Should be checked for consistency.
    4 - (tetra.Header) Header file contiainung meta-information about the
        particle catalog.
    5 - ([]tetra.Particle) Contiguous block of particles. Garuanteed to be
        of size unsafe.Sizeof(Particle{}) * header.Count.
*/
package catalog

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	
	"unsafe"

	tetra "github.com/phil-mansfield/gotetra"
)

const (
	// Endianness used by default when writing catalogs. Catalogs of any
	// endianness can be read.
	DefaultEndiannessFlag int32 = -1
)

// gadgetHeader is the formatting for meta-information used by Gadget 2.
type gadgetHeader struct {
	NPart [6]uint32
	Mass [6]float64
	Time, Redshift float64
	FlagSfr, FlagFeedback int32
	NPartTotal [6]uint32
	FlagCooling, NumFiles int32
	BoxSize, Omega0, OmegaLambda, HubbleParam float64
	FlagStellarAge, HashTabSize int32

	Padding [88]byte
}

// readInt32 returns single 32-bit interger from the given file using the
// given endianness.
func readInt32(r io.Reader, order binary.ByteOrder) int32 {
	var n int32
	if err := binary.Read(r, order, &n); err != nil { panic(err) }
	return n
}

// Standardize returns a tetra.Header that corresponds to the source
// Gadget 2 header.
func (gh *gadgetHeader) Standardize() *tetra.Header {
	h := &tetra.Header{}

	h.Count = int64(gh.NPart[1] + gh.NPart[0] << 32)
	h.TotalCount = int64(gh.NPartTotal[1] + gh.NPartTotal[0] << 32)
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

// ReadGadget reads the gadget particle catalog located at the given location
// and written with the given endianness. Its header and particle sequence
// are returned in a standardized format.
func ReadGadget(path string, order binary.ByteOrder) (*tetra.Header, []tetra.Particle) {
	f, err := os.Open(path)
	if err != nil { panic(err) }
	defer f.Close()

	gh := &gadgetHeader{}

	_ = readInt32(f, order)
	binary.Read(f, binary.LittleEndian, gh)
	_ = readInt32(f, order)
	
	h := gh.Standardize()
	floatBuf := make([]float32, 3 * h.Count)
	ps := make([]tetra.Particle, h.Count)

	_ = readInt32(f, order)
	binary.Read(f, order, floatBuf)
	_ = readInt32(f, order)

	for i := range ps { 
		ps[i].Xs[0] = gh.WrapDistance(float64(floatBuf[3 * i + 0]))
		ps[i].Xs[1] = gh.WrapDistance(float64(floatBuf[3 * i + 1]))
		ps[i].Xs[2] = gh.WrapDistance(float64(floatBuf[3 * i + 2]))
	}

	_ = readInt32(f, order)
	binary.Read(f, order, floatBuf)
	_ = readInt32(f, order)

	rootA := float32(math.Sqrt(float64(gh.Time)))
	for i := range ps { 
		ps[i].Vs[0] = float64(floatBuf[3 * i + 0] * rootA)
		ps[i].Vs[1] = float64(floatBuf[3 * i + 1] * rootA)
		ps[i].Vs[2] = float64(floatBuf[3 * i + 2] * rootA)
	}

	ids := make([]int64, h.Count)

	_ = readInt32(f, order)
	binary.Read(f, order, ids)
	_ = readInt32(f, order)

	for i := range ps { ps[i].Id = ids[i] }

	return h, ps
}

// Write writes the given header and particle sequence to the specified file.
func Write(path string, h *tetra.Header, ps []tetra.Particle) {
	f, err := os.Create(path)
	if err != nil { panic(err) }
	defer f.Close()

	h.Count = int64(len(ps))
	order := endianness(DefaultEndiannessFlag)

	err = binary.Write(f, order, DefaultEndiannessFlag)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, int32(unsafe.Sizeof(tetra.Header{})))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, int32(unsafe.Sizeof(ps[0])))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, h)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, ps)
	if err != nil { panic(err.Error()) }
}

// readHeader reads only the header from the given file.
func readHeader(path string) (*tetra.Header, *os.File, binary.ByteOrder) {
	f, err := os.OpenFile(path, os.O_RDWR, os.ModePerm)
	if err != nil { panic(err) }

	order := endianness(readInt32(f, binary.LittleEndian))

	// Sanity checks:
	headerSize := readInt32(f, order)
	particleSize := readInt32(f, order)
	if int32(unsafe.Sizeof(tetra.Header{})) != headerSize {
		panic(fmt.Sprintf(
			"Size of header in code, %d, does not match size in catalog, %d.",
			unsafe.Sizeof(tetra.Header{}), headerSize,
		))
	} else if int32(unsafe.Sizeof(tetra.Particle{})) != particleSize {
		panic(fmt.Sprintf(
			"Size of header in code, %d, does not match size in catalog, %d.",
			unsafe.Sizeof(tetra.Particle{}), particleSize,
		))
	}

	h := &tetra.Header{}
	err = binary.Read(f, order, h)
	if err != nil { panic(err.Error()) }

	return h, f, order
}

// Append appends a particle sequence to the end of the given file.
func Append(path string, ps []tetra.Particle) {
	h, f, order := readHeader(path)
	defer f.Close()

	idx, err := f.Seek(0, 2)
	if err != nil { panic(err.Error()) }
	println(idx)
	err = binary.Write(f, order, ps)
	if err != nil { panic(err.Error()) }
	_, err = f.Seek(12, 0)
	if err != nil { panic(err.Error()) }

	h.Count += int64(len(ps))
	err = binary.Write(f, order, h)
	if err != nil { panic(err.Error()) }
}

// Read reads a header and particle sequence from the given file.
func Read(path string) (*tetra.Header, []tetra.Particle) {
	h, f, order := readHeader(path)
	defer f.Close()

	ps := make([]tetra.Particle, h.Count)
	err := binary.Read(f, order, ps)
	if err != nil { panic(err.Error()) }

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
