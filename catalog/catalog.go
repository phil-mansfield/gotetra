package catalog

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	
	"unsafe"

	"github.com/phil-mansfield/gotetra"
)

const (
	DefaultEndiannessFlag = -1
)

type GadgetHeader struct {
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

// ReadInt32 returns single 32-bit interger from the given file generated on
// a machine with the given endianness.
func ReadInt32(r io.Reader, order binary.ByteOrder) int32 {
	var n int32
	if err := binary.Read(r, order, &n); err != nil { panic(err) }
	return n
}

func WriteInt32(w io.Writer, order binary.ByteOrder, n int32) {
	if err := binary.Write(w, order, n); err != nil { panic(err) }
}

// Standardize returns a standardized header the corresponds to the source
// Gadget header.
func (gh *GadgetHeader) Standardize() *gotetra.Header {
	h := &gotetra.Header{}

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

// WrapDistance returns a value of x which is inside the box described by h.
func (h *GadgetHeader) WrapDistance(x float64) float64 {
	if x < 0 {
		return x + h.BoxSize
	} else if x >= h.BoxSize {
		return x - h.BoxSize
	}
	return x
}

// ReadGadget reads the gadget particle catalog located at the given location
// generated on a machine with the given endianness and returns the particles
// along with the standardized header file contained within.
func ReadGadget(fileName string, order binary.ByteOrder) (*gotetra.Header, []gotetra.Particle) {
	f, err := os.Open(fileName)
	if err != nil { panic(err) }
	defer f.Close()

	gh := &GadgetHeader{}

	_ = ReadInt32(f, order)
	binary.Read(f, binary.LittleEndian, gh)
	_ = ReadInt32(f, order)
	
	h := gh.Standardize()
	floatBuf := make([]float32, 3 * h.Count)
	ps := make([]gotetra.Particle, h.Count)

	_ = ReadInt32(f, order)
	binary.Read(f, order, floatBuf)
	_ = ReadInt32(f, order)

	for i := range ps { 
		ps[i].Xs[0] = gh.WrapDistance(float64(floatBuf[3 * i + 0]))
		ps[i].Xs[1] = gh.WrapDistance(float64(floatBuf[3 * i + 1]))
		ps[i].Xs[2] = gh.WrapDistance(float64(floatBuf[3 * i + 2]))
	}

	_ = ReadInt32(f, order)
	binary.Read(f, order, floatBuf)
	_ = ReadInt32(f, order)

	rootA := float32(math.Sqrt(float64(gh.Time)))
	for i := range ps { 
		ps[i].Vs[0] = float64(floatBuf[3 * i + 0] * rootA)
		ps[i].Vs[1] = float64(floatBuf[3 * i + 1] * rootA)
		ps[i].Vs[2] = float64(floatBuf[3 * i + 2] * rootA)
	}

	ids := make([]int64, h.Count)

	_ = ReadInt32(f, order)
	binary.Read(f, order, ids)
	_ = ReadInt32(f, order)

	for i := range ps { ps[i].Id = ids[i] }

	return h, ps
}

// Write writes a standardized header and particle slice to the given location
func Write(path string, h *gotetra.Header, ps []gotetra.Particle) {
	f, err := os.Create(path)
	if err != nil { panic(err) }
	defer f.Close()

	order := endianness(DefaultEndiannessFlag)

	err = binary.Write(f, order, DefaultEndiannessFlag)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, int32(unsafe.Sizeof(gotetra.Header{})))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, int32(unsafe.Sizeof(ps[0])))
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, h)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, ps)
	if err != nil { panic(err.Error()) }
}

func readHeader(path string) (*gotetra.Header, *os.File, binary.ByteOrder) {
	f, err := os.Open(path)
	if err != nil { panic(err) }

	order := endianness(ReadInt32(f, binary.LittleEndian))

	// Sanity checks:
	headerSize := ReadInt32(f, order)
	particleSize := ReadInt32(f, order)
	if int32(unsafe.Sizeof(gotetra.Header{})) != headerSize {
		panic(fmt.Sprintf(
			"Size of header in code, %d, does not match size in catalog, %d.",
			unsafe.Sizeof(gotetra.Header{}), headerSize,
		))
	} else if int32(unsafe.Sizeof(gotetra.Particle{})) != particleSize {
		panic(fmt.Sprintf(
			"Size of header in code, %d, does not match size in catalog, %d.",
			unsafe.Sizeof(gotetra.Particle{}), particleSize,
		))
	}

	h := &gotetra.Header{}
	err = binary.Read(f, order, h)
	if err != nil { panic(err.Error()) }

	return h, f, order
}

// Append appends a particle list to the given file containing a standardized
// header and updates the header accordingly.
func Append(path string, ps []gotetra.Particle) {
	h, f, order := readHeader(path)
	defer f.Close()

	_, err := f.Seek(int64(unsafe.Sizeof(ps[0])) * h.Count, 1)
	if err != nil { panic(err.Error()) }
	err = binary.Write(f, order, ps)
	if err != nil { panic(err.Error()) }
	_, err = f.Seek(12, 0)
	if err != nil { panic(err.Error()) }

	h.Count += int64(len(ps))
	err = binary.Write(f, order, h)
	if err != nil { panic(err.Error()) }
}

// Read reads a header and particle list from a standardized catalog header.
func Read(path string) (*gotetra.Header, []gotetra.Particle) {
	h, f, order := readHeader(path)
	defer f.Close()

	ps := make([]gotetra.Particle, h.Count)
	err := binary.Read(f, order, ps)
	if err != nil { panic(err.Error()) }

	return h, ps
}

func endianness(flag int32) binary.ByteOrder {
	if flag == 0 {
		return binary.LittleEndian
	} else if flag == -1 {
		return binary.BigEndian
	} else {
		panic("Unrecognized endianness flag.")
	}
}

func endiannessFlag(order binary.ByteOrder) int32 {
	if order == binary.LittleEndian {
		return 0
	} else {
		return -1
	}
}
