package catalog

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"

	"github.com/phil-mansfield/gotetra/geom"
)

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

// Header describes meta-information about the current catalog.
type Header struct {
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
func (gh *gadgetHeader) Standardize() *Header {
	h := &Header{}

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
func ReadGadgetHeader(path string, order binary.ByteOrder) *Header {
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
func ReadGadget(path string, order binary.ByteOrder) (*Header, []Particle) {
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
