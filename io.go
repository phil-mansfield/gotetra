package gotetra

import (
	"encoding/binary"
	"math"
	"os"
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

func ReadInt32(f *os.File, order binary.ByteOrder) int32 {
	var n int32
	if err := binary.Read(f, order, &n); err != nil { panic(err) }
	return n
}

func (gh *GadgetHeader) Standardize() *Header {
	h := &Header{}

	h.Count = int64(gh.NPart[1] + gh.NPart[0] << 32)
	h.TotalCount = int64(gh.NPartTotal[1] + gh.NPartTotal[0] << 32)
	h.Mass = float64(gh.Mass[1])
	h.BoxSize = float64(gh.BoxSize)

	h.Cosmo.Z = gh.Redshift
	h.Cosmo.OmegaM = gh.Omega0
	h.Cosmo.OmegaL = gh.OmegaLambda
	h.Cosmo.H100 = gh.HubbleParam

	return h
}

func (h *GadgetHeader) WrapDistance(x float64) float64 {
	if x < 0 {
		return x + h.BoxSize
	} else if x >= h.BoxSize {
		return x - h.BoxSize
	}
	return x
}

func ReadGadget(fileName string, order binary.ByteOrder) ([]Particle, *Header) {
	f, err := os.Open(fileName)
	if err != nil { panic(err) }

	gh := &GadgetHeader{}

	_ = ReadInt32(f, order)
	binary.Read(f, order, gh)
	_ = ReadInt32(f, order)

	h := gh.Standardize()
	floatBuf := make([]float32, 3 * h.Count)
	ps := make([]Particle, h.Count)

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

	return ps, h
}

func ReadOnePositionGadget(fileName string, order binary.ByteOrder) (*Header, []float64) {
	f, err := os.Open(fileName)
	if err != nil { panic(err) }

	gh := &GadgetHeader{}
	_ = ReadInt32(f, order)
	binary.Read(f, order, gh)
	_ = ReadInt32(f, order)	

	h := gh.Standardize()

	pos := make([]float64, 3)
	_ = ReadInt32(f, order)
	binary.Read(f, order, pos)
	
	pos[0] = gh.WrapDistance(float64(pos[0]))
	pos[1] = gh.WrapDistance(float64(pos[1]))
	pos[2] = gh.WrapDistance(float64(pos[2]))

	return h, pos
}
