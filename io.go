package gotetra

import (
	"encoding/binary"
	"os"
	"unsafe"
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
}

func ReadInt32(f *os.File, order binary.ByteOrder) int32 {
	var n int32
	if err := binary.Read(f, order, &n); err != nil { panic(err) }
	return n
}

func ReadGadget(fileName string, order binary.ByteOrder) ([]Particle, *Header) {
	f, err := os.Open(fileName)
	if err != nil { panic(err) }

	headerSize := ReadInt32(f, order)	
	header := &GadgetHeader{}

	println(headerSize)
	println(unsafe.Sizeof(*header))

	panic("Whyyyy")
}
