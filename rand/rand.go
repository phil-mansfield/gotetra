package rand

import (
	"math"
	"time"
)

const (
	DefaultBufSize = 1<<10
)

type generatorBackend interface {
	Init(seed uint64)
	Next() float64
	NextSequence(target []float64)
}

type Generator struct {
	backend generatorBackend
	savedGaussian bool
	nextGaussianDx float64
}

type GeneratorType uint8
const (
	Xorshift GeneratorType = iota
	Golang
	Tausworthe

	Default = Tausworthe
)

func NewTimeSeed(gt GeneratorType) *Generator {
	return New(gt, uint64(time.Now().UnixNano()))
}

func New(gt GeneratorType, seed uint64) *Generator {
	var backend generatorBackend

	switch(gt) {
	case Xorshift:
		backend = new(xorshiftGenerator)
	case Golang:
		backend = new(golangGenerator)
	case Tausworthe:
		backend = new(tauswortheGenerator)
	default:
		panic("Unrecognized GeneratorType")
	}

	backend.Init(seed)
	gen := &Generator{ backend, false, -1 }
	return gen
}

func (gen *Generator) UniformInt(low, high int) int {
	f := gen.backend.Next()
	return int(math.Floor(float64(high - low) * f + float64(low)))
	
}

func (gen *Generator) Uniform(low, high float64) float64 {
	if low == 0.0 && high == 1.0 { return gen.backend.Next() }
	return (gen.backend.Next() * (high - low)) + low
}

func (gen *Generator) UniformAt(low, high float64, target []float64) {
	gen.backend.NextSequence(target)
	if low == 0.0 && high == 1.0 { return }
	for i := 0; i < len(target); i ++ {
		target[i] = target[i] * (high - low) + low
	}
}
