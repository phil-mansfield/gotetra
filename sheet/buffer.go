package sheet

import (
	"github.com/phil-mansfield/gotetra/catalog"
	"github.com/phil-mansfield/gotetra/geom"
)

// ParticleBuffer is a wrapper around catalog files which allows floats to be
// appended to files on the fly without much overhead from I/O or without
// excessive memory usage/reallocating.
type ParticleBuffer struct {
	buf  []catalog.Particle
	idx  int
	xs, vs []geom.Vec
}

// NewParticleBuffer creates a ParticleBuffer associated with the given file.
func NewParticleBuffer(xs, vs []geom.Vec, bufSize int) *ParticleBuffer {
	pb := &ParticleBuffer{make([]catalog.Particle, bufSize), 0, xs, vs}
	return pb
}

// Append adds a value to the float buffer, which will eventually be
// written to the target file.
func (pb *ParticleBuffer) Append(ps []catalog.Particle) {
	for _, p := range ps {
		pb.buf[pb.idx] = p
		pb.idx++
		if pb.idx == len(pb.buf) {
			pb.Flush()
		}
	}
}

// Flush writes the contents of the buffer to its target file. This will
// be called automatically whenever the buffer fills.
func (pb *ParticleBuffer) Flush() {
	for i := 0; i < pb.idx; i++ {
		pb.xs[pb.buf[i].Id - 1] = pb.buf[i].Xs
		pb.vs[pb.buf[i].Id - 1] = pb.buf[i].Vs
	}
	pb.idx = 0
}
