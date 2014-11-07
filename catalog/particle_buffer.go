package catalog

import (
	tetra "github.com/phil-mansfield/gotetra"
)

// ParticleBuffer is a wrapper around catalog files which allows floats to be
// appended to files on the fly without much overhead from I/O or without
// excessive memory usage/reallocating.
type ParticleBuffer struct {
	buf []tetra.Particle
	idx int
	path string
}

// NewParticleBuffer creates a ParticleBuffer associated with the given file.
func NewParticleBuffer(path string, bufSize int) *ParticleBuffer {
	pb := &ParticleBuffer{ make([]tetra.Particle, bufSize), 0, path }
	return pb
}

// Append adds a value to the float buffer, which will eventually be
// written to the target file.
func (pb *ParticleBuffer) Append(p tetra.Particle) {
	pb.buf[pb.idx] = p
	pb.idx++
	if pb.idx == len(pb.buf) { pb.Flush() }
}

// Flush writes the contents of the buffer to its target file. This will
// be called automatically whenever the buffer fills.
func (pb *ParticleBuffer) Flush() {
	Append(pb.path, pb.buf[0:pb.idx])
	pb.idx = 0
}
