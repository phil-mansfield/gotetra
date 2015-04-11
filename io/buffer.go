package io

import (
	"github.com/phil-mansfield/gotetra/geom"
)

// ParticleBuffer is a wrapper around catalog files which allows floats to be
// appended to files on the fly without much overhead from I/O or without
// excessive memory usage/reallocating.
type ParticleBuffer struct {
	xBuf, vBuf []geom.Vec
	idBuf []int64
	idx int
	xs, vs []geom.Vec
}

// NewParticleBuffer creates a ParticleBuffer associated with the given file.
func NewParticleBuffer(xs, vs []geom.Vec, bufSize int) *ParticleBuffer {
	pb := &ParticleBuffer{
		make([]geom.Vec, bufSize),
		make([]geom.Vec, bufSize),
		make([]int64, bufSize),
		0, xs, vs,
	}
	return pb
}

// Append adds a value to the float buffer, which will eventually be
// written to the target file.
func (pb *ParticleBuffer) Append(xBuf, vBuf []geom.Vec, idBuf []int64) {
	for pi := range xBuf {
		pb.xBuf[pb.idx] = xBuf[pi]
		pb.vBuf[pb.idx] = vBuf[pi]
		pb.idBuf[pb.idx] = idBuf[pi]
		
		pb.idx++
		if pb.idx == len(pb.xBuf) {
			pb.Flush()
		}
	}
}

// Flush writes the contents of the buffer to its target file. This will
// be called automatically whenever the buffer fills.
func (pb *ParticleBuffer) Flush() {
	for i := 0; i < pb.idx; i++ {
		pb.xs[pb.idBuf[i] - 1] = pb.xBuf[i]
		pb.vs[pb.idBuf[i] - 1] = pb.vBuf[i]
	}
	pb.idx = 0
}
