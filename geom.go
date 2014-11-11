package gotetra

import (
	"math"
)

const (
	eps = 1e-6
)

func (h *Header) compressCoords(x, y, z, dx, dy, dz int64) int64 {
	// Can be send up with & tricks or conditionals.
	newX := (x + dx + h.CountWidth) % h.CountWidth
	newY := (y + dy + h.CountWidth) % h.CountWidth
	newZ := (z + dz + h.CountWidth) % h.CountWidth

	return newX + newY * h.CountWidth + newZ * h.CountWidth * h.CountWidth
}

func (h *Header) TetraCorners(idx int64, out []int64) {
	countArea := h.CountWidth * h.CountWidth

	x := idx % h.CountWidth
	y := (idx % countArea) / h.CountWidth
	z := idx / countArea

	out[0] = h.compressCoords(x, y, z, 0, 0, 1)
	out[1] = h.compressCoords(x, y, z, 0, 1, 0)
	out[2] = h.compressCoords(x, y, z, 1, 0, 0)
}

func (h *Header) wrapDist(x1, x2 float64) float64 {
	var low, high float64

	if x1 < x2 {
		low, high = x1, x2
	} else {
		low, high = x2, x1
	}

	d1 := high - low
	d2 := low + h.TotalWidth - high

	if d1 > d2 {
		return d2
	} else {
		return d1
	}
}

func (h *Header) Distance(p1, p2 *Particle) float64 {
	dx := h.wrapDist(float64(p1.Xs[0]), float64(p2.Xs[0]))
	dy := h.wrapDist(float64(p1.Xs[1]), float64(p2.Xs[1]))
	dz := h.wrapDist(float64(p1.Xs[2]), float64(p2.Xs[2]))

	return math.Sqrt(dx * dx + dy * dy + dz * dz)
}

type VolumeBuffer struct {
	buf1, buf2, buf3, bufCross [3]float64
}

func NewVolumeBuffer() *VolumeBuffer {
	return &VolumeBuffer{ }
}

func (h *Header) Volume(ps []*Particle, vb *VolumeBuffer) float64 {
	h.subX(ps[1], ps[0], &vb.buf1)
	h.subX(ps[2], ps[0], &vb.buf2)
	h.subX(ps[3], ps[0], &vb.buf3)

	cross(&vb.buf2, &vb.buf3, &vb.bufCross)
	return math.Abs(dot(&vb.buf1, &vb.bufCross)) / 6.0
}

func (h *Header) WithinTetra(
	p *Particle,
	ps []*Particle,
	vol float64,
	vb *VolumeBuffer,
) bool {

	buf := []*Particle{ps[0], ps[1], ps[2], p}
	sum := 0.0
	for i := 0; i < 4; i++ {
		sum += h.Volume(buf, vb)
		buf[i] = ps[(i + 3) % 4]
	}

	return (math.Abs(sum - vol) / vol) <= eps
}

func (h *Header) subX(p1, p2 *Particle, out *[3]float64) {
	for i := 0; i < 3; i++ {
		(*out)[i] = float64(p1.Xs[i]) - float64(p2.Xs[i])
	}
}

func cross(v1, v2, out *[3]float64) {
	(*out)[0] = v1[1] * v2[2] - v1[2] * v2[1]
	(*out)[1] = v1[2] * v2[0] - v1[0] * v2[2]
	(*out)[2] = v1[0] * v2[1] - v1[1] * v2[0]
}

func dot(v1, v2 *[3]float64) float64 {
	sum := 0.0
	for i := 0; i < 3; i++ {
		sum += (*v1)[i] * (*v2)[i] 
	}
	return sum
}
