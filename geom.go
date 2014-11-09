package gotetra

import (
	"math"
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
	d2 := low + h.TotalWidth

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
