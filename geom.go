package gotetra

func (h *Header) CompressCoords(x, y, z, dx, dy, dz int64) int64 {
	// Can be spend up with & tricks or conditionals.
	newX := (x + dx + h.CountWidth) % h.CountWidth
	newY := (y + dy + h.CountWidth) % h.CountWidth
	newZ := (z + dz + h.CountWidth) % h.CountWidth

	return newX + newY * h.CountWidth + newZ * h.CountWidth * h.CountWidth
}

func (h *Header) TetraCorners(idx int64) (int64, int64, int64) {
	countArea := h.CountWidth * h.CountWidth

	x := idx % h.CountWidth
	y := (idx % countArea) / h.CountWidth
	z := idx / countArea

	c1 := h.CompressCoords(x, y, z, 0, 0, 1)
	c2 := h.CompressCoords(x, y, z, 0, 1, 0)
	c3 := h.CompressCoords(x, y, z, 1, 0, 0)

	return c1, c2, c3
}
