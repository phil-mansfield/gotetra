package gotetra

import (
	"fmt"
	"math"
)

const (
	eps = 1e-6
)

var (
	dirs = [6][4][3]int64{
		{{1, 0, 0}, {1, 1, 0}},
		{{1, 0, 0}, {1, 0, 1}},
		{{0, 1, 0}, {1, 1, 0}},
		{{0, 0, 1}, {1, 0, 1}},
		{{0, 1, 0}, {0, 1, 1}},
		{{0, 0, 1}, {0, 1, 1}},
	}
)

// Note, should probably make this a [4]Particle
type Tetra []*Particle

func (t Tetra) Valid() bool {
	for i := 0; i < 4; i++ {
		if t[i] == nil { return false }
	}
	return true
}

type Bounds struct {
	MinX, MaxX, MinY, MaxY, MinZ, MaxZ int
}

func (h *Header) compressCoords(x, y, z, dx, dy, dz int64) int64 {
	newX := (x + dx + h.CountWidth) % h.CountWidth
	newY := (y + dy + h.CountWidth) % h.CountWidth
	newZ := (z + dz + h.CountWidth) % h.CountWidth

	return newX + newY * h.CountWidth + newZ * h.CountWidth * h.CountWidth
}

func (h *Header) TetraCorners(idx int64, dir int, out []int64) {
	if dir < 0 || dir >= 6 {
		panic(fmt.Sprintf("Unknown direction %d", dir))
	}

	countArea := h.CountWidth * h.CountWidth

	x := idx % h.CountWidth
	y := (idx % countArea) / h.CountWidth
	z := idx / countArea

	out[0] = h.compressCoords(
		x, y, z, dirs[dir][0][0], dirs[dir][0][1], dirs[dir][0][2],
	)
	out[1] = h.compressCoords(
		x, y, z, dirs[dir][1][0], dirs[dir][1][1], dirs[dir][1][2],
	)
	out[2] = h.compressCoords(x, y, z, 1, 1, 1)
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

func (h *Header) Volume(ps Tetra, vb *VolumeBuffer) float64 {
	h.subX(ps[1], ps[0], &vb.buf1)
	h.subX(ps[2], ps[0], &vb.buf2)
	h.subX(ps[3], ps[0], &vb.buf3)

	cross(&vb.buf2, &vb.buf3, &vb.bufCross)
	return math.Abs(dot(&vb.buf1, &vb.bufCross)) / 6.0
}

func (h *Header) WithinTetra(
	p *Particle,
	ps Tetra,
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

func (t Tetra) Bounds(cellWidth float64, ib *Bounds) *Bounds {
	if ib == nil {
		ib = &Bounds{}
	}

	minX := float64(t[0].Xs[0])
	minY := float64(t[0].Xs[1])
	minZ := float64(t[0].Xs[2])
	maxX, maxY, maxZ := minX, minY, minZ

	for i := 1; i < 4; i++ {
		minX, maxX = minMax(float64(t[i].Xs[0]), minX, maxX)
		minY, maxY = minMax(float64(t[i].Xs[1]), minY, maxY)
		minZ, maxZ = minMax(float64(t[i].Xs[2]), minZ, maxZ)
	}
	
	ib.MinX = int(minX / cellWidth)
	ib.MaxX = int(math.Ceil(maxX / cellWidth))
	ib.MinY = int(minY / cellWidth)
	ib.MaxY = int(math.Ceil(maxY / cellWidth))
	ib.MinZ = int(minZ / cellWidth)
	ib.MaxZ = int(math.Ceil(maxZ / cellWidth))

	return ib
}

func (h *Header) minMax(coord int64) (min, max float64) {
	min = float64(coord) * h.Width
	return min, min + h.Width
}

func minMax(x, oldMin, oldMax float64) (min, max float64) {
	if x > oldMax {
		return oldMin, x
	} else if x < oldMin {
		return x, oldMax
	} else {
		return oldMin, oldMax
	}
}

type Point [3]float64

// source:
// http://vcg.isti.cnr.it/activities/geometryegraphics/pointintetraedro.html
func (h *Header) Sample(t Tetra, randBuf []float64, out []Point) {
	if len(randBuf) * 3 != len(out) {
		panic(fmt.Sprintf("buf len %d not long enough for %d points.",
			len(randBuf), len(out)))
	}

	for i := range out {
		// Three of the four barycentric coordinates
		t1, t2, t3 := randBuf[i], randBuf[i + 1], randBuf[i + 2]

		if t1 + t2 + t3 < 1.0 {
			continue
		} else if t2 + t3 > 1.0 {
			t1, t2, t3 = t1, 1.0 - t3, 1.0 - t1 - t2
		} else {
			t1, t2, t3 = 1.0 - t2 - t3, t2,  t1 + t2 + t3 - 1.0
		}

		t4 := 1.0 - t1 - t2 - t3

		h.Scale(t1, &out[i], &out[i])
		h.Scale(t2, &out[i], &out[i])
		h.Scale(t3, &out[i], &out[i])
		h.Scale(t4, &out[i], &out[i])
	}
}

func (h *tetra.Header) Scale(k float64, v, out *Point) {
	for i := 0; i < 3; i ++ { out[i] = v[i] * 3 }
}
