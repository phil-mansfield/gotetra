package geom

import (
	"fmt"
	"math"

	"github.com/phil-mansfield/num/rand"
)

// Tetra is a tetrahedron with points inside a box with periodic boundary
// conditions.
//
// NOTE: To speed up computations, Tetra contains a number of non-trivial
// private fields. If there is a need to store a large number of tetrahedra,
// it is advised to construct a slice of [4]Points instead and switch to
// Tetra instances only for calculations.
type Tetra struct {
	Corners [4]Vec
	volume  float64
	bary    Vec
	width   float64
	vb      volumeBuffer
	sb      sampleBuffer

	volumeValid, baryValid bool
}

type TetraIdxs struct {
	Corners [4]int64
}

type volumeBuffer struct {
	buf1, buf2, buf3 Vec
}

type sampleBuffer struct {
	d, c [4]Vec
}

const (
	eps = 1e-6

	TetraDirCount = 6
)

var (
	dirs = [TetraDirCount][2][3]int64{
		{{1, 0, 0}, {1, 1, 0}},
		{{1, 0, 0}, {1, 0, 1}},
		{{0, 1, 0}, {1, 1, 0}},
		{{0, 0, 1}, {1, 0, 1}},
		{{0, 1, 0}, {0, 1, 1}},
		{{0, 0, 1}, {0, 1, 1}},
	}
)

// NewTetra creates a new tetrahedron with corners at the specified positions
// within a periodic box of the given width.
func NewTetra(c1, c2, c3, c4 *Vec, width float64) (t *Tetra, ok bool) {
	t = &Tetra{}
	return t, t.Init(c1, c2, c3, c4, width)
}

// Init initializes a tetrahedron to correspond to the given corners. It returns
// true if all the given pointers are valid and the tetrahedron was properly
// initialized and false otherwise. This behavior is chosen so that this
// function interacts nicely with ParticleManagers.
func (t *Tetra) Init(c1, c2, c3, c4 *Vec, width float64) (ok bool) {
	t.volumeValid = false
	t.baryValid = false

	if c1 == nil || c2 == nil || c3 == nil || c4 == nil {
		return false
	}

	c1.ModAt(width, &t.Corners[0])
	c2.ModAt(width, &t.Corners[1])
	c3.ModAt(width, &t.Corners[2])
	c4.ModAt(width, &t.Corners[3])

	t.width = width

	// Remaining fields are buffers and need not be initialized.

	return true
}

// NewTetraIdxs creates a collection of indices corresponding to a tetrhedron
// with an anchor point at the given index. The parameter dir selects a
// particular tetrahedron configuration and must line in the range
// [0, TetraDirCount).
func NewTetraIdxs(idx, gridWidth int64, dir int) *TetraIdxs {
	idxs := &TetraIdxs{}
	return idxs.Init(idx, gridWidth, dir)
}

// Init initializes a TetraIdxs collection using the same rules as NewTetraIdxs.
func (idxs *TetraIdxs) Init(idx, gridWidth int64, dir int) *TetraIdxs {
	return nil
}

// Volume computes the volume of a tetrahedron.
func (t *Tetra) Volume() float64 {
	if t.volumeValid {
		return t.volume
	}

	t.volume = math.Abs(t.signedVolume(
		&t.Corners[0], &t.Corners[1], &t.Corners[2], &t.Corners[3]),
	)

	t.volumeValid = true
	return t.volume
}

// Contains returns true if a tetrahedron contains the given point and false
// otherwise.
func (t *Tetra) Contains(v *Vec) bool {
	vol := t.Volume()

	// (my appologies for the gross code here)
	vi := t.signedVolume(v, &t.Corners[0], &t.Corners[1], &t.Corners[2])
	volSum := math.Abs(vi)
	sign := math.Signbit(vi)
	if volSum > vol * (1 + eps) {
		return false
	}

	vi = t.signedVolume(v, &t.Corners[1], &t.Corners[3], &t.Corners[2])
	if math.Signbit(vi) != sign {
		return false
	}
	volSum += math.Abs(vi)
	if volSum > vol * (1 + eps) {
		return false
	}

	vi = t.signedVolume(v, &t.Corners[0], &t.Corners[3], &t.Corners[1])
	if math.Signbit(vi) != sign {
		return false
	}
	volSum += math.Abs(vi)
	if volSum > vol * (1 + eps) {
		return false
	}

	vi = t.signedVolume(v, &t.Corners[0], &t.Corners[2], &t.Corners[3])
	if math.Signbit(vi) != sign {
		return false
	}
	return true
}

func epsEq(x, y, eps float64) bool {
	return (x == 0 && y == 0) || math.Abs((x - y) / x) <= eps
}

// TODO: Think about whether or not this actually does what you want with the
// sign bit.
func (t *Tetra) signedVolume(c1, c2, c3, c4 *Vec) float64 {
	c2.SubAt(c1, t.width, &t.vb.buf1)
	c3.SubAt(c1, t.width, &t.vb.buf2)
	c4.SubAt(c1, t.width, &t.vb.buf3)

	t.vb.buf2.CrossSelf(&t.vb.buf3)

	return t.vb.buf1.Dot(&t.vb.buf2) / 6.0
}

// CellBounds returns the bounding box of a tetrahedron, aligned to the given
// grid. The indices returned represent a tetrahedron whose barycenter is
// within the fundamental domain of the grid. As such, the returned indices
// may not be within the domain [0, g.Width).
func (t *Tetra) CellBounds(g *Grid) *CellBounds {
	return t.CellBoundsAt(g, &CellBounds{})
}

// CellBoundsAt returns the same quantity as CellBounds, but the result is
// placed at the given location.
func (t *Tetra) CellBoundsAt(g *Grid, out *CellBounds) *CellBounds {
	bary := t.Barycenter()

	for i := 0; i < 4; i++ {
		t.Corners[i].SubAt(bary, t.width, &t.sb.c[i])
	}

	var minDs, maxDs [3]float32

	for i := 0; i < 4; i++ {
		for d := 0; d < 3; d++ {
			if d == 0 {
				minDs[d], maxDs[d] = t.sb.c[i][d], t.sb.c[i][d]
			} else {
				minDs[d], maxDs[d] = minMax(t.sb.c[i][d], minDs[d], maxDs[d])
			}
		}
	}

	mult := float64(g.Width) / t.width
	for d := 0; d < 3; d++ {
		fIdx := float64(bary[d] + minDs[d]) * mult
		out.Min[d] =  int(math.Floor(fIdx)) - g.Origin[d]
		fIdx = float64(bary[d] + maxDs[d]) * mult
		out.Max[d] = int(math.Ceil(fIdx)) - g.Origin[d]
	}

	return out
}

func minMax(x, oldMin, oldMax float32) (min, max float32) {
	if x > oldMax {
		return oldMin, x
	} else if x < oldMin {
		return x, oldMax
	} else {
		return oldMin, oldMax
	}
}

// Sample fills a buffer of vecotrs with points generated uniformly at random
// from within a tetrahedron. The length of randBuf must be three times the
// length of vecBuf.
func (t *Tetra) Sample(gen *rand.Generator, randBuf []float64, vecBuf []Vec) {
	if len(randBuf) != len(vecBuf) *3 {
		panic(fmt.Sprintf("buf len %d not long enough for %d points.",
			len(randBuf), len(vecBuf)))
	}

	gen.UniformAt(0.0, 1.0, randBuf)
	bary := t.Barycenter()

	// Some gross code to prevent allocations. cs are the displacement vectors
	// to the corners and the ds are the barycentric components of the random
	// points.
	for i := 0; i < 4; i++ {
		t.Corners[i].SubAt(bary, t.width, &t.sb.c[i])
		t.sb.d[i].ScaleSelf(0.0)
	}

	for i := range vecBuf {
		// Generate three of the four barycentric coordinates
		t1, t2, t3 := randBuf[i*3], randBuf[i*3+1], randBuf[i*3+2]

		if t1+t2+t3 < 1.0 {
			continue
		} else if t2+t3 > 1.0 {
			t1, t2, t3 = t1, 1.0-t3, 1.0-t1-t2
		} else {
			t1, t2, t3 = 1.0-t2-t3, t2, t1+t2+t3-1.0
		}

		// Solve for the last one.
		t4 := 1.0 - t1 - t2 - t3

		t.sb.c[0].ScaleAt(t1, &t.sb.d[0])
		t.sb.c[1].ScaleAt(t2, &t.sb.d[1])
		t.sb.c[2].ScaleAt(t3, &t.sb.d[2])
		t.sb.c[3].ScaleAt(t4, &t.sb.d[3])

		t.sb.c[0].AddAt(&t.sb.c[1], &vecBuf[i])
		vecBuf[i].AddSelf(&t.sb.c[2]).AddSelf(&t.sb.c[3])
	}
}

// Barycenter computes the barycenter of a tetrahedron.
func (t *Tetra) Barycenter() *Vec {
	if t.baryValid {
		return &t.bary
	}

	buf1, buf2 := &t.vb.buf1, &t.vb.buf2
	center(&t.Corners[0], &t.Corners[1], buf1, t.width)
	center(&t.Corners[2], &t.Corners[3], buf2, t.width)
	center(buf1, buf2, &t.bary, t.width)
	t.bary.ModSelf(t.width)
	t.baryValid = true

	return &t.bary
}

func center(r1, r2, out *Vec, width float64) {
	r1.SubAt(r2, width, out).AddSelf(r2).AddSelf(r2).ScaleSelf(0.5)
}
