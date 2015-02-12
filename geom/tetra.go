package geom

import (
	"fmt"
	"math"

	"github.com/phil-mansfield/gotetra/rand"
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

type TetraIdxs [4]int64

type volumeBuffer struct {
	buf1, buf2, buf3 Vec
}

type sampleBuffer struct {
	d, c [4]Vec
}

const (
	eps = 5e-5

	// TetraDirCount is the number of orientations that a tetrahedron can have
	// within a cube. Should be iterated over when Calling TetraIdxs.Init().
	TetraDirCount = 6
)

var (
	// Yeah, this one was fun to figure out.
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
func NewTetraIdxs(idx, countWidth, skip int64, dir int) *TetraIdxs {
	idxs := &TetraIdxs{}
	return idxs.Init(idx, countWidth, skip, dir)
}

// This wastes an integer multiplication. Oh no!
func compressCoords(x, y, z, dx, dy, dz, countWidth int64) int64 {
	newX := x + dx
	newY := y + dy
	newZ := z + dz

	if newX >= countWidth {
		newX -= countWidth
	} else if newX < 0 {
		newX += countWidth
	}

	if newY >= countWidth {
		newY -= countWidth
	} else if newY < 0 {
		newY += countWidth
	}

	if newZ >= countWidth {
		newZ -= countWidth
	} else if newZ < 0 {
		newZ += countWidth
	}

	return newX + newY*countWidth + newZ*countWidth*countWidth
}

// Init initializes a TetraIdxs collection using the same rules as NewTetraIdxs.
func (idxs *TetraIdxs) Init(idx, countWidth, skip int64, dir int) *TetraIdxs {
	if dir < 0 || dir >= 6 {
		panic(fmt.Sprintf("Unknown direction %d", dir))
	}

	countArea := countWidth * countWidth

	x := idx % countWidth
	y := (idx % countArea) / countWidth
	z := idx / countArea

	idxs[0] = compressCoords(
		x, y, z,
		skip * dirs[dir][0][0], skip * dirs[dir][0][1], skip * dirs[dir][0][2],
		countWidth,
	)
	idxs[1] = compressCoords(
		x, y, z,
		skip * dirs[dir][1][0], skip * dirs[dir][1][1], skip * dirs[dir][1][2],
		countWidth,
	)
	idxs[2] = compressCoords(x, y, z, skip, skip, skip, countWidth)
	idxs[3] = idx

	return idxs
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
	if volSum > vol*(1+eps) {
		return false
	}

	vi = t.signedVolume(v, &t.Corners[1], &t.Corners[3], &t.Corners[2])
	if math.Signbit(vi) != sign {
		return false
	}
	volSum += math.Abs(vi)
	if volSum > vol*(1+eps) {
		return false
	}

	vi = t.signedVolume(v, &t.Corners[0], &t.Corners[3], &t.Corners[1])
	if math.Signbit(vi) != sign {
		return false
	}
	volSum += math.Abs(vi)
	if volSum > vol*(1+eps) {
		return false
	}

	vi = t.signedVolume(v, &t.Corners[0], &t.Corners[2], &t.Corners[3])
	if math.Signbit(vi) != sign {
		return false
	}

	// This last check is neccessary due to periodic boundaries.
	return epsEq(math.Abs(vi)+volSum, vol, eps)
}

func epsEq(x, y, eps float64) bool {
	return (x == 0 && y == 0) ||
		(x != 0 && y != 0 && math.Abs((x-y)/x) <= eps)
}

func (t *Tetra) signedVolume(c1, c2, c3, c4 *Vec) float64 {
	c2.SubAt(c1, t.width, &t.vb.buf1)
	c3.SubAt(c1, t.width, &t.vb.buf2)
	c4.SubAt(c1, t.width, &t.vb.buf3)

	t.vb.buf2.CrossSelf(&t.vb.buf3)
	return t.vb.buf1.Dot(&t.vb.buf2) / 6.0
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

// RandomSample fills a buffer of vecotrs with points generated uniformly at
// random from within a tetrahedron. The length of randBuf must be three times
// the length of vecBuf.
func (tet *Tetra) RandomSample(gen *rand.Generator, randBuf []float64, vecBuf []Vec) {
	N := len(vecBuf)
	if len(randBuf) != N*3 {
		panic(fmt.Sprintf("buf len %d not long enough for %d points.",
			len(randBuf), N))
	}

	gen.UniformAt(0.0, 1.0, randBuf)

	xs := randBuf[0: N]
	ys := randBuf[N: 2*N]
	zs := randBuf[2*N: 3*N]
	tet.Distribute(xs, ys, zs, vecBuf)
}

// Distribute converts a sequences of points generated uniformly within a 
// unit cube to be distributed uniformly within the base tetrahedron. The
// results are placed in vecBuf.
func (tet *Tetra) Distribute(xs, ys, zs []float64, vecBuf []Vec) {
	bary := tet.Barycenter()
	w := float32(tet.width)
	// Some gross code to prevent allocations. cs are the displacement vectors
	// to the corners and the ds are the barycentric components of the random
	// points.
	for i := 0; i < 4; i++ {
		tet.Corners[i].SubAt(bary, tet.width, &tet.sb.c[i])
		tet.sb.d[i].ScaleSelf(0.0)
	}

	// Note: this inner loop is very optimized. Don't try to "fix" it.
	for i := range vecBuf {
		// Find three of the four barycentric coordinates, see
		// C. Rocchini, P. Cignoni, 2001.
		s, t, u := float32(xs[i]), float32(ys[i]), float32(zs[i])

		if s+t > 1 {
			s, t = 1-s, 1-t
		}

		if t+u > 1 {
			t, u = 1-u, 1-s-t
		} else if s+t+u > 1 {
			s, u = 1-t-u, s+t+u-1
		}
		v := 1 - s - t - u

		// Could break loop here, but that flushes the cache and
		// registers.

		for j := 0; j < 3; j++ {
			d0 := tet.sb.c[0][j] * s
			d1 := tet.sb.c[1][j] * t
			d2 := tet.sb.c[2][j] * u
			d3 := tet.sb.c[3][j] * v
			
			val := bary[j] + d0 + d1 + d2 + d3
			if val >= w {
				val -= w
			} else if val < 0 {
				val += w
			}

			vecBuf[i][j] = val
		}
	}
}

// DistributeUnit distributes a set of points in a unit cube across a unit
// tetrahedron and stores the results to vecBuf.
func DistributeUnit(xs, ys, zs []float64, vecBuf []Vec) {
	for i := range vecBuf {
		s, t, u := float32(xs[i]), float32(ys[i]), float32(zs[i])
	
		if s+t > 1 {
			s, t = 1-s, 1-t
		}
	
		if t+u > 1 {
			t, u = 1-u, 1-s-t
		} else if s+t+u > 1 {
			s, u = 1-t-u, s+t+u-1
		}
		
		vecBuf[i][0], vecBuf[i][1], vecBuf[i][2] = s, t, u
	}
}

// DistributeTetra takes a set of points distributed across a unit tetrahedron
// and distributed them across the given tetrahedron through barycentric
// coordinate transformations.
func (tet *Tetra) DistributeTetra(pts []Vec, out []Vec) {
	// TODO: Now that we don't need to worry about boundaries, this can probably
	// be cleaned up.
	bary := tet.Barycenter()

	// Some gross code to prevent allocations. cs are the displacement vectors
	// to the corners and the ds are the barycentric components of the random
	// points.
	for i := 0; i < 4; i++ {
		tet.Corners[i].SubAt(bary, tet.width, &tet.sb.c[i])
		tet.sb.d[i].ScaleSelf(0.0)
	}
	for i := range pts {
		pt := &pts[i]
		s, t, u := pt[0], pt[1], pt[2]
		v := 1 - s - t - u

		for j := 0; j < 3; j++ {
			d0 := tet.sb.c[0][j] * s
			d1 := tet.sb.c[1][j] * t
			d2 := tet.sb.c[2][j] * u
			d3 := tet.sb.c[3][j] * v
			
			val := bary[j] + d0 + d1 + d2 + d3
			out[i][j] = val
		}
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

func (t *Tetra) CellBounds(cellWidth float64) *CellBounds {
	cb := &CellBounds{}
	t.CellBoundsAt(cellWidth, cb)
	return cb
}

func (t *Tetra) CellBoundsAt(cellWidth float64, cb *CellBounds) {
	mins := &t.sb.d[0]
	maxes := &t.sb.d[1]
	for i := 0; i < 3; i++ {
		mins[i] = t.Corners[0][i]
		maxes[i] = t.Corners[0][i]
	}

	for j := 1; j < 4; j++ {
		for i := 0; i < 3; i++ {
			if t.Corners[j][i] < mins[i] {
				mins[i] = t.Corners[j][i]
			} else if t.Corners[j][i] > maxes[i] {
				maxes[i] = t.Corners[j][i]
			}
		}
	}

	for i := 0; i < 3; i++ {
		maxes[i] -= mins[i]
	}
	width := maxes
	origin := mins


	for i := 0; i < 3; i++ {
		cb.Origin[i] = int(math.Floor(float64(origin[i]) / cellWidth))
		cb.Width[i] = 1 + int(math.Floor(
			float64(width[i] + origin[i]) / cellWidth),
		)
		cb.Width[i] -= cb.Origin[i]
	}
}
