/*package density interpolates sequences of particle positions onto a density
grid.
*/
package sheet

import (
	"math"

	"github.com/phil-mansfield/gotetra/rand"
	"github.com/phil-mansfield/gotetra/geom"
)

// Interpolator creates a grid-based density distribution from seqeunces of
// positions.
type Interpolator interface {
	// Interpolate adds the density distribution implied by points to the
	// density grid used by the Interpolator. Particles should all be within
	// the bounds of the bounding grid and points not within the interpolation
	// grid will be ignored.
	Interpolate(g *Grid, mass float64, xs []geom.Vec)
}

type Grid struct {
	Rhos []float64
	BoxWidth, CellWidth, CellVolume float64
	G, BG geom.Grid
}

type Cell struct {
	Width, X, Y, Z int
}

type cic struct { }
type ngp struct { }

type mcarlo struct {
	subIntr Interpolator
	segWidth int64
	steps int

	gen *rand.Generator

	// Buffers
	idxBuf geom.TetraIdxs
	tet geom.Tetra
	randBuf []float64
	vecBuf []geom.Vec
}

func NewGrid(boxWidth float64, gridWidth int, rhos []float64, c *Cell) *Grid {
	g := &Grid{}
	g.Init(boxWidth, gridWidth, rhos, c)
	return g
}

func (g *Grid) Init(boxWidth float64, gridWidth int, rhos []float64, c *Cell) {
	if len(rhos) != c.Width * c.Width * c.Width {
		panic("Length of rhos doesn't match cell width.")
	}

	g.G.Init(&[3]int{c.X*c.Width, c.Y*c.Width, c.Z*c.Width}, c.Width)
	g.BG.Init(&[3]int{0, 0, 0}, c.Width * gridWidth)
	g.BoxWidth = boxWidth
	g.CellWidth = boxWidth / float64(c.Width * gridWidth)
	g.CellVolume = g.CellWidth * g.CellWidth * g.CellWidth
	g.Rhos = rhos
}

func CloudInCell() Interpolator {
	return &cic{}
}

func NearestGridPoint() Interpolator {
	return &ngp{}
}

func MonteCarlo(segWidth int64, gen *rand.Generator, points int) Interpolator {
	return &mcarlo{
		NearestGridPoint(), segWidth, points, gen,
		geom.TetraIdxs{}, geom.Tetra{},
		make([]float64, points * 3),
		make([]geom.Vec, points),
	}
}

// Interpolate interpolates a sequence of particles onto a density grid via a
// nearest grid point scheme.
func (intr *ngp) Interpolate(g *Grid, mass float64, xs []geom.Vec) {
	frac := mass / g.CellVolume
	cw := float32(g.CellWidth)
	for _, pt := range xs {
		i := int(pt[0] / cw)
		j := int(pt[1] / cw)
		k := int(pt[2] / cw)

		if idx, ok := g.G.IdxCheck(i, j, k); ok {
			g.Rhos[idx] += frac
			continue
		}
	}
}

// Interpolate interpolates a sequence of particles onto a density grid via a
// cloud in cell scheme.
func (intr *cic) Interpolate(g *Grid, mass float64, xs []geom.Vec) {
	frac := mass / g.CellVolume
	cw, cw2 := g.CellWidth, g.CellWidth / 2
	for _, pt := range xs {
		
		xp, yp, zp := float64(pt[0])-cw2, float64(pt[1])-cw2, float64(pt[2])-cw2
		xc, yc, zc := cellPoints(xp, yp, zp, g.CellWidth)
		dx, dy, dz := (xp / cw)-xc, (yp / cw)-yc, (zp / cw)-zc
		tx, ty, tz := 1-dx, 1-dy, 1-dz

		i0, i1 := g.nbrs(int(xc))
		j0, j1 := g.nbrs(int(yc))
		k0, k1 := g.nbrs(int(zc))

		over000 := tx*ty*tz*frac
		over100 := dx*ty*tz*frac
		over010 := tx*dy*tz*frac
		over110 := dx*dy*tz*frac
		over001 := tx*ty*dz*frac
		over101 := dx*ty*dz*frac
		over011 := tx*dy*dz*frac
		over111 := dx*dy*dz*frac

		g.incr(i0, j0, k0, over000)
		g.incr(i1, j0, k0, over100)
		g.incr(i0, j1, k0, over010)
		g.incr(i1, j1, k0, over110)
		g.incr(i0, j0, k1, over001)
		g.incr(i1, j0, k1, over101)
		g.incr(i0, j1, k1, over011)
		g.incr(i1, j1, k1, over111)
	}
}

func (g *Grid) nbrs(i int) (i0, i1 int) {
	if i == -1 {
		return g.BG.Width - 1, 0
	}
	if i+1 == g.BG.Width {
		return i, 0
	}
	return i, i + 1
}

func (g *Grid) incr(i, j, k int, frac float64) {
	if idx, ok := g.G.IdxCheck(i, j, k); ok {
		g.Rhos[idx] += frac
	}
}

func cellPoints(x, y, z, cw float64) (xc, yc, zc float64) {
	return math.Floor(x / cw), math.Floor(y / cw), math.Floor(z / cw)
}

func (intr *mcarlo) Interpolate(g *Grid, mass float64, xs []geom.Vec) {
	cb := &geom.CellBounds{}

	segWidth := intr.segWidth
	gridWidth := segWidth + 1
	pts := len(intr.vecBuf)

	for z := int64(0); z < segWidth; z++ {
		for y := int64(0); y < segWidth; y++ {
			for x := int64(0); x < segWidth; x++ {

				idx := x + y * gridWidth + z * gridWidth * gridWidth

				for dir := 0; dir < 6; dir++ {
					intr.idxBuf.Init(idx, gridWidth, dir)
					
					intr.tet.Init(
						&xs[intr.idxBuf[0]],
						&xs[intr.idxBuf[1]],
						&xs[intr.idxBuf[2]],
						&xs[intr.idxBuf[3]],
						g.BoxWidth,
					)
					
					intr.tet.CellBoundsAt(g.CellWidth, cb)
					ptMass := mass / float64(pts) / 6.0
					intr.tet.RandomSample(intr.gen, intr.randBuf, intr.vecBuf)
					intr.subIntr.Interpolate(g, ptMass, intr.vecBuf)
				}

			}
		}
	}
}
