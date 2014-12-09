/*package density interpolates sequences of particle positions onto a density
grid.
*/
package density

import (
	"log"
	"math"

	"github.com/phil-mansfield/num/rand"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/catalog"
)

// Interpolator creates a grid-based density distribution from seqeunces of
// positions.
type Interpolator interface {
	// Interpolate adds the density distribution implied by points to the
	// density grid used by the Interpolator. Particles should all be within
	// the bounds of the bounding grid and points not within the interpolation
	// grid will be ignored.
	Interpolate(gs []Grid, mass float64, ids []int64, xs []geom.Vec)
	// Count(gs []Grid, ids []int64, xs []geom.Vec, buf []bool)
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
	man *catalog.ParticleManager
	countWidth int64
	steps int

	gen *rand.Generator

	// Buffers
	idxBuf geom.TetraIdxs
	tet geom.Tetra
	randBuf []float64
	vecBuf []geom.Vec
}

type cellCenter struct {
	man *catalog.ParticleManager
	countWidth int64
	idxBuf geom.TetraIdxs
	tet geom.Tetra
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

func CellCenter(man *catalog.ParticleManager, countWidth int64) Interpolator {
	return &cellCenter{man, countWidth, geom.TetraIdxs{}, geom.Tetra{}}
}

func MonteCarlo(man *catalog.ParticleManager, countWidth int64,
	gen *rand.Generator, steps int) Interpolator {

	return &mcarlo{
		NearestGridPoint(), man, countWidth, steps,
		gen, geom.TetraIdxs{}, geom.Tetra{}, 
		make([]float64, steps * 3), make([]geom.Vec, steps),
	}
}

func SubMonteCarlo( gs []Grid, countWidth int, steps int) Interpolator {
	panic("Not Yet Implemneted.")
}

// Interpolate interpolates a sequence of particles onto a density grid via a
// nearest grid point scheme.
func (intr *ngp) Interpolate(gs []Grid, mass float64, ids []int64, xs []geom.Vec) {
	frac := mass / gs[0].CellVolume
	for _, pt := range xs {
		xp, yp, zp := float64(pt[0]), float64(pt[1]), float64(pt[2])
		xc, yc, zc := cellPoints(xp, yp, zp, gs[0].CellWidth)
		i, j, k := int(xc), int(yc), int(zc)

		for gIdx := range gs {
			g := &gs[gIdx]
			if idx, ok := g.G.IdxCheck(i, j, k); ok {
				g.Rhos[idx] += frac
				continue
			}
		}
	}
}

// Interpolate interpolates a sequence of particles onto a density grid via a
// cloud in cell scheme.
func (intr *cic) Interpolate(gs []Grid, mass float64, ids []int64, xs []geom.Vec) {
	frac := mass / gs[0].CellVolume
	cw, cw2 := gs[0].CellWidth, gs[0].CellWidth / 2
	for _, pt := range xs {
		
		xp, yp, zp := float64(pt[0])-cw2, float64(pt[1])-cw2, float64(pt[2])-cw2
		xc, yc, zc := cellPoints(xp, yp, zp, gs[0].CellWidth)
		dx, dy, dz := (xp / cw)-xc, (yp / cw)-yc, (zp / cw)-zc
		tx, ty, tz := 1-dx, 1-dy, 1-dz

		i0, i1 := gs[0].nbrs(int(xc))
		j0, j1 := gs[0].nbrs(int(yc))
		k0, k1 := gs[0].nbrs(int(zc))

		over000 := tx*ty*tz*frac
		over100 := dx*ty*tz*frac
		over010 := tx*dy*tz*frac
		over110 := dx*dy*tz*frac
		over001 := tx*ty*dz*frac
		over101 := dx*ty*dz*frac
		over011 := tx*dy*dz*frac
		over111 := dx*dy*dz*frac

		for gIdx := range gs {
			g := &gs[gIdx]

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

func (intr *cellCenter) Interpolate(gs []Grid, mass float64, ids []int64, xs []geom.Vec) {
	cb := &geom.CellBounds{}

	for _, id := range ids {
		for dir := 0; dir < 5; dir++ {
			intr.idxBuf.Init(id, intr.countWidth, dir)

			p0 := intr.man.Get(intr.idxBuf[0])
			p1 := intr.man.Get(intr.idxBuf[1])
			p2 := intr.man.Get(intr.idxBuf[2])
			p3 := intr.man.Get(intr.idxBuf[3])

			if p0 == nil || p1 == nil || p2 == nil || p3 == nil {
				log.Printf("Tetrahedron [%v %v %v %v] not in manager.\n",
					p0, p1, p2, p3)

				if p0 == nil {
					log.Printf("Point at index %d not in manager\n",
						intr.idxBuf[0])
				} else if p1 == nil {
					log.Printf("Point at index %d not in manager\n",
						intr.idxBuf[1])
				} else if p2 == nil {
					log.Printf("Point at index %d not in manager\n",
						intr.idxBuf[2])
				} else if p3 == nil {
					log.Printf("Point at index %d not in manager\n",
						intr.idxBuf[3])
				}

				continue
			}

			intr.tet.Init(&p0.Xs, &p1.Xs, &p2.Xs, &p3.Xs, gs[0].BoxWidth)
			intr.tet.CellBoundsAt(gs[0].CellWidth, cb)

			for i := range gs {
				if gs[i].G.Intersect(cb, &gs[i].BG) {
					intr.intrTetra(mass / 6.0, &gs[i], cb)
				}
			}
		}
	}
}

func (intr *cellCenter) intrTetra(mass float64, g *Grid, cb *geom.CellBounds) {
	minX := maxInt(cb.Min[0], g.BG.Origin[0])
	maxX := minInt(cb.Max[0], g.BG.Origin[0] + g.BG.Width)
	minY := maxInt(cb.Min[1], g.BG.Origin[1])
	maxY := minInt(cb.Max[1], g.BG.Origin[1] + g.BG.Width)
	minZ := maxInt(cb.Min[2], g.BG.Origin[2])
	maxZ := minInt(cb.Max[2], g.BG.Origin[2] + g.BG.Width)

	frac := mass * g.CellVolume / intr.tet.Volume()

	for x := minX; x <= maxX; x++ {
		for y := minY; y <= maxY; y++ {
			for z := minZ; z <= maxZ; z++ {
				xIdx, yIdx, zIdx := g.BG.Wrap(x, y, z)
				idx := g.BG.Idx(xIdx, yIdx, zIdx)
				g.Rhos[idx] += frac
			}
		}
	}
}

func maxInt(x, y int) int {
	if x > y {
		return x
	}
	return y
}

func minInt(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func (intr *mcarlo) Interpolate(gs []Grid, mass float64, ids []int64, xs []geom.Vec) {
	ptMass := mass / float64(intr.steps) / 6.0
	intersectGs := make([]Grid, len(gs))
	cb := &geom.CellBounds{}

	for _, id := range ids {
		for dir := 0; dir < 6; dir++ {
			intr.idxBuf.Init(id, intr.countWidth, dir)
			
			p0 := intr.man.Get(intr.idxBuf[0])
			p1 := intr.man.Get(intr.idxBuf[1])
			p2 := intr.man.Get(intr.idxBuf[2])
			p3 := intr.man.Get(intr.idxBuf[3])

			if p0 == nil || p1 == nil || p2 == nil || p3 == nil {
				log.Printf("Tetrahedron [%v %v %v %v] not in manager.\n",
					p0, p1, p2, p3)

				if p0 == nil {
					log.Printf("Point at index %d not in manager\n",
						intr.idxBuf[0])
				} else if p1 == nil {
					log.Printf("Point at index %d not in manager\n",
						intr.idxBuf[1])
				} else if p2 == nil {
					log.Printf("Point at index %d not in manager\n",
						intr.idxBuf[2])
				} else if p3 == nil {
					log.Printf("Point at index %d not in manager\n",
						intr.idxBuf[3])
				}

				continue
			}

			intr.tet.Init(&p0.Xs, &p1.Xs, &p2.Xs, &p3.Xs, gs[0].BoxWidth)
			intr.tet.Sample(intr.gen, intr.randBuf, intr.vecBuf)
			intr.tet.CellBoundsAt(gs[0].CellWidth, cb)

			intersectNum := 0
			for i := range gs {
				if gs[i].G.Intersect(cb, &gs[i].BG) {
					intersectGs[intersectNum] = gs[i]
					intersectNum++
				}
			}
			
			if intersectNum == 0 {
				log.Panic("Got zero intersections!!")
				log.Panicf("%v and %v\n", cb, gs[0].G.CellBounds())
			}

			intr.subIntr.Interpolate(intersectGs[0: intersectNum],
				ptMass, nil, intr.vecBuf)
		}
	}
}
