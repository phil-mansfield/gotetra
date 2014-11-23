/*package density interpolates sequences of particle positions onto a density
grid.
*/
package density

import (
	"fmt"
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
	Interpolate(mass float64, ids []int64, xs []geom.Vec)
}

// Flag indicates the interpolation scheme that should be used to assign
// densities.
type Flag int

type cic struct {
	g, bg      geom.Grid
	cellWidth  float64
	cellVolume float64
	rhos       []float64
}

type ngp struct {
	g          geom.Grid
	cellWidth  float64
	cellVolume float64
	rhos       []float64
}

type mcarlo struct {
	subIntr Interpolator
	man *catalog.ParticleManager
	totalWidth float64
	countWidth int64
	steps int

	gen *rand.Generator

	// Buffers
	idxBuf geom.TetraIdxs
	tet geom.Tetra
	randBuf []float64
	vecBuf []geom.Vec
}

const (
	CloudInCell Flag = iota
	NearestGridPoint
	MonteCarlo
)

const (
	DefaultMonteCarloSteps = 100
	DefaultMonteCarloGenerator = rand.Tausworthe
)

// Bounds returns a large bounding grid and a smaller interpolation Grid
// which acts as a single subcell of the bounding Grid. cells gives the number
// of cells in the interpolation Grid, and gridWidth gives the number of
// interpolation Grids within the bounding grid [on one side].
func Bounds(cells, gridWidth, gx, gy, gz int) (g, bg *geom.Grid) {
	g = geom.NewGrid(&[3]int{gx * cells, gy * cells, gz * cells}, cells)
	bg = geom.NewGrid(&[3]int{0, 0, 0}, cells*gridWidth)
	return g, bg
}

// NewInterpolator creates an Interpolator instance using the given
// interpolation scheme which adds to the grid rhos. rhos has boundaries
// defined by the Grid g which is embedded in the bounding Grid bg. These two
// grids will almost always be possible to create through a call to Bounds. The
// variable width refers to the interpolation grid, not the bounding grid.
//
// TODO: this initializer is really gross and in some cases requires values that
// don't even get used. Rewrite this later.
func NewInterpolator(
	flag Flag, g, bg *geom.Grid,
	width float64, countWidth int64,
	man *catalog.ParticleManager,
	rhos []float64,
) Interpolator {
	if g.Volume != len(rhos) {
		panic(fmt.Sprintf(
			"Volume of grid, %d, does not equal length of rhos, %d.",
			g.Volume, len(rhos),
		))
	}

	cellWidth := width / float64(g.Width)
	cellVolume := cellWidth * cellWidth * cellWidth

	switch flag {
	case CloudInCell:
		return &cic{*g, *bg, cellWidth, cellVolume, rhos}
	case NearestGridPoint:
		return &ngp{*g, cellWidth, cellVolume, rhos}
	case MonteCarlo:
		return NewMonteCarloInterpolator(
			g, bg, width, countWidth, man, DefaultMonteCarloSteps, rhos,
		)
	}
	panic(fmt.Sprintf("Unknown flag %d", flag))
}

// NewMonteCarloInterpolator creates a new Interpolator using the specified
// number of interpolation points.
func NewMonteCarloInterpolator(
	g, bg *geom.Grid,
	width float64, countWidth int64,
	man *catalog.ParticleManager,
	steps int,
	rhos []float64,
) Interpolator {

	subIntr := NewInterpolator(NearestGridPoint, g, bg, width, countWidth, man, rhos)
	totalWidth := width * float64(bg.Width) / float64(g.Width)

	return &mcarlo{
		subIntr, man, totalWidth, countWidth, steps,
		rand.NewTimeSeed(DefaultMonteCarloGenerator),
		geom.TetraIdxs{}, geom.Tetra{}, 
		make([]float64, steps * 3), make([]geom.Vec, steps),
	}
}

// Interpolate interpolates a sequence of particles onto a density grid via a
// nearest grid point scheme.
func (intr *ngp) Interpolate(mass float64, ids []int64, xs []geom.Vec) {
	frac := mass / intr.cellVolume
	for _, pt := range xs {
		xp, yp, zp := float64(pt[0]), float64(pt[1]), float64(pt[2])
		xc, yc, zc := cellPoints(xp, yp, zp, intr.cellWidth)
		i, j, k := int(xc), int(yc), int(zc)

		if idx, ok := intr.g.IdxCheck(i, j, k); ok {
			intr.rhos[idx] += frac
		}
	}
}

// Interpolate interpolates a sequence of particles onto a density grid via a
// cloud in cell scheme.
func (intr *cic) Interpolate(mass float64, ids []int64, xs []geom.Vec) {
	frac := mass / intr.cellVolume
	cw, cw2 := intr.cellWidth, intr.cellWidth / 2
	for _, pt := range xs {
		
		xp, yp, zp := float64(pt[0])-cw2, float64(pt[1])-cw2, float64(pt[2])-cw2
		xc, yc, zc := cellPoints(xp, yp, zp, intr.cellWidth)
		dx, dy, dz := (xp / cw)-xc, (yp / cw)-yc, (zp / cw)-zc
		tx, ty, tz := 1-dx, 1-dy, 1-dz

		i0, i1 := intr.nbrs(int(xc))
		j0, j1 := intr.nbrs(int(yc))
		k0, k1 := intr.nbrs(int(zc))

		intr.incr(i0, j0, k0, tx*ty*tz*frac)
		intr.incr(i1, j0, k0, dx*ty*tz*frac)
		intr.incr(i0, j1, k0, tx*dy*tz*frac)
		intr.incr(i1, j1, k0, dx*dy*tz*frac)
		intr.incr(i0, j0, k1, tx*ty*dz*frac)
		intr.incr(i1, j0, k1, dx*ty*dz*frac)
		intr.incr(i0, j1, k1, tx*dy*dz*frac)
		intr.incr(i1, j1, k1, dx*dy*dz*frac)
	}
}

func (intr *cic) nbrs(i int) (i0, i1 int) {
	if i == -1 {
		return intr.bg.Width - 1, 0
	}
	if i+1 == intr.bg.Width {
		return i, 0
	}
	return i, i + 1
}

func (intr *cic) incr(i, j, k int, frac float64) {
	if idx, ok := intr.g.IdxCheck(i, j, k); ok {
		intr.rhos[idx] += frac
	}
}

func cellPoints(x, y, z, cw float64) (xc, yc, zc float64) {
	return math.Floor(x / cw), math.Floor(y / cw), math.Floor(z / cw)
}

func (intr *mcarlo) Interpolate(mass float64, ids []int64, xs []geom.Vec) {
	ptMass := mass / float64(intr.steps) / 6.0

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

			intr.tet.Init(&p0.Xs, &p1.Xs, &p2.Xs, &p3.Xs, intr.totalWidth)
			intr.tet.Sample(intr.gen, intr.randBuf, intr.vecBuf)
			intr.subIntr.Interpolate(ptMass, nil, intr.vecBuf)
		}
	}
}
