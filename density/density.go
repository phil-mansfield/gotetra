/*package density interpolates sequences of particle positions onto a density
grid.
*/
package density

import (
	"github.com/phil-mansfield/gotetra/rand"
	"github.com/phil-mansfield/gotetra/geom"
)

const (
	// An arbitrarily chosen prime number indicating the number of precomputed
	// random sequences to use for the MonteCarlo interpolator.
	//
	// Some more thought should be put into this because it gets prohibitively
	// large for high particle counts.
	UnitBuffers = 7919
)

// Interpolator creates a grid-based density distribution from seqeunces of
// positions.
type Interpolator interface {
	// Interpolate adds the density distribution implied by points to the
	// density grid used by the Interpolator. Particles should all be within
	// the bounds of the bounding grid and points not within the interpolation
	// grid will be ignored.
	Interpolate(rhos []float32, cb *geom.CellBounds, mass float64, xs []geom.Vec)
}

type ngp struct { }

// The ordering of these fields makes no goddamned sense.
type mcarlo struct {
	subIntr Interpolator
	segWidth int64
	steps int

	gen *rand.Generator

	skip int64
	// Buffers
	idxBuf geom.TetraIdxs
	tet geom.Tetra

	unitBufs [UnitBuffers][]geom.Vec
	vecBuf []geom.Vec
}

func NearestGridPoint() Interpolator {
	return &ngp{}
}

func MonteCarlo(
	segWidth int64,
	gen *rand.Generator,
	points int,
	skip int64,
) Interpolator {
	mc := &mcarlo{
		NearestGridPoint(), segWidth, points, gen, skip,
		geom.TetraIdxs{}, geom.Tetra{},
		[UnitBuffers][]geom.Vec{},
		make([]geom.Vec, points),
	}

	xs := make([]float64, points)
	ys := make([]float64, points)
	zs := make([]float64, points)

	for i := range mc.unitBufs {
		buf := make([]geom.Vec, points)
		gen.UniformAt(0.0, 1.0, xs)
		gen.UniformAt(0.0, 1.0, ys)
		gen.UniformAt(0.0, 1.0, zs)
		geom.DistributeUnit(xs, ys, zs, buf)
		
		mc.unitBufs[i] = buf
	}

	return mc
}

// Interpolate interpolates a sequence of particles onto a density grid via a
// nearest grid point scheme.
func (intr *ngp) Interpolate(
	rhos []float32, cb *geom.CellBounds, mass float64, xs []geom.Vec,
) {
	length := cb.Width[0]
	area := cb.Width[0] * cb.Width[1]
	m := float32(mass)
	for _, pt := range xs {
		i := int(pt[0])
		j := int(pt[1])
		k := int(pt[2])
		rhos[i + j * length + k * area] += m
	}
}

func (intr *mcarlo) Interpolate(
	rhos []float32, cb *geom.CellBounds, mass float64, xs []geom.Vec,
) {
	segWidth := intr.segWidth
	gridWidth := segWidth + 1
	idxWidth := intr.segWidth / intr.skip

	ptMass := mass / float64(len(intr.vecBuf)) / 6.0 *
		float64(intr.skip * intr.skip * intr.skip)

	for z := int64(0); z < idxWidth; z++ {
		for y := int64(0); y < idxWidth; y++ {
			for x := int64(0); x < idxWidth; x++ {

				idx := (x * intr.skip) +
					(y * intr.skip) * gridWidth +
					(z * intr.skip) * gridWidth * gridWidth

				for dir := 0; dir < 6; dir++ {
					intr.idxBuf.Init(idx, gridWidth, intr.skip, dir)
					
					intr.tet.Init(
						&xs[intr.idxBuf[0]],
						&xs[intr.idxBuf[1]],
						&xs[intr.idxBuf[2]],
						&xs[intr.idxBuf[3]],
						1e6,
					)
					
					intr.tet.DistributeTetra(
						intr.unitBufs[idx % UnitBuffers],
						intr.vecBuf,
					)

					intr.subIntr.Interpolate(rhos, cb, ptMass, intr.vecBuf)
				}
			}
		}
	}
}

// AddBuffer adds the contents of a density buffer constrained by the given
// CellBounds to a periodic grid with the given number of cells.
func AddBuffer(grid, buf []float32, cb *geom.CellBounds, cells int) {
	for z := 0; z < cb.Width[2]; z++ {
		zBufIdx := z * cb.Width[0] * cb.Width[1]
		zGridIdx := ((z + cb.Origin[2]) % cells) * cells * cells
		for y := 0; y < cb.Width[1]; y++ {
			yBufIdx := y * cb.Width[0]
			yGridIdx := ((y + cb.Origin[1]) % cells) * cells

			for x := 0; x < cb.Width[0]; x++ {
				xBufIdx := x
				xGridIdx := x + cb.Origin[0]
				if xGridIdx >= cells { xGridIdx -= cells }

				bufIdx := xBufIdx + yBufIdx + zBufIdx
				gridIdx := xGridIdx + yGridIdx + zGridIdx
				
				grid[gridIdx] += buf[bufIdx]
			}
		}
	}
}

func min(x, y int) int {
	if x < y { return x }
	return y
}

func max(x, y int) int {
	if x < y { return y }
	return x
}
