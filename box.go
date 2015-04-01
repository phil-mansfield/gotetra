package gotetra

import (
	"math"

	"github.com/phil-mansfield/gotetra/density"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/io"
)

////////////////
// Interfaces //
////////////////

type Box interface {
	Overlap(hd *io.SheetHeader) Overlap

	CellSpan() [3]int
	CellOrigin() [3]int
	CellWidth() float64
	Cells() int
	Points() int

	Vals() []float64

	ProjectionAxis() (dim int, ok bool)
}


type Overlap interface {
	// BufferSize calculates buffer size required to represent the
	// underlying grid.
	BufferSize() int

	// ScaleVecs converts a vector array into the overlap's code units.
	ScaleVecs(vs []geom.Vec)

	// Interpolate performs a nearest grid point interpolation on the given
	// position vectors onto the grid self.
	Interpolate(
		rhos []float64, rhoCb *geom.CellBounds,
		xs []geom.Vec, xCb *geom.CellBounds,
		ptVal float64, low, high int,
	)

	// Add adds the contents of buf to grid where buf is the overlap grid and
	// grid is the domain grid. The domain grid is contained within the given
	// cell bounds.
	Add(buf, grid []float64, cb *geom.CellBounds)
}

/////////////////////////
// Box implementations //
/////////////////////////

type baseBox struct {
	cb geom.CellBounds
	vals []float64
	pts, cells int
	cellWidth float64
}

func (b *baseBox) CellSpan() [3]int { return b.cb.Width }
func (b *baseBox) CellOrigin() [3]int { return b.cb.Origin }
func (b *baseBox) CellWidth() float64 { return b.cellWidth }
func (b *baseBox) Cells() int {return b.cells }
func (b *baseBox) Points() int { return b.pts }
func (b *baseBox) Vals() []float64 { return b.vals }

type box2D struct {
	baseBox
	proj int
}

func (b *box2D) Overlap(hd *io.SheetHeader) Overlap {
	seg := &segmentOverlap2D{ }
	dom := &domainOverlap2D{ }

	seg.boxWidth = b.cellWidth * float64(b.cells)
	dom.boxWidth = b.cellWidth * float64(b.cells)
	seg.cells = b.cells
	dom.cells = b.cells
	seg.proj = b.proj
	dom.proj = b.proj

	seg.cb = *hd.CellBounds(b.cells)
	dom.cb = b.cb

	if seg.BufferSize() <= dom.BufferSize() {
		return seg
	} else {
		return dom
	}
}

func (b *box2D) ProjectionAxis() (int, bool) { return b.proj, true }

type box3D struct {
	baseBox
}

func (b *box3D) Overlap(hd *io.SheetHeader) Overlap {
	seg := &segmentOverlap2D{ }
	dom := &domainOverlap2D{ }

	seg.boxWidth = b.cellWidth * float64(b.cells)
	dom.boxWidth = b.cellWidth * float64(b.cells)
	seg.cells = b.cells
	dom.cells = b.cells

	seg.cb = *hd.CellBounds(b.cells)
	dom.cb = b.cb

	if seg.BufferSize() <= dom.BufferSize() {
		return seg
	} else {
		return dom
	}
}

func (b *box3D) ProjectionAxis() (int, bool) { return -1, false }


// NewBox creates a grid and a wrapper for the redering box defined by the
// given config file, and which lives inside a simulation box with the given
// width and pixel count.
func NewBox(boxWidth float64, pts, cells int, config *io.BoxConfig) Box {
	if config.IsProjection() {
		return newBox2D(boxWidth, pts, cells, config)
	} else {
		return newBox3D(boxWidth, pts, cells, config)
	}
}

func newBox2D(boxWidth float64, pts, cells int, config *io.BoxConfig) Box {
	// TODO: Rewrite for code reuse.

	b := new(box2D)

	if config.ProjectionAxis == "X" {
		b.proj = 0
	} else if config.ProjectionAxis == "Y" {
		b.proj = 1
	} else if config.ProjectionAxis == "Z" {
		b.proj = 2
	} else {
		panic("Internal flag inconsistency.")
	}

	cellWidth := boxWidth / float64(cells)
	origin := [3]float64{ config.X, config.Y, config.Z }
	width := [3]float64{ config.XWidth, config.YWidth, config.ZWidth }

	for j := 0; j < 3; j++ {
		b.cb.Origin[j] = int(math.Floor(origin[j] / cellWidth))
		b.cb.Width[j] = 1 + int(
			math.Floor((width[j] + origin[j]) / cellWidth),
		)
		b.cb.Width[j] -= b.cb.Origin[j]
	}

	b.cells = cells
	b.pts = pts
	b.cellWidth = cellWidth

	iDim, jDim := 0, 1
	if b.proj == 0 { iDim, jDim = 1, 2 }
	if b.proj == 1 { iDim, jDim = 0, 2 }
	b.vals = make([]float64, b.cb.Width[iDim] * b.cb.Width[jDim])

	return b
}

func newBox3D(boxWidth float64, pts, cells int, config *io.BoxConfig) Box {
	// TODO: Rewrite for code reuse.

	b := new(box3D)

	cellWidth := boxWidth / float64(cells)
	origin := [3]float64{ config.X, config.Y, config.Z }
	width := [3]float64{ config.XWidth, config.YWidth, config.ZWidth }

	for j := 0; j < 3; j++ {
		b.cb.Origin[j] = int(math.Floor(origin[j] / cellWidth))
		b.cb.Width[j] = 1 + int(
			math.Floor((width[j] + origin[j]) / cellWidth),
		)
		b.cb.Width[j] -= b.cb.Origin[j]
	}

	b.cells = cells
	b.pts = pts
	b.cellWidth = cellWidth

	b.vals = make([]float64, b.cb.Width[0] * b.cb.Width[1] * b.cb.Width[2])

	return b
}

/////////////////////////////
// Overlap Implementations //
/////////////////////////////

type baseOverlap struct {
	cb geom.CellBounds
	cells int
	boxWidth float64
}

func (w *baseOverlap) ScaleVecs(vs []geom.Vec) {
	w.cb.ScaleVecs(vs, w.cells, w.boxWidth)
}

type baseOverlap2D struct {
	baseOverlap
	proj int
}

func (w *baseOverlap2D) BufferSize() int {
	// I'm not doing this the obvious way (with a division) to avoid
	// integer overflow.
	prod := 1
	for i := 0; i < 3; i++ {
		if i == w.proj { continue }
		prod *= w.cb.Width[i]
	}

	return prod
}

func (w *baseOverlap2D) Add(buf, grid []float64, cb *geom.CellBounds) {
	// Since the majority of the program's serial time is spent here, there are
	// some completely stupid-looking optimizations here.

	// TODO: Change variables names for *Work* to *Over*, to reflect the name
	// change fo the type. I just don't want to do this right now...

	iDim, jDim := 0, 1
	if w.proj == 0 { iDim, jDim = 1, 2 }	
	if w.proj == 1 { iDim, jDim = 0, 2 }

	diOrigin := w.cb.Origin[iDim] - cb.Origin[iDim]
	djOrigin := w.cb.Origin[jDim] - cb.Origin[jDim]
	if diOrigin >= w.cells { diOrigin -= w.cells }
	if djOrigin >= w.cells { djOrigin -= w.cells }
	iWorkLow, iWorkHigh := diOrigin, diOrigin + w.cb.Width[iDim] - 1
	jWorkLow, jWorkHigh := djOrigin, djOrigin + w.cb.Width[jDim] - 1

	// TODO: Speed this up by moving the conditionals out of the loops. It
	// shouldn't be that big of a deal since there will be exactly one branch
	// misprediction per loop.
	jWork := 0
	for j := jWorkLow; j <= jWorkHigh; j++ {
		jDomain := j
		if jDomain >= w.cells { jDomain -= w.cells }
		jWorkFlat   := jWork   * w.cb.Width[iDim]
		jDomainFlat := jDomain * cb.Width[iDim]

		iWork := 0
		for i := iWorkLow; i <= iWorkHigh; i++ {
			iDomain := i
			if iDomain >= w.cells { iDomain -= w.cells }
			grid[iDomain + jDomainFlat] = buf[iWork + jWorkFlat]

			iWork++
		}
		jWork++
	}
}

type baseOverlap3D struct {
	baseOverlap
}

func (w *baseOverlap3D) BufferSize() int {
	return w.cb.Width[0] * w.cb.Width[1] * w.cb.Width[2]
}

func (w *baseOverlap3D) Add(buf, grid []float64, cb *geom.CellBounds) {
	// Since the majority of the program's serial time is spent here, there are
	// some completely stupid-looking optimizations here.

	dxOrigin := w.cb.Origin[0] - cb.Origin[0]
	dyOrigin := w.cb.Origin[1] - cb.Origin[1]
	dzOrigin := w.cb.Origin[2] - cb.Origin[2]
	if dxOrigin >= w.cells { dxOrigin -= w.cells }
	if dyOrigin >= w.cells { dyOrigin -= w.cells }
	if dzOrigin >= w.cells { dzOrigin -= w.cells }
	xWorkLow, xWorkHigh := dxOrigin, dxOrigin + w.cb.Width[0] - 1
	yWorkLow, yWorkHigh := dyOrigin, dyOrigin + w.cb.Width[1] - 1
	zWorkLow, zWorkHigh := dzOrigin, dzOrigin + w.cb.Width[2] - 1

	// TODO: Speed this up by moving the conditionals out of the loops. It
	// shouldn't be that big of a deal since there will be exactly one branch
	// misprediction.
	zWork := 0
	for z := zWorkLow; z <= zWorkHigh; z++ {
		zDomain := z
		if zDomain >= w.cells { zDomain -= w.cells }
		zWorkFlat   := zWork   * w.cb.Width[1] * w.cb.Width[0]
		zDomainFlat := zDomain * cb.Width[1]   * cb.Width[0]

		yWork := 0
		for y := yWorkLow; y <= yWorkHigh; y++ {
			yDomain := y
			if yDomain >= w.cells { yDomain -= w.cells }
			yWorkFlat   := yWork   * w.cb.Width[0]
			yDomainFlat := yDomain * cb.Width[0]

			xWork := 0
			for x := xWorkLow; x <= xWorkHigh; x++ {
				xDomain := x
				if xDomain >= w.cells { xDomain -= w.cells }

				grid[xDomain + yDomainFlat + zDomainFlat] = 
					buf[xWork + yWorkFlat + zWorkFlat]
				
				xWork++
			}
			yWork++
		}
		zWork++
	}
}

type segmentOverlap2D struct {
	baseOverlap2D
}

func (w *segmentOverlap2D) Interpolate(
	buf []float64, bufCb *geom.CellBounds,
	xs []geom.Vec, xCb *geom.CellBounds,
	ptVal float64, low, high int,
) {
	xDiff, yDiff, zDiff := cbSubtr(bufCb, xCb)
	iDim, jDim, kDim := 0, 1, 2
	iDiff, jDiff, kDiff := xDiff, yDiff, zDiff
	if w.proj == 0 {
		iDim, jDim, kDim = 1, 2, 0
		iDiff, jDiff, kDiff = yDiff, zDiff, xDiff
	} else if w.proj == 1 {
		iDim, jDim, kDim = 0, 2, 1
		iDiff, jDiff, kDiff = xDiff, zDiff, yDiff
	}

	ptVal /= float64(bufCb.Width[kDim])
	length := xCb.Width[iDim]

	for idx := low; idx < high; idx++ {
		pt := xs[idx]
		
		i := int(pt[iDim])
		ii := bound(i - iDiff, w.cells)
		if ii < bufCb.Width[iDim] && ii >= 0 {
			j := int(pt[jDim])
			jj := bound(j - jDiff, w.cells)
			if jj < bufCb.Width[jDim] && jj >= 0 {
				k := int(pt[kDim])
				kk := bound(k - kDiff, w.cells)
				if kk < bufCb.Width[kDim] && kk >= 0 {
					buf[i + j * length] += ptVal
				}
			}
		}
	}
}

type segmentOverlap3D struct {
	baseOverlap3D
}

// TODO: change from 'buf' and 'x' to 'dom' and 'seg'

func (w *segmentOverlap3D) Interpolate(
	buf []float64, bufCb *geom.CellBounds,
	xs []geom.Vec, xCb *geom.CellBounds,
	ptVal float64, low, high int,
) {
	length := xCb.Width[0]
	area :=   xCb.Width[0] * xCb.Width[1]

	xDiff, yDiff, zDiff := cbSubtr(bufCb, xCb)

	for idx := low; idx < high; idx++ {
		pt := xs[idx]
		
		x := int(pt[0])
		i := bound(x - xDiff, w.cells)
		if i < bufCb.Width[0] && i >= 0 {
			y := int(pt[1])
			j := bound(y - yDiff, w.cells)
			if j < bufCb.Width[1] && j >= 0 {
				z := int(pt[2])
				k := bound(z - zDiff, w.cells)
				if k < bufCb.Width[2] && k >= 0 {
					buf[x + y * length + z * area] += ptVal
				}
			}
		}
	}
}

type domainOverlap2D struct {
	baseOverlap2D
}

func (w *domainOverlap2D) Interpolate(
	buf []float64, bufCb *geom.CellBounds,
	xs []geom.Vec, cb *geom.CellBounds,
	ptVal float64, low, high int,
) {
	iDim, jDim, kDim := 0, 1, 2
	if w.proj == 0 { iDim, jDim, kDim = 1, 2, 0 }
	if w.proj == 1 { iDim, jDim, kDim = 0, 2, 1 }

	ptVal /= float64(bufCb.Width[kDim])
	length := bufCb.Width[iDim]

	for idx := low; idx < high; idx++ {
		pt := xs[idx]

		i := bound(int(pt[iDim]), w.cells)
		if i < bufCb.Width[iDim] {
			j := bound(int(pt[jDim]), w.cells)
			if j < bufCb.Width[jDim] {
				k := bound(int(pt[kDim]), w.cells)
				if k < bufCb.Width[kDim] {
					buf[i + j * length] += ptVal
				}
			}
		}
	}
}

type domainOverlap3D struct {
	baseOverlap3D
}

func (w *domainOverlap3D) Interpolate(
	buf []float64, bufCb *geom.CellBounds,
	xs []geom.Vec, xCb *geom.CellBounds,
	ptVal float64, low, high int,
) {
	// bufCb should be redundant with w.Cb.

	length := bufCb.Width[0]
	area :=   bufCb.Width[0] * bufCb.Width[1]

	for idx := low; idx < high; idx++ {
		pt := xs[idx]
		
		// If there's a bug in this function, try replacing int() with
		// int(math.Floor())

		i := bound(int(pt[0]), w.cells)
		if i < bufCb.Width[0] {
			j := bound(int(pt[1]), w.cells)
			if j < bufCb.Width[1] {
				k := bound(int(pt[2]), w.cells)
				if k < bufCb.Width[2] {
					buf[i + j * length + k * area] += ptVal
				}
			}
		}
	}
}

func cbSubtr(cb1, cb2 *geom.CellBounds) (i, j, k int) {
	i = cb1.Origin[0] - cb2.Origin[0]
	j = cb1.Origin[1] - cb2.Origin[1]
	k = cb1.Origin[2] - cb2.Origin[2]
	return i, j, k
}

func bound(x, cells int) int {
	if x < 0 {
		return x + cells
	} else if x >= cells {
		return x - cells
	}
	return x
}

// Typechecking
var (
	_ Overlap = &segmentOverlap2D{ }
	_ Overlap = &domainOverlap3D{ }
	_ Overlap = &segmentOverlap2D{ }
	_ Overlap = &domainOverlap3D{ }

	_ density.Interpolator = &segmentOverlap2D{ }
	_ density.Interpolator = &domainOverlap3D{ }
	_ density.Interpolator = &segmentOverlap2D{ }
	_ density.Interpolator = &domainOverlap3D{ }

	_ Box = &box2D{ }
	_ Box = &box3D{ }
)
