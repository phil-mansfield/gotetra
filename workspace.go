package gotetra

import (
	"math"

	"github.com/phil-mansfield/gotetra/density"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/io"
)

// The workflow here is that the user gives us Boxes that we need to
// interpolate onto and which represent a unique memory footprint. Then, for
// every sheet segment which intersects with this box, we create an immutable
// thread-safe Workspace which abstracts over 

////////////////
// Interfaces //
////////////////

type Workspace interface {
	// BufferSize calculates buffer size required to represent the
	// underlying grid.
	BufferSize() int

	// ScaleVecs converts a vector array into the workspace's code units.
	ScaleVecs(vs []geom.Vec)

	// Interpolate performs a nearest grid point interpolation on the given
	// position vectors onto the grid self.
	Interpolate(
		buf []float64,
		xs []geom.Vec, cb *geom.CellBounds,
		ptVal float64, low, high int,
	)

	// Add adds the contents of buf to grid where buf is the workspace grid and
	// grid is the domain grid. The domain grid is contained within the given
	// cell bounds.
	Add(buf, grid []float64, cb *geom.CellBounds)
}

type Box interface {
	Workspace(hd *io.SheetHeader) Workspace

	CellSpan() [3]int
	CellOrigin() [3]int
	CellWidth() float64

	Vals() []float64

	ProjectionAxis() (dim int, ok bool)
}

/////////////////////////
// Box implementations //
/////////////////////////

// TODO: Standardize the usage of 'cells' vs. 'pixels' or 'totalPixels' or
// whatever. Also 'cellWidth' vs. 'boxWidth'
type baseBox struct {
	cb geom.CellBounds
	vals []float64
	cells int
	cellWidth float64
}

func (b *baseBox) CellSpan() [3]int { return b.cb.Width }
func (b *baseBox) CellOrigin() [3]int { return b.cb.Origin }
func (b *baseBox) CellWidth() float64 { return b.cellWidth }
func (b *baseBox) Vals() []float64 { return b.vals }

type box2D struct {
	baseBox
	proj int
}

func (b *box2D) Workspace(hd *io.SheetHeader) Workspace {
	seg := &segmentWorkspace2D{ }
	dom := &domainWorkspace2D{ }

	seg.boxWidth = b.cellWidth * float64(b.cells)
	dom.boxWidth = b.cellWidth * float64(b.cells)
	seg.totalPixels = b.cells
	dom.totalPixels = b.cells
	seg.proj = b.proj
	dom.proj = b.proj

	seg.cb = *hd.CellBounds(b.cells)
	dom.cb = b.cb

	if seg.BufferSize() < dom.BufferSize() {
		return seg
	} else {
		return dom
	}
}

func (b *box2D) ProjectionAxis() (int, bool) { return b.proj, true }

type box3D struct {
	baseBox
}

func (b *box3D) Workspace(hd *io.SheetHeader) Workspace {
	seg := &segmentWorkspace2D{ }
	dom := &domainWorkspace2D{ }

	seg.boxWidth = b.cellWidth * float64(b.cells)
	dom.boxWidth = b.cellWidth * float64(b.cells)
	seg.totalPixels = b.cells
	dom.totalPixels = b.cells

	seg.cb = *hd.CellBounds(b.cells)
	dom.cb = b.cb

	if seg.BufferSize() < dom.BufferSize() {
		return seg
	} else {
		return dom
	}
}

func (b *box3D) ProjectionAxis() (int, bool) { return -1, false }


// NewBox creates a grid and a wrapper for the redering box defined by the
// given config file, and which lives inside a simulation box with the given
// width and pixel count.
func NewBox(boxWidth float64, totalPixels int, config *io.BoxConfig) Box {
	if config.IsProjection() {
		return newBox2D(boxWidth, totalPixels, config)
	} else {
		return newBox3D(boxWidth, totalPixels, config)
	}
}

func newBox2D(boxWidth float64, totalPixels int, config *io.BoxConfig) Box {
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

	cellWidth := boxWidth / float64(totalPixels)
	origin := [3]float64{ config.X, config.Y, config.Z }
	width := [3]float64{ config.XWidth, config.YWidth, config.ZWidth }

	for j := 0; j < 3; j++ {
		b.cb.Origin[j] = int(math.Floor(origin[j] / cellWidth))
		b.cb.Width[j] = 1 + int(
			math.Floor((width[j] + origin[j]) / cellWidth),
		)
		b.cb.Width[j] -= b.cb.Origin[j]
	}

	b.cells = totalPixels
	b.cellWidth = cellWidth

	iDim, jDim := 0, 1
	if b.proj == 0 { iDim, jDim = 1, 2 }
	if b.proj == 1 { iDim, jDim = 0, 2 }
	b.vals = make([]float64, b.cb.Width[iDim] * b.cb.Width[jDim])

	return b
}

func newBox3D(boxWidth float64, totalPixels int, config *io.BoxConfig) Box {
	// TODO: Rewrite for code reuse.

	b := new(box3D)

	cellWidth := boxWidth / float64(totalPixels)
	origin := [3]float64{ config.X, config.Y, config.Z }
	width := [3]float64{ config.XWidth, config.YWidth, config.ZWidth }

	for j := 0; j < 3; j++ {
		b.cb.Origin[j] = int(math.Floor(origin[j] / cellWidth))
		b.cb.Width[j] = 1 + int(
			math.Floor((width[j] + origin[j]) / cellWidth),
		)
		b.cb.Width[j] -= b.cb.Origin[j]
	}

	b.cells = totalPixels
	b.cellWidth = cellWidth

	b.vals = make([]float64, b.cb.Width[0] * b.cb.Width[1] * b.cb.Width[2])

	return b
}

///////////////////////////////
// Workspace Implementations //
///////////////////////////////

type baseWorkspace struct {
	cb geom.CellBounds
	totalPixels int
	boxWidth float64
}

func (w *baseWorkspace) ScaleVecs(vs []geom.Vec) {
	w.cb.ScaleVecs(vs, w.totalPixels, w.boxWidth)
}

type baseWorkspace2D struct {
	baseWorkspace
	proj int
}

func (w *baseWorkspace2D) BufferSize() int {
	// I'm not doing this the obvious way (with a division) to avoid
	// integer overflow.
	prod := 1
	for i := 0; i < 3; i++ {
		if i == w.proj { continue }
		prod *= w.cb.Width[i]
	}

	return prod
}

func (w *baseWorkspace2D) Add(buf, grid []float64, cb *geom.CellBounds) {
	// Since the majority of the program's serial time is spent here, there are
	// some completely stupid-looking optimizations here.

	iDim, jDim := 0, 1
	if w.proj == 0 { iDim, jDim = 1, 2 }	
	if w.proj == 1 { iDim, jDim = 0, 2 }

	diOrigin := w.cb.Origin[iDim] - cb.Origin[iDim]
	djOrigin := w.cb.Origin[jDim] - cb.Origin[jDim]
	if diOrigin >= w.totalPixels { diOrigin -= w.totalPixels }
	if djOrigin >= w.totalPixels { djOrigin -= w.totalPixels }
	iWorkLow, iWorkHigh := diOrigin, diOrigin + w.cb.Width[iDim] - 1
	jWorkLow, jWorkHigh := djOrigin, djOrigin + w.cb.Width[jDim] - 1

	// TODO: Speed this up by moving the conditionals out of the loops. It
	// shouldn't be that big of a deal since there will be exactly one branch
	// misprediction per loop.
	jWork := 0
	for j := jWorkLow; j <= jWorkHigh; j++ {
		jDomain := j
		if jDomain >= w.totalPixels { jDomain -= w.totalPixels }
		jWorkFlat   := jWork   * w.cb.Width[iDim]
		jDomainFlat := jDomain * cb.Width[iDim]

		iWork := 0
		for i := iWorkLow; i <= iWorkHigh; i++ {
			iDomain := i
			if iDomain >= w.totalPixels { iDomain -= w.totalPixels }
			grid[iDomain + jDomainFlat] = buf[iWork + jWorkFlat]

			iWork++
		}
		jWork++
	}
}

type baseWorkspace3D struct {
	baseWorkspace
}

func (w *baseWorkspace3D) BufferSize() int {
	return w.cb.Width[0] * w.cb.Width[1] * w.cb.Width[2]
}

func (w *baseWorkspace3D) Add(buf, grid []float64, cb *geom.CellBounds) {
	// Since the majority of the program's serial time is spent here, there are
	// some completely stupid-looking optimizations here.

	dxOrigin := w.cb.Origin[0] - cb.Origin[0]
	dyOrigin := w.cb.Origin[1] - cb.Origin[1]
	dzOrigin := w.cb.Origin[2] - cb.Origin[2]
	if dxOrigin >= w.totalPixels { dxOrigin -= w.totalPixels }
	if dyOrigin >= w.totalPixels { dyOrigin -= w.totalPixels }
	if dzOrigin >= w.totalPixels { dzOrigin -= w.totalPixels }
	xWorkLow, xWorkHigh := dxOrigin, dxOrigin + w.cb.Width[0] - 1
	yWorkLow, yWorkHigh := dyOrigin, dyOrigin + w.cb.Width[1] - 1
	zWorkLow, zWorkHigh := dzOrigin, dzOrigin + w.cb.Width[2] - 1

	// TODO: Speed this up by moving the conditionals out of the loops. It
	// shouldn't be that big of a deal since there will be exactly one branch
	// misprediction.
	zWork := 0
	for z := zWorkLow; z <= zWorkHigh; z++ {
		zDomain := z
		if zDomain >= w.totalPixels { zDomain -= w.totalPixels }
		zWorkFlat   := zWork   * w.cb.Width[1] * w.cb.Width[0]
		zDomainFlat := zDomain * cb.Width[1]   * cb.Width[0]

		yWork := 0
		for y := yWorkLow; y <= yWorkHigh; y++ {
			yDomain := y
			if yDomain >= w.totalPixels { yDomain -= w.totalPixels }
			yWorkFlat   := yWork   * w.cb.Width[0]
			yDomainFlat := yDomain * cb.Width[0]

			xWork := 0
			for x := xWorkLow; x <= xWorkHigh; x++ {
				xDomain := x
				if xDomain >= w.totalPixels { xDomain -= w.totalPixels }

				grid[xDomain + yDomainFlat + zDomainFlat] = 
					buf[xWork + yWorkFlat + zWorkFlat]
				
				xWork++
			}
			yWork++
		}
		zWork++
	}
}

type segmentWorkspace2D struct {
	baseWorkspace2D
}

func (w *segmentWorkspace2D) Interpolate(
	buf []float64,
	xs []geom.Vec, cb *geom.CellBounds,
	ptVal float64, low, high int,
) {
	// We're garuanteed that every point in xs must fall in bounds. cb is
	// not used.
	iDim, jDim := 0, 1
	if w.proj == 0 { iDim, jDim = 1, 2 }
	if w.proj == 1 { iDim, jDim = 0, 2 }

	length := w.cb.Width[iDim]
	for idx := low; idx < high; idx++ {
		pt := xs[idx]
		i, j := int(pt[iDim]), int(pt[jDim])
		buf[i + j * length] += ptVal
	}
}

type segmentWorkspace3D struct {
	baseWorkspace3D
}

func (w *segmentWorkspace3D) Interpolate(
	buf []float64,
	xs []geom.Vec, cb *geom.CellBounds,
	ptVal float64, low, high int,
) {
	// We're garuanteed that every point in xs must fall in bounds. cb is
	// not used.
	length := w.cb.Width[0]
	area :=   w.cb.Width[0] * w.cb.Width[1]
	for idx := low; idx < high; idx++ {
		pt := xs[idx]
		i, j, k := int(pt[0]), int(pt[1]), int(pt[2])
		buf[i + j * length + k * area] += ptVal
	}
}

type domainWorkspace2D struct {
	baseWorkspace2D
}

func (w *domainWorkspace2D) Interpolate(
	buf []float64,
	xs []geom.Vec, cb *geom.CellBounds,
	ptVal float64, low, high int,
) {
	// We are not garuanteed that every point in xs must fall in bounds.

	xDiff, yDiff, zDiff := cbSubtr(cb, &w.cb)

	iDim, jDim, kDim := 0, 1, 2
	iDiff, jDiff, kDiff := xDiff, yDiff, zDiff
	if w.proj == 0 {
		iDim, jDim, kDim = 1, 2, 0
		iDiff, jDiff, kDiff = yDiff, zDiff, xDiff
	} else if w.proj == 1 {
		iDim, jDim, kDim = 0, 2, 1
		iDiff, jDiff, kDiff = xDiff, zDiff, yDiff
	}

	length := cb.Width[iDim]
	for idx := low; idx < high; idx++ {
		pt := xs[idx]

		i := bound(int(pt[iDim]) - iDiff, w.totalPixels)
		if i < cb.Width[iDim] {
			j := bound(int(pt[jDim]) - jDiff, w.totalPixels)
			if j < cb.Width[jDim] {
				k := bound(int(pt[kDim]) - kDiff, w.totalPixels)
				if k < cb.Width[kDim] {
					buf[i + j * length] += ptVal
				}
			}
		}
	}
}

type domainWorkspace3D struct {
	baseWorkspace3D
}

func (w *domainWorkspace3D) Interpolate(
	buf []float64,
	xs []geom.Vec, cb *geom.CellBounds,
	ptVal float64, low, high int,
) {
	// We are not garuanteed that every point in xs must fall in bounds.
	length := cb.Width[0]
	area :=   cb.Width[0] * cb.Width[1]

	iDiff, jDiff, kDiff := cbSubtr(cb, &w.cb)

	for idx := low; idx < high; idx++ {
		pt := xs[idx]
		
		i := bound(int(pt[0]) - iDiff, w.totalPixels)
		if i < cb.Width[0] {
			j := bound(int(pt[1]) - jDiff, w.totalPixels)
			if j < cb.Width[1] {
				k := bound(int(pt[2]) - kDiff, w.totalPixels)
				if k < cb.Width[2] {
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
	_ Workspace = &segmentWorkspace2D{ }
	_ Workspace = &domainWorkspace3D{ }
	_ Workspace = &segmentWorkspace2D{ }
	_ Workspace = &domainWorkspace3D{ }

	_ density.Interpolator = &segmentWorkspace2D{ }
	_ density.Interpolator = &domainWorkspace3D{ }
	_ density.Interpolator = &segmentWorkspace2D{ }
	_ density.Interpolator = &domainWorkspace3D{ }

	_ Box = &box2D{ }
	_ Box = &box3D{ }
)
