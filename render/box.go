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

	Vals() density.Buffer

	ProjectionAxis() (dim int, ok bool)
}


type Overlap interface {
	density.Interpolator

	// BufferSize calculates buffer size required to represent the
	// underlying grid.
	BufferSize() int

	// ScaleVecs converts a vector array into the overlap's code units.
	ScaleVecs(vs []geom.Vec)

	// Add adds the contents of buf to grid where buf is the overlap grid and
	// grid is the domain grid. The domain grid is contained within the given
	// cell bounds.
	Add(buf, grid density.Buffer)
}

/////////////////////////
// Box implementations //
/////////////////////////

type baseBox struct {
	cb geom.CellBounds
	bvals density.Buffer
	pts, cells int
	cellWidth float64
	q density.Quantity
}

func (b *baseBox) CellOrigin() [3]int { return b.cb.Origin }
func (b *baseBox) CellSpan() [3]int { return b.cb.Width }
func (b *baseBox) CellWidth() float64 { return b.cellWidth }
func (b *baseBox) Cells() int {return b.cells }
func (b *baseBox) Points() int { return b.pts }
func (b *baseBox) Vals() density.Buffer { return b.bvals }

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

	dom.domCb = b.cb
	seg.domCb = b.cb
	dom.bufCb = b.cb
	seg.bufCb = *hd.CellBounds(b.cells)

	if seg.BufferSize() <= dom.BufferSize() || 
		(dom.bufCb.Width[0] >= b.cells / 2 ||
		dom.bufCb.Width[1] >= b.cells / 2 ||
		dom.bufCb.Width[2] >= b.cells / 2) {
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
	seg := &segmentOverlap3D{ }
	dom := &domainOverlap3D{ }

	seg.boxWidth = b.cellWidth * float64(b.cells)
	dom.boxWidth = b.cellWidth * float64(b.cells)
	seg.cells = b.cells
	dom.cells = b.cells

	dom.domCb = b.cb
	seg.domCb = b.cb
	dom.bufCb = b.cb
	seg.bufCb = *hd.CellBounds(b.cells)

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
func NewBox(
	boxWidth float64, pts, cells int, q density.Quantity, config *io.BoxConfig,
) Box {
	if config.IsProjection() && q.CanProject() {
		return newBox2D(boxWidth, pts, cells, q, config)
	} else {
		return newBox3D(boxWidth, pts, cells, q, config)
	}
}

func newBox2D(
	boxWidth float64, pts, cells int, q density.Quantity, config *io.BoxConfig,
) Box {
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

	len := b.cb.Width[iDim] * b.cb.Width[jDim]

	g := geom.NewGridLocation(b.cb.Origin, b.cb.Width, boxWidth, cells)
	b.bvals = density.NewBuffer(q, len, b.pts, g)
	b.q = q

	return b
}

func newBox3D(
	boxWidth float64, pts, cells int, q density.Quantity, config *io.BoxConfig,
) Box {
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

	len := b.cb.Width[0] * b.cb.Width[1] * b.cb.Width[2]
	g := geom.NewGridLocation(b.cb.Origin, b.cb.Width, boxWidth, cells)
	b.bvals = density.NewBuffer(q, len, b.pts, g)

	return b
}

/////////////////////////////
// Overlap Implementations //
/////////////////////////////

type baseOverlap struct {
	bufCb, domCb geom.CellBounds
	cells int
	boxWidth float64
}

func (b *baseOverlap) BufferCellBounds() *geom.CellBounds { return &b.bufCb }
func (b *baseOverlap) DomainCellBounds() *geom.CellBounds { return &b.domCb }
func (w *baseOverlap) ScaleVecs(vs []geom.Vec) {
	w.bufCb.ScaleVecs(vs, w.cells, w.boxWidth)
}
func (b *baseOverlap) Cells() int { return b.cells }

type baseOverlap2D struct {
	baseOverlap
	proj int
}

func (w *baseOverlap2D) BufferSize() int {
	// I'm not doing this the obvious way (with a division) to avoid
	// integer overflow. (However unlikely that might be.)
	prod := 1
	for i := 0; i < 3; i++ {
		if i == w.proj { continue }
		prod *= w.bufCb.Width[i]
	}

	return prod
}


func (w *domainOverlap2D) Add(bbuf, bgrid density.Buffer) {
	bufNum, valid := bbuf.CountBuffer()
	gridNum, _ := bgrid.CountBuffer()

	if buf, ok := bbuf.ScalarBuffer(); ok {
		grid, _ := bgrid.ScalarBuffer()

		for i, val := range buf {
			grid[i] += val
			if valid { gridNum[i] += bufNum[i] }
		}
	} else if buf, ok := bbuf.VectorBuffer(); ok {
		grid, _ := bgrid.VectorBuffer()

		for i, val := range buf {
			grid[i][0] += val[0]
			grid[i][1] += val[1]
			grid[i][2] += val[2]
			if valid { gridNum[i] += bufNum[i] }
		}
	}
}

func (w *domainOverlap3D) Add(bbuf, bgrid density.Buffer) {
	bufNum, valid := bbuf.CountBuffer()
	gridNum, _ := bgrid.CountBuffer()

	if buf, ok := bbuf.ScalarBuffer(); ok {
		grid, _ := bgrid.ScalarBuffer()
		for i, val := range buf { 
			grid[i] += val
			if valid { gridNum[i] += bufNum[i] }
		}
	} else if buf, ok := bbuf.VectorBuffer(); ok {
		grid, _ := bgrid.VectorBuffer()
		for i, val := range buf {
			grid[i][0] += val[0]
			grid[i][1] += val[1]
			grid[i][2] += val[2]
			if valid { gridNum[i] += bufNum[i] }
		}
	}
}

type baseOverlap3D struct {
	baseOverlap
}

func (w *baseOverlap3D) BufferSize() int {
	return w.bufCb.Width[0] * w.bufCb.Width[1] * w.bufCb.Width[2]
}

func (w *segmentOverlap2D) Add(bbuf, bgrid density.Buffer) {
	bufNum, valid := bbuf.CountBuffer()
	gridNum, _ := bgrid.CountBuffer()

	if buf, ok := bbuf.ScalarBuffer(); ok {
		grid, _ := bgrid.ScalarBuffer()

		iDim, jDim := 0, 1
		if w.proj == 0 { iDim, jDim = 1, 2 }	
		if w.proj == 1 { iDim, jDim = 0, 2 }
		
		for jBuf := 0; jBuf < w.bufCb.Width[jDim]; jBuf++ {
			jDom := jBuf + (w.bufCb.Origin[jDim] - w.domCb.Origin[jDim])
			if jDom < 0 {
				jDom += w.cells
			} else if jDom >= w.cells {
				jDom -= w.cells
			}
			
			if jDom >= w.domCb.Width[jDim] { continue }
			flatDomJ := jDom * w.domCb.Width[iDim]
			flatBufJ := jBuf * w.bufCb.Width[iDim]
			
			for iBuf := 0; iBuf < w.bufCb.Width[iDim]; iBuf++ {
				iDom := iBuf + (w.bufCb.Origin[iDim] - w.domCb.Origin[iDim])
				if iDom < 0 {
					iDom += w.cells
				} else if iDom >= w.cells {
					iDom -= w.cells
				}
				if iDom >= w.domCb.Width[iDim] { continue }
				
				domIdx := iDom + flatDomJ
				bufIdx := iBuf + flatBufJ
				grid[domIdx] += buf[bufIdx]
				if valid { gridNum[domIdx] += bufNum[bufIdx] }
			}
		}
	} else if buf, ok := bbuf.VectorBuffer(); ok {
		grid, _ := bgrid.VectorBuffer()
		
		iDim, jDim := 0, 1
		if w.proj == 0 { iDim, jDim = 1, 2 }	
		if w.proj == 1 { iDim, jDim = 0, 2 }
		
		for jBuf := 0; jBuf < w.bufCb.Width[jDim]; jBuf++ {
			jDom := jBuf + (w.bufCb.Origin[jDim] - w.domCb.Origin[jDim])
			if jDom < 0 {
				jDom += w.cells
			} else if jDom >= w.cells {
				jDom -= w.cells
			}
			
			if jDom >= w.domCb.Width[jDim] { continue }
			flatDomJ := jDom * w.domCb.Width[iDim]
			flatBufJ := jBuf * w.bufCb.Width[iDim]
			
			for iBuf := 0; iBuf < w.bufCb.Width[iDim]; iBuf++ {
				iDom := iBuf + (w.bufCb.Origin[iDim] - w.domCb.Origin[iDim])
				if iDom < 0 {
					iDom += w.cells
				} else if iDom >= w.cells {
					iDom -= w.cells
				}
				if iDom >= w.domCb.Width[iDim] { continue }
				
				domIdx := iDom + flatDomJ
				bufIdx := iBuf + flatBufJ
				grid[domIdx][0] += buf[bufIdx][0]
				grid[domIdx][1] += buf[bufIdx][1]
				grid[domIdx][2] += buf[bufIdx][2]
				if valid { gridNum[domIdx] += bufNum[bufIdx] }
			}
		}
	}
}

func (w *segmentOverlap3D) Add(bbuf, bgrid density.Buffer) {
	bufNum, valid := bbuf.CountBuffer()
	gridNum, _ := bgrid.CountBuffer()

	if buf, ok := bbuf.ScalarBuffer(); ok {
		grid, _ := bgrid.ScalarBuffer()

		for zBuf := 0; zBuf < w.bufCb.Width[2]; zBuf++ {
			zDom := zBuf + (w.bufCb.Origin[2] - w.domCb.Origin[2])
			if zDom < 0 {
				zDom += w.cells
			} else if zDom >= w.cells {
				zDom -= w.cells
			}
			if zDom >= w.domCb.Width[2] { continue }
			flatDomZ := zDom * w.domCb.Width[1] * w.domCb.Width[0]
			flatBufZ := zBuf * w.bufCb.Width[1] * w.bufCb.Width[0]
			
			for yBuf := 0; yBuf < w.bufCb.Width[1]; yBuf++ {
				yDom := yBuf + (w.bufCb.Origin[1] - w.domCb.Origin[1])
				if yDom < 0 {
					yDom += w.cells
				} else if yDom >= w.cells {
					yDom -= w.cells
				}
				
				if yDom >= w.domCb.Width[1] { continue }
				flatDomY := yDom * w.domCb.Width[0]
				flatBufY := yBuf * w.bufCb.Width[0]
				
				for xBuf := 0; xBuf < w.bufCb.Width[0]; xBuf++ {
					xDom := xBuf + (w.bufCb.Origin[0] - w.domCb.Origin[0])
					if xDom < 0 {
					xDom += w.cells
					} else if xDom >= w.cells {
						xDom -= w.cells
					}
					if xDom >= w.domCb.Width[0] { continue }
					
					domIdx := xDom + flatDomY + flatDomZ
					bufIdx := xBuf + flatBufY + flatBufZ
					if valid { grid[domIdx] += buf[bufIdx] }
				}
			}
		}
	}  else if buf, ok := bbuf.VectorBuffer(); ok {
		grid, _ := bgrid.VectorBuffer()

		for zBuf := 0; zBuf < w.bufCb.Width[2]; zBuf++ {
			zDom := zBuf + (w.bufCb.Origin[2] - w.domCb.Origin[2])
			if zDom < 0 {
				zDom += w.cells
			} else if zDom >= w.cells {
				zDom -= w.cells
			}
			if zDom >= w.domCb.Width[2] { continue }
			flatDomZ := zDom * w.domCb.Width[1] * w.domCb.Width[0]
			flatBufZ := zBuf * w.bufCb.Width[1] * w.bufCb.Width[0]
			
			for yBuf := 0; yBuf < w.bufCb.Width[1]; yBuf++ {
				yDom := yBuf + (w.bufCb.Origin[1] - w.domCb.Origin[1])
				if yDom < 0 {
					yDom += w.cells
				} else if yDom >= w.cells {
					yDom -= w.cells
				}
				
				if yDom >= w.domCb.Width[1] { continue }
				flatDomY := yDom * w.domCb.Width[0]
				flatBufY := yBuf * w.bufCb.Width[0]
				
				for xBuf := 0; xBuf < w.bufCb.Width[0]; xBuf++ {
					xDom := xBuf + (w.bufCb.Origin[0] - w.domCb.Origin[0])
					if xDom < 0 {
					xDom += w.cells
					} else if xDom >= w.cells {
						xDom -= w.cells
					}
					if xDom >= w.domCb.Width[0] { continue }
					
					domIdx := xDom + flatDomY + flatDomZ
					bufIdx := xBuf + flatBufY + flatBufZ
					grid[domIdx][0] += buf[bufIdx][0]
					grid[domIdx][1] += buf[bufIdx][1]
					grid[domIdx][2] += buf[bufIdx][2]
					if valid { gridNum[domIdx] += bufNum[bufIdx] }
				}
			}
		}
	}
}

type segmentOverlap2D struct {
	baseOverlap2D
}

func (w *segmentOverlap2D) Interpolate(
	bbuf density.Buffer, xs, vs []geom.Vec,
	ptVal float64, bweights density.Buffer,
	low, high, jump int,
) {
	iDim, jDim, kDim := 0, 1, 2
	if w.proj == 0 { iDim, jDim, kDim = 1, 2, 0 }	
	if w.proj == 1 { iDim, jDim, kDim = 0, 2, 1 }
	kOffset := w.domCb.Origin[kDim] - w.bufCb.Origin[kDim]

	length := w.bufCb.Width[iDim]
	ptVal /= float64(w.domCb.Width[kDim])

	wlen := bweights.Length()
	buf, scalar := bbuf.ScalarBuffer()
	wbuf, _ := bweights.ScalarBuffer()
	vbuf, _ := bbuf.VectorBuffer()
	vwbuf, _ := bweights.VectorBuffer()
	counts, countValid := bbuf.CountBuffer()

	for idx := low; idx < high; idx += jump {
		pt := xs[idx]
		i, j := int(pt[iDim]), int(pt[jDim])
		
		k := intFloor(pt[kDim])
		kk := bound(k - kOffset, w.cells)
		
		if kk < w.domCb.Width[kDim] && kk >= 0 {
			bufIdx := i + j * length
			if !scalar {
				for dim := 0; dim < 3; dim++ {
					vbuf[bufIdx][dim] += ptVal * vwbuf[idx][dim]
				}
				if countValid { counts[bufIdx]++ }
			} else if wlen == 0 {
				buf[bufIdx] += ptVal
				if countValid { counts[bufIdx]++ }
			} else {
				buf[bufIdx] += ptVal * wbuf[idx]
				if countValid { counts[bufIdx]++ }
			}
		}
	}
}

type segmentOverlap3D struct {
	baseOverlap3D
}

func (w *segmentOverlap3D) Interpolate(
	bbuf density.Buffer, xs, vs []geom.Vec,
	ptVal float64, bweights density.Buffer,
	low, high, jump int,
) {
	length := w.bufCb.Width[0]
	area :=   w.bufCb.Width[0] * w.bufCb.Width[1]

	wlen := bweights.Length()
	buf, scalar := bbuf.ScalarBuffer()
	wbuf, _ := bweights.ScalarBuffer()
	vbuf, _ := bbuf.VectorBuffer()
	vwbuf, _ := bweights.VectorBuffer()
	counts, countValid := bbuf.CountBuffer()

	if !scalar {
		for idx := low; idx < high; idx += jump {
			pt := xs[idx]
			x, y, z := int(pt[0]), int(pt[1]), int(pt[2])
			bufIdx := x + y * length + z * area
			for dim := 0; dim < 3; dim++ {
				vbuf[bufIdx][dim] += ptVal * vwbuf[idx][dim]
			}
			if countValid { counts[bufIdx]++ }
		}
	} else if wlen == 0 {
		for idx := low; idx < high; idx += jump {
			pt := xs[idx]
			x, y, z := int(pt[0]), int(pt[1]), int(pt[2])
			bufIdx := x + y * length + z * area
			buf[bufIdx] += ptVal
			if countValid { counts[bufIdx]++ }
		}
	} else {
		for idx := low; idx < high; idx += jump {
			pt := xs[idx]
			x, y, z := int(pt[0]), int(pt[1]), int(pt[2])
			bufIdx := x + y * length + z * area
			buf[bufIdx] += ptVal * wbuf[idx]
			if countValid { counts[bufIdx]++ }
		}
	}
}

type domainOverlap2D struct {
	baseOverlap2D
}

func (w *domainOverlap2D) Interpolate(
	bbuf density.Buffer, xs, vs []geom.Vec,
	ptVal float64, bweights density.Buffer,
	low, high, jump int,
) {
	iDim, jDim, kDim := 0, 1, 2
	if w.proj == 0 {
		iDim, jDim, kDim = 1, 2, 0
	} else if w.proj ==1 {
		iDim, jDim, kDim = 0, 2, 1
	}

	length := w.bufCb.Width[iDim]
	ptVal /= float64(w.bufCb.Width[kDim])

	wlen := bweights.Length()
	buf, scalar := bbuf.ScalarBuffer()
	wbuf, _ := bweights.ScalarBuffer()
	vbuf, _ := bbuf.VectorBuffer()
	vwbuf, _ := bweights.VectorBuffer()
	counts, countValid := bbuf.CountBuffer()

	for idx := low; idx < high; idx += jump {
		pt := xs[idx]
		
		i := bound(intFloor(pt[iDim]), w.cells)
		if i < w.bufCb.Width[iDim] {
			j := bound(intFloor(pt[jDim]), w.cells)
			if j < w.bufCb.Width[jDim] {
				k := bound(intFloor(pt[kDim]), w.cells)
				if k < w.bufCb.Width[kDim] {
					bufIdx := i + j * length
					if !scalar {
						for dim := 0; dim < 3; dim++ {
							vbuf[bufIdx][dim] += ptVal * vwbuf[idx][dim]
						}
						if countValid { counts[bufIdx]++ }
					} else if wlen == 0 {
						buf[bufIdx] += ptVal
						if countValid { counts[bufIdx]++ }
					} else {
						buf[bufIdx] += ptVal * wbuf[idx]
						if countValid { counts[bufIdx]++ }
					}
				}
			}
		}
	}
}

type domainOverlap3D struct {
	baseOverlap3D
}

func (w *domainOverlap3D) Interpolate(
	bbuf density.Buffer, xs, vs []geom.Vec,
	ptVal float64, bweights density.Buffer,
	low, high, jump int,
) {
	length := w.bufCb.Width[0]
	area   := w.bufCb.Width[0] * w.bufCb.Width[1]

	wlen := bweights.Length()
	buf, scalar := bbuf.ScalarBuffer()
	wbuf, _ := bweights.ScalarBuffer()
	vbuf, _ := bbuf.VectorBuffer()
	vwbuf, _ := bweights.VectorBuffer()
	counts, countValid := bbuf.CountBuffer()

	for idx := low; idx < high; idx += jump {
		pt := xs[idx]
		
		i := bound(intFloor(pt[0]), w.cells)
		if i < w.bufCb.Width[0] {
			j := bound(intFloor(pt[1]), w.cells)
			if j < w.bufCb.Width[1] {
				k := bound(intFloor(pt[2]), w.cells)
				if k < w.bufCb.Width[2] {
					bufIdx := i + j * length + k * area
					if !scalar {	
						for dim := 0; dim < 3; dim++ {
							vbuf[bufIdx][dim] += (ptVal * vwbuf[idx][dim])
						}
						if countValid { counts[bufIdx]++ }
					} else if wlen == 0 {
						buf[bufIdx] += ptVal
						if countValid { counts[bufIdx]++ }
					} else {
						buf[bufIdx] += ptVal * wbuf[idx]
						if countValid { counts[bufIdx]++ }
					}
				}
			}
		}
	}
}

func intFloor(x float32) int {
	if x > 0 {
		return int(x)
	} else {
		return int(x - 1)
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
