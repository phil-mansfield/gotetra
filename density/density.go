/*package density interpolates sequences of particle positions onto a density
grid.
*/
package density

import (
	"fmt"

	"github.com/phil-mansfield/gotetra/rand"
	"github.com/phil-mansfield/gotetra/geom"
)

type Quantity int64
const (
	Density Quantity = iota
	DensityGradient
	Velocity
	VelocityDivergence
	VelocityCurl
	EndQuantity
)

func (q Quantity) String() string {
	if q <= 0 || q >= EndQuantity {
		panic(fmt.Sprintf("Value %d out of range for Quantity type.", q))
	}

	switch q {
	case Density:
		return "Density"
	case DensityGradient:
		return "DensityGradient"
	case Velocity:
		return "Velocity"
	case VelocityDivergence:
		return "VelocityDivergence"
	case VelocityCurl:
		return "VelocityCurl"
	}

	panic("Quantity.String() missing a switch clause.")
}

func QuantityFromString(str string) (q Quantity, ok bool) {
	switch str {
	case "Density":
		return Density, true
	case "DensityGradient":
		return DensityGradient, true
	case "Velocity":
		return Velocity, true
	case "VelocityDivergence":
		return VelocityDivergence, true
	case "VelocityCurl":
		return VelocityCurl, true
	}
	return 0, false
}

func (q Quantity) CanProject() bool {
	if q <= 0 || q >= EndQuantity {
		panic(fmt.Sprintf("Value %d out of range for Quantity type.", q))
	}

	switch q {
	case Density, Velocity:
		return true
	case DensityGradient, VelocityCurl, VelocityDivergence:
		return false
	}

	panic("Quantity.String() missing a switch clause.")
}

type Interpolator interface {
	Interpolate(
		buf Buffer, xs []geom.Vec,
		ptVal float64, weights Buffer,
		low, high, jump int,
	)

	// The bounding box around the box being written to.
	DomainCellBounds() *geom.CellBounds
	// The bounding box around the box being written from.
	BufferCellBounds() *geom.CellBounds

	Cells() int
}

type Buffer interface {
	Slice(low, high int)
	Length() int
	Quantity() Quantity
	Clear()

	SetGridLocation(g *geom.GridLocation)

	ScalarBuffer() (vals []float64, ok bool)
	VectorBuffer() (vals [][3]float64, ok bool)
	FinalizedScalarBuffer() (vals []float32, ok bool)
	FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool)
}

var NilBuffer = &scalarBuffer{ []float64{} }

type scalarBuffer struct { vals []float64 }
type vectorBuffer struct { vecs [][3]float64 }
type densityBuffer struct { scalarBuffer }
type gradientBuffer struct {
	scalarBuffer
	g *geom.GridLocation
}
type velocityBuffer struct {
	vectorBuffer
	weights *vectorBuffer
}
type divergenceBuffer struct {
	vectorBuffer
	weights *vectorBuffer
	g *geom.GridLocation
}
type curlBuffer struct {
	vectorBuffer
	weights *vectorBuffer
	g *geom.GridLocation
}

func (b *scalarBuffer) SetGridLocation(g *geom.GridLocation) { }
func (b *vectorBuffer) SetGridLocation(g *geom.GridLocation) { }
func (b *gradientBuffer) SetGridLocation(g *geom.GridLocation) { b.g = g }
func (b *curlBuffer) SetGridLocation(g *geom.GridLocation) { b.g = g }
func (b *divergenceBuffer) SetGridLocation(g *geom.GridLocation) { b.g = g }

func (b *scalarBuffer) Quantity() Quantity {
	panic("Qunatity() called on raw scalarBuffer.")
}
func (b *vectorBuffer) Quantity() Quantity {
	panic("Qunatity() called on raw vectorBuffer.")
}
func (b *densityBuffer) Quantity() Quantity { return Density }
func (b *gradientBuffer) Quantity() Quantity { return DensityGradient }
func (b *velocityBuffer) Quantity() Quantity { return Velocity }
func (b *divergenceBuffer) Quantity() Quantity { return VelocityDivergence }
func (b *curlBuffer) Quantity() Quantity { return VelocityCurl }

func (buf *vectorBuffer) Slice(low, high int) {
	buf.vecs = buf.vecs[low: high]
}

func (buf *velocityBuffer) Slice(low, high int) {
	buf.vectorBuffer.Slice(low, high)
	buf.weights.Slice(low, high)
}

func (buf *divergenceBuffer) Slice(low, high int) {
	buf.vectorBuffer.Slice(low, high)
	buf.weights.Slice(low, high)
}

func (buf *curlBuffer) Slice(low, high int) {
	buf.vectorBuffer.Slice(low, high)
	buf.weights.Slice(low, high)
}

func (buf *vectorBuffer) Length() int {
	return len(buf.vecs)
}

func (buf *vectorBuffer) Clear() {
	for i := range buf.vecs {
		buf.vecs[i][0], buf.vecs[i][1], buf.vecs[i][2] = 0, 0, 0
	}
}

func (buf *vectorBuffer) ScalarBuffer() (vals []float64, ok bool) {
	return nil, false
}

func (buf *vectorBuffer) VectorBuffer() (vals [][3]float64, ok bool) {
	return buf.vecs, true
}

func (buf *curlBuffer) FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool) {
	xs, oxs := make([]float32, len(buf.vecs)), make([]float32, len(buf.vecs))
	ys, oys := make([]float32, len(buf.vecs)), make([]float32, len(buf.vecs))
	zs, ozs := make([]float32, len(buf.vecs)), make([]float32, len(buf.vecs))
	for i, vec := range buf.vecs {
		xs[i], ys[i], zs[i] = float32(vec[0]), float32(vec[1]), float32(vec[2])
	}
	vecs := [3][]float32{ xs, ys, zs }
	out := [3][]float32{ oxs, oys, ozs }
	buf.g.Curl(vecs, out, &geom.DerivOptions{ true, geom.None, 4})
	return oxs, oys, ozs, true
}

func (buf *vectorBuffer) FinalizedScalarBuffer() (vals []float32, ok bool) {
	return nil, false
}

func (buf *divergenceBuffer) FinalizedScalarBuffer() (vals []float32, ok bool) {
	out := make([]float32, len(buf.vecs))
	vecs := [3][]float32 {
		make([]float32, len(buf.vecs)),
		make([]float32, len(buf.vecs)),
		make([]float32, len(buf.vecs)),
	}

	buf.g.Divergence(vecs, out, &geom.DerivOptions{ true, geom.None, 4 })

	return out, true
}

func (buf *divergenceBuffer) FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool) {
	return nil, nil, nil, false
}

func (buf *vectorBuffer) FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool) {
	xs = make([]float32, len(buf.vecs))
	ys = make([]float32, len(buf.vecs))
	zs = make([]float32, len(buf.vecs))
	for i, vec := range buf.vecs {
		xs[i], ys[i], zs[i] = float32(vec[0]), float32(vec[1]), float32(vec[2])
	}
	return xs, ys, zs, true
}

func (buf *scalarBuffer) Slice(low, high int) {
	buf.vals = buf.vals[low: high]
}

func (buf *scalarBuffer) Length() int {
	return len(buf.vals)
}

func (buf *scalarBuffer) Clear() {
	for i := range buf.vals { buf.vals[i] = 0 }
}

func (buf *scalarBuffer) ScalarBuffer() (vals []float64, ok bool) {
	return buf.vals, true
}

func (buf *scalarBuffer) VectorBuffer() (vals [][3]float64, ok bool) {
	return nil, false
}

func (buf *scalarBuffer) FinalizedScalarBuffer() (vals []float32, ok bool) {
	vals32 := make([]float32, len(buf.vals))
	for i, x := range buf.vals { vals32[i] = float32(x) }
	return vals32, true
}

func (buf *scalarBuffer) FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool) {
	return nil, nil, nil, false
}

func (buf *gradientBuffer) FinalizedScalarBuffer() (vals []float32, ok bool) {
	return nil, false
}

func (buf *gradientBuffer) FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool) {
	vals := make([]float32, len(buf.vals))
	for i, x := range buf.vals { vals[i] = float32(x) }
	out := [3][]float32 {
		make([]float32, len(buf.vals)),
		make([]float32, len(buf.vals)),
		make([]float32, len(buf.vals)),
	}

	buf.g.Gradient(vals, out, &geom.DerivOptions{ true, geom.None, 4 })
	return out[0], out[1], out[2], true
}


func NewBuffer(q Quantity, len int, g *geom.GridLocation) Buffer {
	switch q {
	case Density:
		return &densityBuffer{
			scalarBuffer{ make([]float64, len) },
		}
	case DensityGradient:
		return &gradientBuffer{
			scalarBuffer{ make([]float64, len) }, g,
		}
	case Velocity:
		return &velocityBuffer{
			vectorBuffer{ make([][3]float64, len) },
			&vectorBuffer{ make([][3]float64, len) },
		}
	case VelocityDivergence:
		return &divergenceBuffer{
			vectorBuffer{ make([][3]float64, len) },
			&vectorBuffer{ make([][3]float64, len) },
			g,
		}
	case VelocityCurl:
		return &curlBuffer{
			vectorBuffer{ make([][3]float64, len) },
			&vectorBuffer{ make([][3]float64, len) },
			g,
		}
	default:
		panic(fmt.Sprintf("Unrecognized Quantity %v", q))
	}
	panic(":3")
}

type ngp struct { }

// The ordering of these fields makes no goddamned sense.
type mcarlo struct {
	subIntr Interpolator
	segWidth int64
	points, cells int

	gen *rand.Generator

	skip int64
	// Buffers
	idxBuf geom.TetraIdxs
	tet geom.Tetra

	unitBufs [][]geom.Vec
	vecBuf []geom.Vec
}

func MonteCarlo(
	segWidth int64,
	points, cells int,
	skip int64,
	unitBufs [][]geom.Vec,
	subIntr Interpolator,
) Interpolator {
	mc := &mcarlo{
		subIntr, segWidth, points, cells,
		rand.NewTimeSeed(rand.Xorshift), skip,
		geom.TetraIdxs{}, geom.Tetra{},
		unitBufs,
		make([]geom.Vec, points),
	}

	return mc
}

func (intr *mcarlo) DomainCellBounds() *geom.CellBounds {
	return intr.subIntr.DomainCellBounds()
}

func (intr *mcarlo) BufferCellBounds() *geom.CellBounds {
	return intr.subIntr.BufferCellBounds()
}

func (intr *mcarlo) Cells() int {
	return intr.subIntr.Cells()
}

func (intr *mcarlo) Interpolate(
	buf Buffer, xs []geom.Vec,
	ptVal float64, weights Buffer,
	low, high, jump int,
) {
	segWidth := intr.segWidth
	gridWidth := segWidth + 1
	idxWidth := intr.segWidth / intr.skip

	ptVal = ptVal / float64(intr.points) / 6.0 *
		float64(intr.skip * intr.skip * intr.skip)

	tetCb := &geom.CellBounds{}

	relCb := &geom.CellBounds{}
	relCb.Width = intr.DomainCellBounds().Width
	relCb.Origin[0], relCb.Origin[1], relCb.Origin[2] =
		cbSubtr(intr.DomainCellBounds(), intr.BufferCellBounds())

	jump64 := int64(jump)

	// I hate this so much. I hate it so much:
	//
	// (Also, I have no idea how it works.)
	neverMod := *intr.subIntr.DomainCellBounds() !=
		*intr.subIntr.BufferCellBounds() ||
		(intr.Cells() / 2 > relCb.Width[0] &&
		intr.Cells() / 2 > relCb.Width[1] &&
		intr.Cells() / 2 > relCb.Width[2])
	mods := 0

	for idx := int64(low); idx < int64(high); idx += jump64 {
		x, y, z := coords(idx, idxWidth)
		gridIdx := index(x * intr.skip, y * intr.skip, z * intr.skip, gridWidth)

		for dir := 0; dir < 6; dir++ {
			intr.idxBuf.Init(gridIdx, gridWidth, intr.skip, dir)
			
			intr.tet.Init(
				&xs[intr.idxBuf[0]],
				&xs[intr.idxBuf[1]],
				&xs[intr.idxBuf[2]],
				&xs[intr.idxBuf[3]],
			)
			
			intr.tet.CellBoundsAt(1.0, tetCb)

			if !tetCb.Intersect(relCb, intr.cells) { continue }

			bufIdx := intr.gen.UniformInt(0, len(intr.unitBufs))
			intr.tet.DistributeTetra(
				intr.unitBufs[bufIdx],
				intr.vecBuf,
			)

			// Lol, whatever.
			if !neverMod {
				for j := 0; j < 3; j++ {
					if coordNeg(&intr.tet, j) {
						mods++
						modCoord(intr.vecBuf, j, float32(intr.Cells()))
					}
				}
			}

			var bweights Buffer
			switch b := buf.(type) {
			case *densityBuffer, *gradientBuffer:
				bweights = NilBuffer
			case *velocityBuffer:
				wbuf, ok := b.weights.VectorBuffer()
				if !ok { panic("buf is non-vector when vector is required.") }
				intr.tet.DistributeTetra64(intr.unitBufs[bufIdx], wbuf)
				bweights = b.weights
			case *curlBuffer:
				wbuf, ok := b.weights.VectorBuffer()
				if !ok { panic("buf is non-vector when vector is required.") }
				intr.tet.DistributeTetra64(intr.unitBufs[bufIdx], wbuf)
				bweights = b.weights
			case *divergenceBuffer:
				wbuf, ok := b.weights.VectorBuffer()
				if !ok { panic("buf is non-vector when vector is required.") }
				intr.tet.DistributeTetra64(intr.unitBufs[bufIdx], wbuf)
				bweights = b.weights
			}

			intr.subIntr.Interpolate(
				buf, intr.vecBuf, ptVal, bweights, 0, intr.points, 1,
			)
		}
	}
}

func coordNeg(tet *geom.Tetra, j int) bool {
	for i := 0; i < 4; i++ {
		if tet.Corners[i][j] < 0 { return true }
	}
	return false
}

func modCoord(buf []geom.Vec, j int, width float32) {
	for i := range buf {
		if buf[i][j] < 0 { buf[i][j] += width }
	}
}

func index(x, y, z, cells int64) int64 {
	return x + y * cells + z * cells * cells
}

func coords(idx, cells int64) (x, y, z int64) {
	x = idx % cells
	y = (idx % (cells * cells)) / cells
	z = idx / (cells * cells)
	return x, y, z
}

func cbSubtr(cb1, cb2 *geom.CellBounds) (i, j, k int) {
	i = cb1.Origin[0] - cb2.Origin[0]
	j = cb1.Origin[1] - cb2.Origin[1]
	k = cb1.Origin[2] - cb2.Origin[2]
	return i, j, k
}
