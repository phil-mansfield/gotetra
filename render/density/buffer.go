package density

import (
	"fmt"

	"github.com/phil-mansfield/gotetra/geom"
)

// This is a giant clusterfuck.
type Buffer interface {
	// Array Management
	Slice(low, high int)
	Length() int
	Clear()

	// Getters and Setters
	Quantity() Quantity
	SetGridLocation(g *geom.GridLocation)
	SetVectors(vecs []geom.Vec) bool

	// Buffer Retrieval
	CountBuffer() (num []int, ok bool)
	ScalarBuffer() (vals []float64, ok bool)
	VectorBuffer() (vals [][3]float64, ok bool)
	FinalizedScalarBuffer() (vals []float32, ok bool)
	FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool)
}

var NilBuffer = &scalarBuffer{ []float64{} }

func NewBuffer(q Quantity, len, wlen int, g *geom.GridLocation) Buffer {
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
			&vectorBuffer{ make([][3]float64, wlen) },
			make([]int, len),
		}
	case VelocityDivergence:
		return &divergenceBuffer{
			vectorBuffer{ make([][3]float64, len) },
			&vectorBuffer{ make([][3]float64, wlen) },
			make([]int, len),
			g,
		}
	case VelocityCurl:
		return &curlBuffer{
			vectorBuffer{ make([][3]float64, len) },
			&vectorBuffer{ make([][3]float64, wlen) },
			make([]int, len),
			g,
		}
	default:
		panic(fmt.Sprintf("Unrecognized Quantity %v", q))
	}
	panic(":3")
}


/////////////////////////////////
// scalarBuffer implementation //
/////////////////////////////////

type scalarBuffer struct { vals []float64 }

// Array Manipulation //

func (buf *scalarBuffer) Slice(low, high int) {
	buf.vals = buf.vals[low: high]
}

func (buf *scalarBuffer) Length() int {
	return len(buf.vals)
}

func (buf *scalarBuffer) Clear() {
	for i := range buf.vals { buf.vals[i] = 0 }
}

// Getters and Setters //

func (b *scalarBuffer) Quantity() Quantity {
	panic("Qunatity() called on raw scalarBuffer.")
}

func (b *scalarBuffer) SetGridLocation(g *geom.GridLocation) { }

func (b *scalarBuffer) SetVectors(vecs []geom.Vec) bool { return false }


// Buffer Retreival //

func (buf *scalarBuffer) CountBuffer() (num []int, ok bool) {
	return nil, false
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

/////////////////////////////////
// vectorBuffer implementation //
/////////////////////////////////

type vectorBuffer struct { vecs [][3]float64 }

// Array Manipulation //

func (buf *vectorBuffer) Slice(low, high int) {
	buf.vecs = buf.vecs[low: high]
}

func (buf *vectorBuffer) Length() int {
	return len(buf.vecs)
}

func (buf *vectorBuffer) Clear() {
	for i := range buf.vecs {
		buf.vecs[i][0], buf.vecs[i][1], buf.vecs[i][2] = 0, 0, 0
	}
}

// Getters and Setters //

func (b *vectorBuffer) Quantity() Quantity {
	panic("Qunatity() called on raw vectorBuffer.")
}

func (b *vectorBuffer) SetGridLocation(g *geom.GridLocation) { }

func (b *vectorBuffer) SetVectors(vecs []geom.Vec) bool {
	for i := range vecs {
		for j := 0; j < 3; j++ { b.vecs[i][j] = float64(vecs[i][j]) }
	}
	return true
}

// Buffer Retrieval //

func (buf *vectorBuffer) CountBuffer() (num []int, ok bool) {
	return nil, false
}

func (buf *vectorBuffer) ScalarBuffer() (vals []float64, ok bool) {
	return nil, false
}

func (buf *vectorBuffer) VectorBuffer() (vals [][3]float64, ok bool) {
	return buf.vecs, true
}

func (buf *vectorBuffer) FinalizedScalarBuffer() (vals []float32, ok bool) {
	return nil, false
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

//////////////////////////////////
// densityBuffer implementation //
//////////////////////////////////

type densityBuffer struct { scalarBuffer }

// Array Manipulation //

// scalarBuffer.Length

// Getters and Setters //

func (b *densityBuffer) Quantity() Quantity { return Density }

// scalarBuffer.SetGridLocation

// scalarBuffer.SetVectors

// Buffer Retrieval //

///////////////////////////////////
// gradientBuffer implementation //
///////////////////////////////////

type gradientBuffer struct {
	scalarBuffer
	g *geom.GridLocation
}

// Array Manipulation //

// vectorBuffer.Length

// Getters and Setters //

func (b *gradientBuffer) Quantity() Quantity { return DensityGradient }

func (b *gradientBuffer) SetGridLocation(g *geom.GridLocation) { b.g = g }

// vectorBuffer.SetGridLocation

// vectorBuffer.SetVectors

// Buffer Retrieval //

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


/////////////////////////////////
// velocityBuffer implemtation //
/////////////////////////////////

type velocityBuffer struct {
	vectorBuffer
	weights *vectorBuffer
	num []int
}

// Array Manipulation //

func (buf *velocityBuffer) Slice(low, high int) {
	buf.vectorBuffer.Slice(low, high)
	buf.num = buf.num[low: high]
}

func (buf *velocityBuffer) Clear() {
	for i := range buf.vecs {
		buf.vecs[i][0], buf.vecs[i][1], buf.vecs[i][2] = 0, 0, 0
		buf.num[i] = 0
	}
}

// vectorBuffer.Length

// Getters and Setters //

func (b *velocityBuffer) Quantity() Quantity { return Velocity }

// vectorBuffer.SetGridLocation

// vectorBuffer.SetVectors

// Buffer Retrieval //

func (buf *velocityBuffer) CountBuffer() (num []int, ok bool) {
	return buf.num, true
}

func (buf *velocityBuffer) FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool) {
	xs = make([]float32, len(buf.vecs))
	ys = make([]float32, len(buf.vecs))
	zs = make([]float32, len(buf.vecs))
	for i, vec := range buf.vecs {
		if i % (len(buf.vecs) / 30) == 0 {
			fmt.Println(i, buf.num[i], vec[0], vec[1], vec[2])
		}
		n := float32(buf.num[i])
		if buf.num[i] == 0 { continue }
		xs[i], ys[i], zs[i] = float32(vec[0])/n, float32(vec[1])/n, float32(vec[2])/n
	}
	return xs, ys, zs, true
}

///////////////////////////////////
// divergenceBuffer implemtation //
///////////////////////////////////

type divergenceBuffer struct {
	vectorBuffer
	weights *vectorBuffer
	num []int
	g *geom.GridLocation
}

// Array Manipulation //

func (buf *divergenceBuffer) Slice(low, high int) {
	buf.vectorBuffer.Slice(low, high)
	buf.num = buf.num[low: high]
}

func (buf *divergenceBuffer) Clear() {
	for i := range buf.vecs {
		buf.vecs[i][0], buf.vecs[i][1], buf.vecs[i][2] = 0, 0, 0
		buf.num[i] = 0
	}
}

// scalarBuffer.Length

// Getters and Setters //

func (b *divergenceBuffer) Quantity() Quantity { return VelocityDivergence }

func (b *divergenceBuffer) SetGridLocation(g *geom.GridLocation) { b.g = g }

// scalarBuffer.SetVectors

// Buffer Retrieval //

func (buf *divergenceBuffer) FinalizedScalarBuffer() (vals []float32, ok bool) {
	out := make([]float32, len(buf.vecs))
	vecs := [3][]float32 {
		make([]float32, len(buf.vecs)),
		make([]float32, len(buf.vecs)),
		make([]float32, len(buf.vecs)),
	}

	for i, vec := range buf.vecs {
		n := float32(buf.num[i])
		if buf.num[i] == 0 { continue }
		vecs[0][i] = float32(vec[0])/n
		vecs[1][i] = float32(vec[0])/n
		vecs[2][i] = float32(vec[2])/n
	} 

	buf.g.Divergence(vecs, out, &geom.DerivOptions{ true, geom.None, 4 })

	return out, true
}

func (buf *divergenceBuffer) FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool) {
	return nil, nil, nil, false
}

func (buf *divergenceBuffer) CountBuffer() (num []int, ok bool) {
	return buf.num, true
}

/////////////////////////////
// curlBuffer implemtation //
/////////////////////////////

type curlBuffer struct {
	vectorBuffer
	weights *vectorBuffer
	num []int
	g *geom.GridLocation
}

// Array Manipulation //

func (buf *curlBuffer) Slice(low, high int) {
	buf.vectorBuffer.Slice(low, high)
	buf.num = buf.num[low: high]
}

func (buf *curlBuffer) Clear() {
	for i := range buf.vecs {
		buf.vecs[i][0], buf.vecs[i][1], buf.vecs[i][2] = 0, 0, 0
		buf.num[i] = 0
	}
}

// vectorBuffer.Length

// Getters and Setters //

func (b *curlBuffer) Quantity() Quantity { return VelocityCurl }

func (b *curlBuffer) SetGridLocation(g *geom.GridLocation) { b.g = g }

// vectorBuffer.SetVectors

// Buffer Retrieval //

func (buf *curlBuffer) CountBuffer() (num []int, ok bool) {
	return buf.num, true
}

func (buf *curlBuffer) FinalizedVectorBuffer() (xs, ys, zs []float32, ok bool) {
	xs, oxs := make([]float32, len(buf.vecs)), make([]float32, len(buf.vecs))
	ys, oys := make([]float32, len(buf.vecs)), make([]float32, len(buf.vecs))
	zs, ozs := make([]float32, len(buf.vecs)), make([]float32, len(buf.vecs))
	for i, vec := range buf.vecs {
		n := float32(buf.num[i])
		if buf.num[i] == 0 { continue }
		xs[i], ys[i], zs[i] = float32(vec[0])/n, float32(vec[1])/n, float32(vec[2])/n
	}
	vecs := [3][]float32{ xs, ys, zs }
	out := [3][]float32{ oxs, oys, ozs }
	buf.g.Curl(vecs, out, &geom.DerivOptions{ true, geom.None, 4})
	return oxs, oys, ozs, true
}
