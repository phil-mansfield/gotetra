package interpolate

import (
	"fmt"
)

type BiCubic struct {
	xs, ys []float64
	vals []float64
	nx int

	lastY float64
	ySplines []*Spline
	xSplineVals []float64
	xSpline *Spline
}

func NewBiCubic(xs, ys, vals []float64) *BiCubic {
	if len(xs) * len(ys) != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but len(xs) = %d and len(ys) = %d",
			len(vals), len(xs), len(ys),
		))
	}

	bi := &BiCubic{}
	bi.nx = len(xs)
	bi.vals = vals

	bi.xs, bi.ys = xs, ys

	bi.initSplines()

	return bi
}

func NewUniformBiCubic(
	x0, dx float64, nx int,
	y0, dy float64, ny int,
	vals []float64,
) *BiCubic {
	if nx*ny != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but len(xs) = %d and len(ys) = %d",
			len(vals), nx, ny,
		))
	}

	bi := &BiCubic{}
	bi.nx = nx
	bi.vals = vals

	bi.xs = make([]float64, nx)
	bi.ys = make([]float64, ny)
	for i := range bi.xs { bi.xs[i] = x0 + float64(i)*dx }
	for i := range bi.ys { bi.ys[i] = y0 + float64(i)*dy }

	bi.initSplines()

	return bi
}

func (bi *BiCubic) initSplines() {
	bi.ySplines = make([]*Spline, len(bi.xs))

	for xi := range bi.xs {
		yVals := make([]float64, len(bi.ys))
		for yi := range bi.ys {
			yVals[yi] = bi.vals[bi.nx * yi + xi]
		}

		bi.ySplines[xi] = NewSpline(bi.ys, yVals)
	}

	bi.lastY = bi.ys[0]
	bi.xSplineVals = make([]float64, len(bi.xs))
	for i := range bi.xSplineVals {
		bi.xSplineVals[i] = bi.ySplines[i].Eval(bi.lastY)
	}

	bi.xSpline = NewSpline(bi.xs, bi.xSplineVals)
}

func (bi *BiCubic) Eval(x, y float64) float64 {
	if y != bi.lastY {
		bi.lastY = y
		for i := range bi.xSplineVals {
			bi.xSplineVals[i] = bi.ySplines[i].Eval(y)
		}

		bi.xSpline.Init(bi.xs, bi.xSplineVals)
	}

	return bi.xSpline.Eval(x)
}

func (bi *BiCubic) EvalAll(xs, ys []float64, out ...[]float64) []float64 {
	if len(out) == 0 { out = [][]float64{ make([]float64, len(xs)) } }
	for i := range xs { out[0][i] = bi.Eval(xs[i], ys[i]) }
	return out[0]
}

type TriCubic struct {
}


func NewTriCubic(xs, ys, zs, vals []float64) *TriCubic {
	panic("NYI")
}

func NewUniformTriCubic(
	x0, dx float64, nx int,
	y0, dy float64, ny int,
	z0, dz float64, nz int,
	vals []float64,
) *TriCubic {
	panic("NYI")
}

func (tri *TriCubic) Eval(x, y, z float64) float64 {
	panic("NYI")
}

func (tri *TriCubic) EvalAll(xs, ys, zs []float64, out ...[]float64) []float64 {
	if len(out) == 0 { out = [][]float64{ make([]float64, len(xs)) } }
	for i := range xs { out[0][i] = tri.Eval(xs[i], ys[i], zs[i]) }
	return out[0]
}
