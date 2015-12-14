package interpolate

import (
	"fmt"
)

////////////////////////////
// BiCubic Implementation //
////////////////////////////

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

/////////////////////////////
// TriCubic Implementation //
/////////////////////////////

type TriCubic struct {
	xs, ys, zs []float64
	vals []float64
	nx, ny int

	lastY, lastZ float64
	zSplines []*Spline
	ySplineVals [][]float64
	ySplines []*Spline
	xSplineVals []float64
	xSpline *Spline
}


func NewTriCubic(xs, ys, zs, vals []float64) *TriCubic {
	if len(xs)*len(ys)*len(zs) != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but len(xs) = %d, len(ys) = %d, and len(zs) = %d",
			len(vals), len(xs), len(ys), len(zs),
		))
	}

	tri := &TriCubic{}
	tri.nx = len(xs)
	tri.ny = len(ys)
	tri.vals = vals

	tri.xs, tri.ys, tri.zs = xs, ys, zs

	tri.initSplines()

	return tri
}

func NewUniformTriCubic(
	x0, dx float64, nx int,
	y0, dy float64, ny int,
	z0, dz float64, nz int,
	vals []float64,
) *TriCubic {
	if nx*ny*nz != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but len(xs) = %d, len(ys) = %d, len(zs) = %d",
			len(vals), nx, ny, nz,
		))
	}

	tri := &TriCubic{}
	tri.nx = nx
	tri.ny = ny
	tri.vals = vals

	tri.xs = make([]float64, nx)
	tri.ys = make([]float64, ny)
	tri.ys = make([]float64, nz)
	for i := range tri.xs { tri.xs[i] = x0 + float64(i)*dx }
	for i := range tri.ys { tri.ys[i] = y0 + float64(i)*dy }
	for i := range tri.zs { tri.zs[i] = z0 + float64(i)*dz }

	tri.initSplines()

	return tri
}

func (tri *TriCubic) initSplines() {
	tri.zSplines = make([]*Spline, len(tri.xs) * len(tri.ys))
	for xi := range tri.xs {
		for yi := range tri.ys {
			zVals := make([]float64, len(tri.zs))
			for zi := range tri.zs {
				zVals[yi] = tri.vals[tri.nx*tri.ny*zi + tri.nx*yi + xi]
			}

			tri.zSplines[yi*tri.nx + xi] = NewSpline(tri.zs, zVals)
		}
	}

	tri.lastZ = tri.zs[0]

	tri.ySplineVals = make([][]float64, len(tri.xs))
	for xi := range tri.xs {
		tri.ySplineVals[xi] = make([]float64, len(tri.ys))
		for yi := range tri.ys {
			tri.ySplineVals[xi][yi] =
				tri.zSplines[yi*tri.nx + xi].Eval(tri.lastZ)
		}
		tri.ySplines[xi] = NewSpline(tri.xs, tri.ySplineVals[xi])
	}

	tri.lastY = tri.ys[0]

	tri.xSplineVals = make([]float64, len(tri.xs))
	for xi := range tri.xSplineVals {
		tri.xSplineVals[xi] = tri.ySplines[xi].Eval(tri.lastY)
	}

	tri.xSpline = NewSpline(tri.xs, tri.xSplineVals)
}

func (tri *TriCubic) Eval(x, y, z float64) float64 {
	if y != tri.lastY || z != tri.lastZ {
		if z != tri.lastZ {
			tri.lastZ = z
			
			tri.ySplineVals = make([][]float64, len(tri.xs))
			for xi := range tri.xs {
				tri.ySplineVals[xi] = make([]float64, len(tri.ys))
				for yi := range tri.ys {
					tri.ySplineVals[xi][yi] =
						tri.zSplines[yi*tri.nx + xi].Eval(tri.lastZ)
				}
				tri.ySplines[xi].Init(tri.xs, tri.ySplineVals[xi])
			}			
		}

		tri.lastY = y
		
		tri.xSplineVals = make([]float64, len(tri.xs))
		for xi := range tri.xSplineVals {
			tri.xSplineVals[xi] = tri.ySplines[xi].Eval(tri.lastY)
		}
		
		tri.xSpline.Init(tri.xs, tri.xSplineVals)

	}

	return tri.xSpline.Eval(x)
}

func (tri *TriCubic) EvalAll(xs, ys, zs []float64, out ...[]float64) []float64 {
	if len(out) == 0 { out = [][]float64{ make([]float64, len(xs)) } }
	for i := range xs { out[0][i] = tri.Eval(xs[i], ys[i], zs[i]) }
	return out[0]
}
