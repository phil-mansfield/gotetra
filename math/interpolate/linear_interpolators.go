package interpolate

import (
	"fmt"
)

///////////////////////////
// Linear Implementation //
///////////////////////////

// Linear is a linear interpolator.
type Linear struct {
	xs searcher
	vals []float64
}

// NewLinear creates a linear interpolator for a sequence of strictly increasing
// or strictly decreasing point, xs, which take on the values given by vals.
//
// Lookups will occur in O(log |xs|), possibly faster depending on the access
// pattern and data layout.
func NewLinear(xs, vals []float64) *Linear {
	if len(xs) != len(vals) {
		panic("Length of input slices are not equal.")
	}
	lin := &Linear{}
	lin.xs.init(xs)
	lin.vals = vals
	return lin
}

// NewUniformLinear creates a linear interplator where a uniformly spaced
// sequence of x values starting at x0 and separated by dx and whose values are
// given by vals.
//
// Lookups will be O(1).
func NewUniformLinear(x0, dx float64, vals []float64) *Linear {
	lin := &Linear{}
	lin.xs.unifInit(x0, dx, len(vals))
	lin.vals = vals
	return lin
}

// Eval returns the interpolated value at x.
//
// Eval panics if called on a values outside the supplied range on inputs.
func (lin *Linear) Eval(x float64) float64 {
	i1 := lin.xs.search(x)
	i2 := i1 + 1
	x1, x2 := lin.xs.val(i1), lin.xs.val(i2)
	v1, v2 := lin.vals[i1], lin.vals[i2]

	return ((v2 - v1) / (x2 - x1)) * (x - x1) + v1
}

// EvalAll evaluates the interpolator at all the given x values. If an output
// array is given, the output is written to that array (the array is still
// returned as a convenience).
//
// If more than one output array is provided, only the first is used.
func (lin *Linear) EvalAll(xs []float64, out ...[]float64) []float64 {
	if len(out) == 0 { out = [][]float64{ make([]float64, len(xs)) } }
	for i, x := range xs { out[0][i] = lin.Eval(x) }
	return out[0]
}


/////////////////////////////
// BiLinear Implementation //
/////////////////////////////

// BiLinear is a bi-linear interpolator.
type BiLinear struct {
	xs, ys searcher
	vals []float64
	nx int
}

// NewBiLinear creates a bi-linear interpolator on top of a grid with the
// values given by vals. The values of the x and y grid lines are given by
// xs ans ys. The vals grid is indexed in the usual way:
// vals(ix, iy) -> vals[ix + iy*nx].
//
// Panics if len(xs) * len(ys) != len(vals).
func NewBiLinear(xs, ys, vals []float64) *BiLinear {
	bi := &BiLinear{}
	bi.xs.init(xs)
	bi.ys.init(ys)
	bi.nx = len(xs)
	bi.vals = vals

	if len(xs) * len(ys) != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but len(xs) = %d and len(ys) = %d",
			len(vals), len(xs), len(ys),
		))
	}

	return bi
}

// NewUniformBiLinear creates a bi-linear interpolator on top of a uniform
// grid with the values given by vals. The values of the x and y grid lines
// start at x0 and y0 and increase with steps of dx and dy, respectively.
// The vals grid is indexed in the usual way: vals(ix, iy) -> vals[ix + iy*nx].
//
// Panics if len(xs) * len(ys) != len(vals).
func NewUniformBiLinear(
	x0, dx float64, nx int,
	y0, dy float64, ny int,
	vals []float64,
) *BiLinear {

	bi := &BiLinear{}

	bi.xs.unifInit(x0, dx, nx)
	bi.ys.unifInit(y0, dy, ny)
	bi.nx = nx
	bi.vals = nil

	if nx * ny  != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but nx = %d and ny = %d",
			len(vals), nx, ny,
		))
	}

	return bi
}

// Eval evaluates the bi-linear interpolator at the coordinate (x, y).
//
// Panics if (x, y) is outside the range of the starting grid.
func (bi *BiLinear) Eval(x, y float64) float64 {
	ix1 := bi.xs.search(x)
	iy1 := bi.ys.search(y)
	ix2, iy2 := ix1 + 1, iy1 + 1

	x1, x2 := bi.xs.val(ix1), bi.xs.val(ix2)
	y1, y2 := bi.ys.val(iy1), bi.ys.val(iy2)
	i11, i12 := ix1 + bi.nx*iy1, ix1 + bi.nx*iy2
	i21, i22 := ix2 + bi.nx*iy1, ix2 + bi.nx*iy2

	v21, v22 := bi.vals[i11], bi.vals[i12]
	v11, v12 := bi.vals[i21], bi.vals[i22]

	dx, dy := x2 - x1, y2 - y1
	dx1, dx2 := x - x1, x2 - x
	dy1, dy2 := y - y1, y2 - y

	return (v11*dx2*dy2 + v12*dx1*dy2 +
		v21*dx2*dy1 + v22*dx1*dy1) / (dx*dy)
}

// EvalAll evaluates the interpolator at all the given (x, y) values. If an
// output array is given, the output is written to that array (the array is
// still returned as a convenience).
//
// If more than one output array is provided, only the first is used.
func (bi *BiLinear) EvalAll(xs, ys []float64, out ...[]float64) []float64 {
	if len(out) == 0 { out = [][]float64{ make([]float64, len(xs)) } }
	for i := range xs { out[0][i] = bi.Eval(xs[i], ys[i]) }
	return out[0]
}

//////////////////////////////
// TriLinear Implementation //
//////////////////////////////

type TriLinear struct {
	xs, ys, zs searcher
	vals []float64
	nx, ny int
}

func NewTriLinear(xs, ys, zs, vals []float64) *TriLinear {
	tri := &TriLinear{}
	tri.xs.init(xs)
	tri.ys.init(ys)
	tri.zs.init(zs)
	tri.nx = len(xs)
	tri.ny = len(ys)
	tri.vals = vals

	if len(xs) * len(ys) * len(zs) != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but len(xs) = %d, len(ys) = %d, and len(zs) = %d",
			len(vals), len(xs), len(ys), len(zs),
		))
	}

	return tri
}

func NewUniformTriLinear(
	x0, dx float64, nx int,
	y0, dy float64, ny int,
	z0, dz float64, nz int,
	vals []float64,
) *TriLinear {

	tri := &TriLinear{}

	tri.xs.unifInit(x0, dx, nx)
	tri.ys.unifInit(y0, dy, ny)
	tri.zs.unifInit(z0, dz, nz)
	tri.nx = nx
	tri.ny = ny
	tri.vals = nil

	if nx * ny * nz != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but nx = %d, ny = %d, and nz = %d",
			len(vals), nx, ny, nz,
		))
	}

	return tri
}

func (tri *TriLinear) Eval(x, y, z float64) float64 {
	ix1 := tri.xs.search(x)
	iy1 := tri.ys.search(y)
	iz1 := tri.ys.search(z)

	ix2, iy2, iz2 := ix1 + 1, iy1 + 1, iz1 + 1

	x1, x2 := tri.xs.val(ix1), tri.xs.val(ix2)
	y1, y2 := tri.ys.val(iy1), tri.ys.val(iy2)
	z1, z2 := tri.ys.val(iz1), tri.ys.val(iz2)

	v111 := tri.vals[ix1 + iy1*tri.nx + iz1*tri.nx*tri.ny]
	v112 := tri.vals[ix1 + iy1*tri.nx + iz2*tri.nx*tri.ny]
	v121 := tri.vals[ix1 + iy2*tri.nx + iz1*tri.nx*tri.ny]
	v122 := tri.vals[ix1 + iy2*tri.nx + iz2*tri.nx*tri.ny]

	v211 := tri.vals[ix2 + iy1*tri.nx + iz1*tri.nx*tri.ny]
	v212 := tri.vals[ix2 + iy1*tri.nx + iz2*tri.nx*tri.ny]
	v221 := tri.vals[ix2 + iy2*tri.nx + iz1*tri.nx*tri.ny]
	v222 := tri.vals[ix2 + iy2*tri.nx + iz2*tri.nx*tri.ny]

	xd := (x - x1) / (x2 - x1)
	yd := (y - y1) / (y2 - y1)
	zd := (z - z1) / (z2 - z1)

	c11 := v111*(1 - xd) + v211*xd
	c12 := v121*(1 - xd) + v221*xd
	c21 := v112*(1 - xd) + v212*xd
	c22 := v122*(1 - xd) + v222*xd

	c1 := c11*(1 - y2) + c21*yd
	c2 := c12*(1 - y2) + c22*yd

	return c1*(1 - zd) + c2*zd
}

func (tri *TriLinear) EvalAll(xs, ys, zs []float64, out ...[]float64) []float64 {
	if len(out) == 0 { out = [][]float64{ make([]float64, len(xs)) } }
	for i := range xs { out[0][i] = tri.Eval(xs[i], ys[i], zs[i]) }
	return out[0]	
}
