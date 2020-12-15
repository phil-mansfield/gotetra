package interpolate

import (
	"fmt"
)

///////////////////////////
// Linear Implementation //
///////////////////////////

// Linear is a linear interpolator.
type Linear struct {
	xs   searcher
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

	return ((v2-v1)/(x2-x1))*(x-x1) + v1
}

// EvalAll evaluates the interpolator at all the given x values. If an output
// array is given, the output is written to that array (the array is still
// returned as a convenience).
//
// If more than one output array is provided, only the first is used.
func (lin *Linear) EvalAll(xs []float64, out ...[]float64) []float64 {
	if len(out) == 0 {
		out = [][]float64{make([]float64, len(xs))}
	}
	for i, x := range xs {
		out[0][i] = lin.Eval(x)
	}
	return out[0]
}

func (lin *Linear) Ref() Interpolator {
	panic("NYI")
}

type linearRef struct {
}

func (lin *linearRef) Eval(x float64) float64 {
	panic("NYI")
}

func (lin *linearRef) EvalAll(xs []float64, out ...[]float64) []float64 {
	panic("NYI")
}

func (lin *linearRef) Ref() Interpolator {
	panic("NYI")
}

/////////////////////////////
// BiLinear Implementation //
/////////////////////////////

// BiLinear is a bi-linear interpolator.
type BiLinear struct {
	xs, ys searcher
	vals   []float64
	nx     int
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

	if len(xs)*len(ys) != len(vals) {
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
	bi.vals = vals

	if nx*ny != len(vals) {
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
	ix2, iy2 := ix1+1, iy1+1
	if ix2 == bi.xs.n {
		ix1--
		ix2--
	}
	if iy2 == bi.ys.n {
		iy1--
		iy2--
	}

	x1, x2 := bi.xs.val(ix1), bi.xs.val(ix2)
	y1, y2 := bi.ys.val(iy1), bi.ys.val(iy2)

	i11, i12 := ix1+bi.nx*iy1, ix1+bi.nx*iy2
	i21, i22 := ix2+bi.nx*iy1, ix2+bi.nx*iy2

	v11, v12 := bi.vals[i11], bi.vals[i12]
	v21, v22 := bi.vals[i21], bi.vals[i22]

	dx, dy := x2-x1, y2-y1
	dx1, dx2 := x-x1, x2-x
	dy1, dy2 := y-y1, y2-y

	return (v11*dx2*dy2 + v12*dx2*dy1 +
		v21*dx1*dy2 + v22*dx1*dy1) / (dx * dy)
}

// EvalAll evaluates the interpolator at all the given (x, y) values. If an
// output array is given, the output is written to that array (the array is
// still returned as a convenience).
//
// If more than one output array is provided, only the first is used.
func (bi *BiLinear) EvalAll(xs, ys []float64, out ...[]float64) []float64 {
	if len(out) == 0 {
		out = [][]float64{make([]float64, len(xs))}
	}
	for i := range xs {
		out[0][i] = bi.Eval(xs[i], ys[i])
	}
	return out[0]
}

func (bi *BiLinear) Ref() BiInterpolator {
	panic("NYI")
}

type biLinearRef struct {
}

func (bi *biLinearRef) Eval(x, y float64) float64 {
	panic("NYI")
}

func (bi *biLinearRef) EvalAll(xs, ys []float64, out ...[]float64) []float64 {
	panic("NYI")
}

func (bi *biLinearRef) Ref() BiInterpolator {
	panic("NYI")
}

//////////////////////////////
// TriLinear Implementation //
//////////////////////////////

// TriLinear is a tri-linear interpolator.
type TriLinear struct {
	xs, ys, zs searcher
	vals       []float64
	nx, ny     int
}

// NewTriLinear creates a tri-linear interpolator on top of a grid with the
// values given by vals. The values of the x, y, and z grid lines are given by
// xs, ys, and zs respectively. The vals grid is indexed in the usual way:
// vals(ix, iy, iz) -> vals[ix + iy*nx + iz*nx*ny].
//
// Panics if len(xs) * len(ys) * len(zs) != len(vals).
func NewTriLinear(xs, ys, zs, vals []float64) *TriLinear {
	tri := &TriLinear{}
	tri.xs.init(xs)
	tri.ys.init(ys)
	tri.zs.init(zs)
	tri.nx = len(xs)
	tri.ny = len(ys)
	tri.vals = vals

	if len(xs)*len(ys)*len(zs) != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but len(xs) = %d, len(ys) = %d, and len(zs) = %d",
			len(vals), len(xs), len(ys), len(zs),
		))
	}

	return tri
}

// NewUniformTriLinear creates a tri-linear interpolator on top of a uniform
// grid with the values given by vals. The values of the x, y, and z grid lines
// start at x0, y0, and z0 and increase with steps of dx, dy, and dz,
// respectively. The vals grid is indexed in the usual way: vals(ix, iy, iz) ->
// vals[ix + iy*nx + iz*nx*ny].
//
// Panics if len(xs) * len(ys) != len(vals).
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
	tri.vals = vals

	if nx*ny*nz != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but nx = %d, ny = %d, and nz = %d",
			len(vals), nx, ny, nz,
		))
	}

	return tri
}

func (tri *TriLinear) Eval(x, y, z float64) float64 {
	ix := tri.xs.search(x)
	iy := tri.ys.search(y)
	iz := tri.ys.search(z)

	dix, diy, diz := 1, tri.nx, tri.nx*tri.ny
	i := ix + iy*tri.nx + iz*tri.ny*tri.nx

	v111 := tri.vals[i+0+0+0]
	v112 := tri.vals[i+0+0+diz]
	v121 := tri.vals[i+0+diy+0]
	v122 := tri.vals[i+0+diy+diz]
	v211 := tri.vals[i+dix+0+0]
	v212 := tri.vals[i+dix+0+diz]
	v221 := tri.vals[i+dix+diy+0]
	v222 := tri.vals[i+dix+diy+diz]

	x1, x2 := tri.xs.val(ix), tri.xs.val(ix+1)
	y1, y2 := tri.ys.val(iy), tri.ys.val(iy+1)
	z1, z2 := tri.ys.val(iz), tri.ys.val(iz+1)

	xd := (x - x1) / (x2 - x1)
	yd := (y - y1) / (y2 - y1)
	zd := (z - z1) / (z2 - z1)

	c11 := v111*(1-xd) + v211*xd
	c21 := v121*(1-xd) + v221*xd
	c12 := v112*(1-xd) + v212*xd
	c22 := v122*(1-xd) + v222*xd

	c1 := c11*(1-yd) + c21*yd
	c2 := c12*(1-yd) + c22*yd

	return c1*(1-zd) + c2*zd
}

func (tri *TriLinear) EvalAll(xs, ys, zs []float64, out ...[]float64) []float64 {
	if len(out) == 0 {
		out = [][]float64{make([]float64, len(xs))}
	}
	for i := range xs {
		out[0][i] = tri.Eval(xs[i], ys[i], zs[i])
	}
	return out[0]
}

func (tri *TriLinear) Ref() TriInterpolator {
	panic("NYI")
}

type triLinearRef struct {
}

func (tri *triLinearRef) Eval(x, y, z float64) float64 {
	panic("NYI")
}

func (tri *triLinearRef) EvalAll(xs, ys, zs []float64, out ...[]float64) []float64 {
	panic("NYI")
}

func (tri *triLinearRef) Ref() TriInterpolator {
	panic("NYI")
}
func NewLinearInterpolator(xs, vals []float64) Interpolator {
	return NewLinear(xs, vals)
}
func NewUniformLinearInterpolator(
	x0, dx float64, vals []float64,
) Interpolator {
	return NewUniformLinear(x0, dx, vals)
}
func NewBiLinearInterpolator(xs, ys, vals []float64) BiInterpolator {
	return NewBiLinear(xs, ys, vals)
}
func NewUniformBiLinearInterpolator(
	x0, dx float64, nx int,
	y0, dy float64, ny int, vals []float64,
) BiInterpolator {
	return NewUniformBiLinear(x0, dx, nx, y0, dy, ny, vals)
}
func NewTriLinearInterpolator(xs, ys, zs, vals []float64) TriInterpolator {
	return NewTriLinear(xs, ys, zs, vals)
}
func NewUniformTriLinearInterpolator(
	x0, dx float64, nx int,
	y0, dy float64, ny int,
	z0, dz float64, nz int, vals []float64,
) TriInterpolator {
	return NewUniformTriLinear(x0, dx, nx, y0, dy, ny, z0, dz, nz, vals)
}

//////////////////////////////
// TriLinearMulti Implementation //
//////////////////////////////

// TriLinearMulti is a tri-linear interpolator with multi-dimensional output
type TriLinearMulti struct {
	xs, ys, zs searcher
	nx, ny     int
	vals       [][]float64
	valsDim    int
}

// NewTriLinearMulti creates a tri-linear interpolator on top of a grid with the
// values given by vals. The values of the x, y, and z grid lines are given by
// xs, ys, and zs respectively. The vals grid is indexed in the usual way:
// vals(ix, iy, iz) -> vals[ix + iy*nx + iz*nx*ny].
//
// Panics if len(xs) * len(ys) * len(zs) != len(vals).
func NewTriLinearMulti(xs, ys, zs []float64, vals [][]float64) *TriLinearMulti {
	tri := &TriLinearMulti{}
	tri.xs.init(xs)
	tri.ys.init(ys)
	tri.zs.init(zs)
	tri.nx = len(xs)
	tri.ny = len(ys)
	tri.vals = vals
	tri.valsDim = len(vals[0])

	if len(xs)*len(ys)*len(zs) != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but len(xs) = %d, len(ys) = %d, and len(zs) = %d",
			len(vals), len(xs), len(ys), len(zs),
		))
	}

	return tri
}

// NewUniformTriLinearMulti creates a tri-linear interpolator on top of a uniform
// grid with the values given by vals. The values of the x, y, and z grid lines
// start at x0, y0, and z0 and increase with steps of dx, dy, and dz,
// respectively. The vals grid is indexed in the usual way: vals(ix, iy, iz) ->
// vals[ix + iy*nx + iz*nx*ny].
//
// Panics if len(xs) * len(ys) != len(vals).
func NewUniformTriLinearMulti(
	x0, dx float64, nx int,
	y0, dy float64, ny int,
	z0, dz float64, nz int,
	vals [][]float64,
) *TriLinearMulti {

	tri := &TriLinearMulti{}

	tri.xs.unifInit(x0, dx, nx)
	tri.ys.unifInit(y0, dy, ny)
	tri.zs.unifInit(z0, dz, nz)
	tri.nx = nx
	tri.ny = ny
	tri.vals = vals
	tri.valsDim = len(vals[0])

	if nx*ny*nz != len(vals) {
		panic(fmt.Sprintf(
			"len(vals) = %d, but nx = %d, ny = %d, and nz = %d",
			len(vals), nx, ny, nz,
		))
	}

	return tri
}

// Eval returns the interpolated values
func (tri *TriLinearMulti) Eval(x, y, z float64) []float64 {
	ix := tri.xs.search(x)
	iy := tri.ys.search(y)
	iz := tri.ys.search(z)
	x1, x2 := tri.xs.val(ix), tri.xs.val(ix+1)
	y1, y2 := tri.ys.val(iy), tri.ys.val(iy+1)
	z1, z2 := tri.ys.val(iz), tri.ys.val(iz+1)
	xd := (x - x1) / (x2 - x1)
	yd := (y - y1) / (y2 - y1)
	zd := (z - z1) / (z2 - z1)
	dix, diy, diz := 1, tri.nx, tri.nx*tri.ny
	i := ix + iy*tri.nx + iz*tri.ny*tri.nx

	out := make([]float64, tri.valsDim)
	for j := 0; j < tri.valsDim; j++ {
		v111 := tri.vals[i+0+0+0][j]
		v112 := tri.vals[i+0+0+diz][j]
		v121 := tri.vals[i+0+diy+0][j]
		v122 := tri.vals[i+0+diy+diz][j]
		v211 := tri.vals[i+dix+0+0][j]
		v212 := tri.vals[i+dix+0+diz][j]
		v221 := tri.vals[i+dix+diy+0][j]
		v222 := tri.vals[i+dix+diy+diz][j]

		c11 := v111*(1-xd) + v211*xd
		c21 := v121*(1-xd) + v221*xd
		c12 := v112*(1-xd) + v212*xd
		c22 := v122*(1-xd) + v222*xd

		c1 := c11*(1-yd) + c21*yd
		c2 := c12*(1-yd) + c22*yd
		out[j] = c1*(1-zd) + c2*zd
	}
	return out
}
