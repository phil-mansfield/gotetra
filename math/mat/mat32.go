package mat

import (
	"math"
)

// Matrix32 represents a matrix of float32 values.
type Matrix32 struct {
	Vals []float32
	Width, Height int
}

// LUFactors32 contains data fields neccessary for a number of matrix
// operations. Exporting this type allows calling routines to better manage
// their memory consumption and to prevent recomputing the same decomposition
// many times.
type LUFactors32 struct {
	lu Matrix32
	pivot []int
	d float32
}

// NewMatrix32 creates a matrix with the specified values and dimensions.
func NewMatrix32(vals []float32, width, height int) *Matrix32 {
	m := &Matrix32{}
	m.Init(vals, width, height)
	return m
}

// Init initializes a matrix with the specified values and dimensions.
func (m *Matrix32) Init(vals []float32, width, height int) {
	if width <= 0 {
		panic("width must be positive.")
	} else if height <= 0 {
		panic("height must be positive.")
	} else if width * height != len(vals) {
		panic("height * width must equal len(vals).")
	}

	m.Vals = vals
	m.Width, m.Height = width, height
}

// Mult multiplies two matrices together.
func (m1 *Matrix32) Mult(m2 *Matrix32) *Matrix32 {
	h, w := m1.Height, m2.Width
	out := NewMatrix32(make([]float32, h*w), w, h)
	return m1.MultAt(m2, out)
}

// Mult multiplies to matrices together and writes the result to the 
// specified matrix.
func (m1 *Matrix32) MultAt(m2, out *Matrix32) *Matrix32 {
	if m1.Width != m2.Height {
		panic("Multiplication of incompatible matrix sizes.")
	}

	for i := range out.Vals { out.Vals[i] = 0 }
	for i := 0; i < m1.Height; i++ {
		off := i*m1.Width
		for j := 0; j < m2.Width; j++ {
			outIdx := off + j
			for k := 0; k < m1.Width; k++ {
				m1Idx := off + k
				m2Idx := k*m2.Width + j
				out.Vals[outIdx] += m1.Vals[m1Idx] * m2.Vals[m2Idx]
			}
		}
	}

	return out
}

// Invert computes the inverse of a matrix.
func (m *Matrix32) Invert() *Matrix32 {
	lu := m.LU()
	inv := NewMatrix32(make([]float32, len(m.Vals)), m.Width, m.Height)
	return lu.InvertAt(inv)
}

// Determinant computes the determinant of a matrix.
func (m *Matrix32) Determinant() float32 {
	lu := m.LU()
	return lu.Determinant()
}

// SolveVector solves the equation m * xs = bs for xs.
func (m *Matrix32) SolveVector(bs []float32) []float32 {
	xs := make([]float32, len(bs))
	lu := m.LU()
	return lu.SolveVector(bs, xs)
}

// SolveMatrix solves the equation m * x = b for x.
func (m *Matrix32) SolveMatrix(b *Matrix32) *Matrix32 {
	x := NewMatrix32(make([]float32, len(m.Vals)), m.Width, m.Height)
	lu := m.LU()
	return lu.SolveMatrix(b, x)
}

// NewLUFactors32 creates an LUFactors32 instance of the requested dimensions.
func NewLUFactors32(n int) *LUFactors32 {
	luf := new(LUFactors32)

	luf.lu.Vals, luf.lu.Width, luf.lu.Height = make([]float32, n*n), n, n
	luf.pivot = make([]int, n)
	luf.d = 1

	return luf
}

// LU returns the LU decomposition of a matrix.
func (m *Matrix32) LU() *LUFactors32 {
	if m.Width != m.Height { panic("m is non-square.") }

	lu := NewLUFactors32(m.Width)
	m.LUFactorsAt(lu)
	return lu
}

// LUFactorsAt stores the LU decomposition of a matrix at the specified
// location.
func (m *Matrix32) LUFactorsAt(luf *LUFactors32) {
	if luf.lu.Width != m.Width || luf.lu.Height != m.Height {
		panic("luf has different dimenstions than m.")
	}
	copy(luf.lu.Vals, m.Vals)
	luf.factorizeInPlace()
}

func (lu *LUFactors32) factorizeInPlace() {
	m, n := &lu.lu, lu.lu.Width
	vv := make([]float32, n)
	lu.d = 1
	for i := 0; i < n; i++ {
		big := float32(0.0)
		for j := 0; j < n; j++ {
			tmp := float32(math.Abs(float64(m.Vals[i*n + j])))
			if tmp > big { big = tmp }
		}
		if big == 0 { panic("Singular Matrix.") }
		vv[i] = 1 / big
	}

	var imax int
	for k := 0; k < n; k++ {
		big := float32(0.0)
		for i := k; i < n; i++ {
			tmp := vv[i] * float32(math.Abs(float64(m.Vals[i*n + k])))
			if tmp > big {
				big = tmp
				imax = i
			}
		}
		if k != imax {
			for j := 0; j < n; j++ {
				m.Vals[imax*n + j], m.Vals[k*n + j] =
					m.Vals[k*n + j], m.Vals[imax*n + j]
			}
			lu.d = -lu.d
			vv[imax] = vv[k]
		}
		lu.pivot[k] = imax
		if m.Vals[k*n + k] == 0 { m.Vals[k*n + k] = 1e-20 }
		for i := k + 1; i < n; i++ {
			m.Vals[i*n + k] /= m.Vals[k*n + k]
			tmp := m.Vals[i*n + k]
			for j := k + 1; j < n; j++ {
				m.Vals[i*n + j] -= tmp*m.Vals[k*n + j]
			}
		}
	}
}

// SolveVector solves M * xs = bs for xs.
//
// bs and xs may poin to the same physical memory.
func (luf *LUFactors32) SolveVector(bs, xs []float32) []float32 {
	n := luf.lu.Width
	if n != len(bs) {
		panic("len(b) != luf.Width")
	} else if n != len(xs) {
		panic("len(x) != luf.Width")
	}

	// A x = b -> (L U) x = b -> L (U x) = b -> L y = b
	ys := xs
	if &bs[0] == &ys[0] {
		bs = make([]float32, n)
		copy(bs, ys)
	}

	// Solve L * y = b for y.
	forwardSubst32(n, luf.pivot, luf.lu.Vals, bs, ys)
	// Solve U * x = y for x.
	backSubst32(n, luf.lu.Vals, ys, xs)

	return xs
}

// Solves L * y = b for y.
// y_i = (b_i - sum_j=0^i-1 (alpha_ij y_j)) / alpha_ij
func forwardSubst32(n int, pivot []int, lu, bs, ys []float32) {
	for i := 0; i < n; i++ {
		ys[pivot[i]] = bs[i]
	}
	for i := 0; i < n; i++ {
		sum := float32(0.0)
		for j := 0; j < i; j++ {
			sum += lu[i*n + j] * ys[j]
		}
		ys[i] = (ys[i] - sum)
	}
}

// Solves U * x = y for x.
// x_i = (y_i - sum_j=i+^N-1 (beta_ij x_j)) / beta_ii
func backSubst32(n int, lu, ys, xs []float32) {
	for i := n - 1; i >= 0; i-- {
		sum := float32(0.0)
		for j := i + 1; j < n; j++ {
			sum += lu[i*n + j] * xs[j]
		}
		xs[i] = (ys[i] - sum) / lu[i*n + i]
	}
}

// SolveMatrix solves the equation m * x = b.
// 
// x and b may point to the same physical memory.
func (luf *LUFactors32) SolveMatrix(b, x *Matrix32) *Matrix32 {
	xs := x.Vals
	n := luf.lu.Width

	if b.Width != b.Height {
		panic("b matrix is non-square.")
	} else if x.Width != x.Height {
		panic("x matrix is non-square.") 
	} else if n != b.Width {
		panic("b matrix different size than m matrix.")
	} else if n != x.Width {
		panic("x matrix different size than m matrix.")
	}

	col := make([]float32, n)

	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			col[i] = xs[i*n + j]
		}
		luf.SolveVector(col, col)
		for i := 0; i < n; i++ {
			xs[i*n + j] = col[i]
		}
	}

	return x
}

// InvertAt inverts the matrix represented by the given LU decomposition
// and writes the results into the specified out matrix.
func (luf *LUFactors32) InvertAt(out *Matrix32) *Matrix32 {
	n := luf.lu.Width
	if out.Width != out.Height {
		panic("out matrix is non-square.")
	} else if n != out.Width {
		panic("out matrix different size than m matrix.")
	}

	for i := range out.Vals {
		out.Vals[i] = 0
	}
	for i := 0; i < n; i++ {
		out.Vals[i*n + i] = 1
	}

	luf.SolveMatrix(out, out)
	return out
}

// Determinant compute the determinant of of the matrix represented by the
// given LU decomposition.
func (luf *LUFactors32) Determinant() float32 {
	d := luf.d
	lu := luf.lu.Vals
	n := luf.lu.Width

	for i := 0; i < luf.lu.Width; i++ {
		d *= lu[i*n + i]
	}
	return d
}
