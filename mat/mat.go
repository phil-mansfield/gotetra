package mat

import (
	"math"
)

type Matrix struct {
	Vals []float64
	Width, Height int
}

type LUFactors struct {
	lu Matrix
	pivot []int
	d float64
}

func NewMatrix(vals []float64, width, height int) *Matrix {
	if width <= 0 {
		panic("width must be positive.")
	} else if height <= 0 {
		panic("height must be positive.")
	} else if width * height != len(vals) {
		panic("height * width must equal len(vals).")
	}

	return &Matrix{Vals: vals, Width: width, Height: height}
}

func NewLUFactors(n int) *LUFactors {
	luf := new(LUFactors)

	luf.lu.Vals, luf.lu.Width, luf.lu.Height = make([]float64, n*n), n, n
	luf.pivot = make([]int, n)
	luf.d = 1

	return luf
}

func (m *Matrix) LU() *LUFactors {
	if m.Width != m.Height { panic("m is non-square.") }

	lu := NewLUFactors(m.Width)
	m.LUFactorsAt(lu)
	return lu
}

func (m *Matrix) LUFactorsAt(luf *LUFactors) {
	if luf.lu.Width != m.Width || luf.lu.Height != m.Height {
		panic("luf has different dimenstions than m.")
	}

	n := m.Width
	scale := make([]float64, n)
	lu := luf.lu.Vals
	luf.d = 1
	copy(lu, m.Vals)

	for i := 0; i < n; i++ {
		iOffset := i * n

		max := 0.0
		for j := 0; j < n; j++ {
			tmp := math.Abs(lu[iOffset + j])
			if tmp > max { max = tmp }
		}
		if max == 0 { panic("m is singular") }
		scale[i] = 1 / max
	}

	for k := 0; k < n; k++ {
		max := 0.0
		maxi := 0
		for i := k; i < n; i++ {
			tmp := scale[i] * math.Abs(lu[i*n + k])
			if tmp > max {
				max = tmp
				maxi = i
			}
		}

		if k != maxi {
			kOffset, maxiOffset := n*k, n*maxi
			for j := 0; j < n; j++ {
				idx1, idx2 := kOffset + j, maxiOffset + j
				lu[idx1], lu[idx2] = lu[idx2], lu[idx1]
			}
			luf.d = -luf.d
			scale[maxi] = scale[k]
		}
		luf.pivot[k] = maxi

		if lu[n*k + k] == 0 { lu[n*k + k] = 1e-40 }

		kOffset := k*n
		for i := k + 1; i < n; i++ {
			iOffset := i*n
			lu[iOffset + k] /= lu[kOffset + k]
			tmp := lu[iOffset + k]
			for j := k + 1; j < n; j++ {
				lu[iOffset + j] -= tmp * lu[kOffset + j]
			}
		}
	}
}

// SolveVector solves M * xs = bs for xs.
//
// bs and xs may poin to the same physical memory.
func (luf *LUFactors) SolveVector(bs, xs []float64) {
	n := luf.lu.Width
	if n != len(bs) {
		panic("len(b) != luf.Width")
	} else if n != len(xs) {
		panic("len(x) != luf.Width")
	}

	copy(xs, bs)
	lu := luf.lu.Vals
	nzIdx := 0

	for i := 0; i < n ; i++ {
		piv := luf.pivot[i]
		sum := xs[piv]
		xs[piv] = xs[i]

		if nzIdx != 0 {
			iOffset := i*n
			for j := nzIdx - 1; j < i; j++ {
				sum -= lu[iOffset + j] * xs[j]
			}
		} else if sum != 0 {
			nzIdx = i+1
		}

		xs[i] = sum
	}

	for i := n - 1; i >= 0; i-- {
		sum := xs[i]
		iOffset := n * i
		for j := i + 1; j < n; j++ {
			sum -= lu[iOffset + j] * xs[j]
		}
		xs[i] = sum / lu[iOffset + i]
	}
}

// SolveMatrix solves the equation m * x = b.
// 
// x and b may point to the same physical memory.
func (luf *LUFactors) SolveMatrix(b, x *Matrix) {
	lu := luf.lu.Vals
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

	col := make([]float64, n)

	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			col[i] = lu[i*n + j]
		}
		luf.SolveVector(col, col)
		for i := 0; i < n; i++ {
			lu[i*n + j] = col[i]
		}
	}
}

func (luf *LUFactors) Inverse(out *Matrix) {
	n := lu.lu.Width
	if out.Width != out.Height {
		panic("out matrix is non-square.")
	} else if n != out.Width {
		panic("out matrix different size than m matrix.")
	}

	for i := range out.Vals {
		out.Vals[i] = 0
	}
	for i := 0; i < n; i++ {
		out.Vals[i*n + i]
	}

	luf.SolveMatrix(out, out)
}

func (luf *LUFactors) Determinant() float64 {
	d := luf.d
	lu := luf.lu.Vals
	for i := 0; i < luf.lu.Width; i++ {
		d *= lu[i*n + i]
	}
	return d
}
