package interpolate

import (
	"math"
	"github.com/phil-mansfield/gotetra/math/mat"
)

type Kernel struct {
	cs []float64
	center int
}

type BoundaryCondition int
const (
	Periodic BoundaryCondition = iota
	Reflection
	ZeroPad
	Extension
)

// Get returns the value in xs that corresponds to the given index for a
// particular choice of bounday conditions.
func (b BoundaryCondition) Get(xs []float64, i int) float64 {
	switch {
	case i < 0:
		switch b {
		case Periodic: return xs[len(xs) - 1 + i]
		case Reflection: return xs[-(i + 1)]
		case ZeroPad: return 0
		case Extension: return xs[0]
		}
		panic("Impossible")
	case i >= len(xs):
		switch b {
		case Periodic: return xs[i - len(xs)]
		case Reflection: return xs [(len(xs)-1) - (i - len(xs)-1)]
		case ZeroPad: return 0
		case Extension: return xs[len(xs) - 1]
		}
		panic("Impossible")
	default: return xs[i]
	}
}

func (b BoundaryCondition) posGet(xs []float64, i int) float64 {
	switch b {
	case Periodic: return xs[i - len(xs)]
	case Reflection: return xs[len(xs) - 1 - i]
	case ZeroPad: return 0
	case Extension: return xs[0]
	}
	panic("Impossible")
}

func (b BoundaryCondition) negGet(xs []float64, i int) float64 {
	switch b {
	case Periodic: return xs[i + len(xs)]
	case Reflection: return xs [-1 - i]
	case ZeroPad: return 0
	case Extension: return xs[len(xs) - 1]
	}
	panic("Impossible")
}

// Convolve convolves a 1d data set according to the filter f. Boundary
// conditions are specified with b.
func (k *Kernel) Convolve(xs []float64, b BoundaryCondition) []float64 {
	out := make([]float64, len(xs))
	k.ConvolveAt(xs, b, out)
	return out
}

// ConvolveAt convolves a 1d data set according to the filter f. Boundary
// conditions are specified with b and the output is written to out.
func (k *Kernel) ConvolveAt(xs []float64, b BoundaryCondition, out []float64) {
	n := len(xs)
	nl, nr := k.center, len(k.cs) - 1 - k.center
	var x float64

	for i := 0; i <= nl; i++ {
		sum := 0.0
		for j, c := range k.cs {
			idx := i + j - k.center
			if idx < 0 {
				switch b {
				case Periodic: x = xs[(n - 1) + idx]
				case Reflection: x = xs[-(idx - 1)]
				case ZeroPad: x = 0
				case Extension: x = xs[0]
				}
			} else {
				x = xs[idx]
			}
			sum += x * c
		}
		out[i] = sum
	}
	for i := nl + 1; i < n - nr; i++ {
		sum := 0.0
		for j, c := range k.cs {
			idx := i + j - k.center
			sum += xs[idx] * c
		}
		out[i] = sum
	}
	for i := n - nr; i < n; i++ {
		sum := 0.0
		for j, c := range k.cs {
			idx := i + j - k.center
			if idx >= n {
				switch b {
				case Periodic: x = xs[idx - n]
				case Reflection: x = xs[2*(n-1) - idx]
				case ZeroPad: x = 0
				case Extension: x = xs[n - 1]
				}
			} else {
				x = xs[idx]
			}
			sum += x * c
		}
		out[i] = sum
	}
}

func (k *Kernel) normalize() {
	sum := 0.0
	for _, c := range k.cs { sum += c }
	for i := range k.cs { k.cs[i] /= sum }
}

func NewGaussianKernel(width int, sigma, dx float64) *Kernel {
	if width % 2 != 1 { panic("Kernel width must be odd.") }

	k := new(Kernel)
	k.cs = make([]float64, width)
	k.center = width / 2

	for i := 0; i <= k.center; i++ {
		x := float64(i - k.center) * dx
		k.cs[i] = math.Exp(-x*x / (2*sigma*sigma))
	}
	// Gaussians are symmetric: no need to compute again.
	for i := k.center + 1; i < len(k.cs); i++ {
		k.cs[i] = k.cs[len(k.cs) - 1 -  i]
	}

	k.normalize()
	return k
}

func NewTophatKernel(width int) *Kernel {
	if width % 2 != 1 { panic("Kernel width must be odd.") }
	
	k := new(Kernel)
	k.cs = make([]float64, width)
	k.center = width / 2

	for i := range k.cs { k.cs[i] = 1 }

	k.normalize()
	return k
}

func NewSavGolKernel(order, width int) *Kernel {
	if width % 2 != 1 {
		panic("Kernel width must be odd.")
	} else if width <= order {
		panic("Kernel width cannot be smaller than pOrder.")
	}

	k := new(Kernel)
	k.cs = make([]float64, width)
	k.center = width / 2

	k.savgol(order, 0)
	return k
}

func NewSavGolDerivKernel(dx float64, dOrder, pOrder, width int) *Kernel {
	if width % 2 != 1 {
		panic("Kernel width must be odd.")
	} else if dOrder > pOrder {
		panic("dOrder cannot be larger than pOrder.")
	} else if width <= pOrder {
		panic("Kernel width cannot be smaller than pOrder.")
	}

	k := new(Kernel)
	k.cs = make([]float64, width)
	k.center = width / 2

	k.savgol(pOrder, dOrder)
	fact := float64(factorial(dOrder))
	for i := range k.cs {
		k.cs[i] /= fact*dx
	}
	return k
}

func (k *Kernel) savgol(m, ld int) {
	n := len(k.cs) / 2

	aBuf := make([]float64, (m + 1)*(m + 1))
	a := mat.NewMatrix(aBuf, m + 1, m + 1)

	// "ipj" -> "i + j".
	for ipj := 0; ipj <= m * 2; ipj++ {
		ipj64 := float64(ipj)

		sum := 0.0
		if ipj == 0 { sum = 1.0 }

		for k := 1; k <= n; k++ { sum += math.Pow(float64(k), ipj64) }
		for k := 1; k <= n; k++ { sum += math.Pow(float64(-k), ipj64) }
		mm := 2*m - ipj
		if mm > ipj { mm = ipj }
		for imj := -mm; imj <= mm; imj += 2 {
			// (i+j) + (i-j) -> 2i, (i+j) - (i-j) -> 2j
			i, j := (ipj+imj) / 2, (ipj-imj) / 2
			a.Vals[j*a.Width + i] = sum
		}
	}
	
	lu := a.LU()
	b := make([]float64, m + 1)
	b[ld] = 1
	lu.SolveVector(b, b)

	for i := -n; i <= n; i++ {
		sum, fac, i64 := b[0], 1.0, float64(i)
		for mm := 1; mm < m+1; mm++ {
			fac *= i64
			sum += b[mm]*fac
		}
		k.cs[i + n] = sum
	}
}

func factorial(x int) int {
	prod := 1
	for i := 2; i <= x; i++ { prod *= i }
	return prod
}
