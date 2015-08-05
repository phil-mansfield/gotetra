package interpolate

import (
	"log"
)

type splineCoeff struct {
	a, b, c, d float64
}

// Spline represents a 1D cubic spline which can be used to interpolate between
// points.
type Spline struct {
	xs, ys, y2s, sqrs []float64
	coeffs []splineCoeff

	incr bool

	// Usually the input data is uniform. This is our estimate of the point
	// spacing.
	dx float64
}

// NewSpline creates a spline based off a table of x and y values. The values
// must be sorted in increasing or decreasing order in x.
//
// xs and ys must not be modified throughout the lifetime of the Spline.
func NewSpline(xs, ys []float64) *Spline {
	if len(xs) != len(ys) {
		log.Fatalf(
			"Table given to NewSpline() has len(xs) = %d but len(ys) = %d.",
			len(xs), len(ys),
		)
	} else if len(xs) <= 1 {
		log.Fatalf("Table given to NewSpline() has length of %d.", len(xs))
	}

	sp := new(Spline)

	sp.xs, sp.ys = xs, ys
	sp.xs = make([]float64, len(xs))
	sp.ys = make([]float64, len(xs))
	sp.y2s = make([]float64, len(xs))
	sp.coeffs = make([]splineCoeff, len(xs)-1)

	if xs[0] < xs[1] {
		sp.incr = true
		for i := 0; i < len(xs)-1; i++ {
			if xs[i+1] < xs[i] {
				log.Fatal("Table given to NewSpline() not sorted.")
			}
		}
	} else {
		sp.incr = false
		for i := 0; i < len(xs)-1; i++ {
			if xs[i+1] > xs[i] {
				log.Fatal("Table given to NewSpline() not sorted.")
			}
		}
	}

	sp.dx = (xs[len(xs)-1] - xs[0]) / float64(len(xs)-1)

	copy(sp.xs, xs)
	copy(sp.ys, ys)
	sp.calcY2s()
	sp.calcCoeffs()
	return sp
}

// Eval computes the value of the spline at the given point.
//
// x must be within the range of x values given to NewSpline().
func (sp *Spline) Eval(x float64) float64 {
	if x < sp.xs[0] == sp.incr || x > sp.xs[len(sp.xs)-1] == sp.incr {
		log.Fatalf("Point %g given to Spline.Eval() out of bounds [%g, %g].",
			x, sp.xs[0], sp.xs[len(sp.xs) - 1])
	}

	i := sp.bsearch(x)
	dx := x - sp.xs[i]
	a, b, c, d := sp.coeffs[i].a, sp.coeffs[i].b, sp.coeffs[i].c, sp.coeffs[i].d
	return a*dx*dx*dx + b*dx*dx + c*dx + d
}

// Diff computes the derivative of spline at the given point to the
// specified order.
//
// x must be within the range of x values given to NewSpline().
func (sp *Spline) Diff(x float64, order int) float64 {
	if x < sp.xs[0] == sp.incr || x > sp.xs[len(sp.xs)-1] == sp.incr {
		log.Fatalf("Point %g given to Spline.Differentiate() out of bounds.", x)
	}

	i := sp.bsearch(x)
	dx := x - sp.xs[i]
	a, b, c, d := sp.coeffs[i].a, sp.coeffs[i].b, sp.coeffs[i].c, sp.coeffs[i].d
	switch order {
	case 0:
		return a*dx*dx*dx + b*dx*dx + c*dx + d
	case 1:
		return 3*a*dx*dx + 2*b*dx + c
	case 2:
		return 6*a*dx + 2*b
	case 3:
		return 6*a
	default:
		return 0
	}
}


// bsearch returns the the index of the largest element in xs which is smaller
// than x.
func (sp *Spline) bsearch(x float64) int {
	// Guess under the assumption of uniform spacing.
	guess := int((x - sp.xs[0]) / sp.dx)
	if guess >= 0 && guess < len(sp.xs)-1 &&
		(sp.xs[guess] <= x == sp.incr) &&
		(sp.xs[guess+1] >= x == sp.incr) {

		return guess
	}

	// Binary search.
	lo, hi := 0, len(sp.xs)-1
	for hi-lo > 1 {
		mid := (lo + hi) / 2
		if sp.incr == (x >= sp.xs[mid]) {
			lo = mid
		} else {
			hi = mid
		}
	}

	if lo == len(sp.xs) - 1 { 
		log.Fatalf("Point %g out of Spline bounds [%g, %g].",
			x, sp.xs[0], sp.xs[len(sp.xs) - 1])
	}
	return lo
}

// secondDerivative computes the second derivative at every point in the table
// given in NewSpline.
func (sp *Spline) calcY2s() {
	// These arrays do not escape to the heap.
	n := len(sp.xs)
	as, bs := make([]float64, n-2), make([]float64, n-2)
	cs, rs := make([]float64, n-2), make([]float64, n-2)

	// Solve for everything but the boundaries. The boundaries will be set to
	// zero. Better yet, they could be set to something computed via finite
	// differences.
	sp.y2s[0], sp.y2s[n-1] = 0, 0

	xs, ys := sp.xs, sp.ys
	for i := range rs {
		// j indexes into xs and ys.
		j := i + 1

		as[i] = (xs[j] - xs[j-1]) / 6
		bs[i] = (xs[j+1] - xs[j-1]) / 3
		cs[i] = (xs[j+1] - xs[j]) / 6
		rs[i] = ((ys[j+1] - ys[j]) / (xs[j+1] - xs[j])) -
			((ys[j] - ys[j-1]) / (xs[j] - xs[j-1]))
	}

	TriDiagAt(as, bs, cs, rs, sp.y2s[1: n-1])
}

func (sp *Spline) calcCoeffs() {
	coeffs, xs, ys, y2s := sp.coeffs, sp.xs, sp.ys, sp.y2s
	for i := range sp.coeffs {
		coeffs[i].a = (y2s[i+1] - y2s[i]) / (xs[i+1] - xs[i])
		coeffs[i].b = y2s[i] / 2
		coeffs[i].c = (ys[i+1] - ys[i]) / (xs[i+1] - xs[i]) -
			(xs[i+1] - xs[i])*(y2s[i]/3 + y2s[i+1]/5)
		coeffs[i].d = ys[i]
	}
}
