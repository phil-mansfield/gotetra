/*package interpolate provides routines for creating smooth analytic funcitons
through sparse or noisy data.
*/
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
	xs, ys, y2s []float64
	coeffs []splineCoeff

	incr bool
	// Usually the input data is uniform. This is our estimate of the point
	// spacing.
	dx float64
}

// NewSpline creates a spline based off a table of x and y values. The values
// must be sorted in increasing or decreasing order in x.
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

	sp.y2s = make([]float64, len(xs))
	sp.coeffs = make([]splineCoeff, len(xs)-1)
	sp.xs, sp.ys = xs, ys
	sp.Init(xs, ys)

	return sp
}

// Init reinitializes a spline to use a new sequence of points without doing
// any additional heap allocations. |xs| and |ys| must be the same as the
// previous point set.
func (sp *Spline) Init(xs, ys []float64) {
	if len(xs) != len(sp.xs) || len(ys) != len(sp.ys) {
		panic("Length of input arrays do not equal internal spline arrays.")
	}
	sp.xs, sp.ys = xs, ys

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
	sp.calcY2s()
	sp.calcCoeffs()
}

// Eval computes the value of the spline at the given point.
//
// x must be within the range of x values given to NewSpline().
func (sp *Spline) Eval(x float64) float64 {
	if x <= sp.xs[0] == sp.incr || x >= sp.xs[len(sp.xs)-1] == sp.incr {
		if x == sp.xs[0] { return sp.ys[0] }
		if x == sp.xs[len(sp.xs) - 1] { return sp.ys[len(sp.ys) - 1] }

		log.Fatalf("Point %g given to Spline.Eval() out of bounds [%g, %g].",
			x, sp.xs[0], sp.xs[len(sp.xs) - 1])
	}

	i := sp.bsearch(x)
	dx := x - sp.xs[i]
	a, b, c, d := sp.coeffs[i].a, sp.coeffs[i].b, sp.coeffs[i].c, sp.coeffs[i].d
	return a*dx*dx*dx + b*dx*dx + c*dx + d
}

// Deriv computes the derivative of spline at the given point to the
// specified order.
//
// x must be within the range of x values given to NewSpline().
func (sp *Spline) Deriv(x float64, order int) float64 {
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

// Integrate integrates the spline from lo to hi.
func (sp *Spline) Integrate(lo, hi float64) float64 {
	if lo > hi { return -sp.Integrate(hi, lo) }
	if lo < sp.xs[0] == sp.incr || lo > sp.xs[len(sp.xs)-1] == sp.incr {
		log.Fatalf("Low bound %g in Spline.Integrate() out of bounds.", lo)
	} else if hi < sp.xs[0] == sp.incr || hi > sp.xs[len(sp.xs)-1] == sp.incr {
		log.Fatalf("High bound %g in Spline.Integrate() out of bounds.", hi)
	}

	iLo, iHi := sp.bsearch(lo), sp.bsearch(hi)
	if iLo == iHi {
		return integTerm(&sp.coeffs[iLo], lo, hi)
	}
	sum := integTerm(&sp.coeffs[iLo], lo, sp.xs[iLo+1]) +
		integTerm(&sp.coeffs[iHi], sp.xs[iHi], hi)

	for i := iLo + 1; i < iHi; i++ {
		sum += integTerm(&sp.coeffs[i], sp.xs[i], sp.xs[i + 1])
	}
	return sum
}

func integTerm(coeff *splineCoeff, lo, hi float64) float64 {
	a, b, c, d := coeff.a, coeff.b, coeff.c, coeff.d
	dx := hi - lo
	return a*dx*dx*dx*dx/4 + b*dx*dx*dx/3 + c*dx*dx/2 + d*dx
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
		dx := xs[i+1] - xs[i]
		coeffs[i].a = (-y2s[i]/6 + y2s[i+1]/6) / dx
		coeffs[i].b = y2s[i] / 2
		coeffs[i].c = (ys[i+1] - ys[i])/dx + dx*(-y2s[i]/3 - y2s[i+1]/6)
		coeffs[i].d = ys[i]
	}
}

// TriTiagAt solves the system of equations
//
// | b0 c0 ..    |   | out0 |   | r0 |
// | a1 a2 c2 .. |   | out1 |   | r1 |
// | ..          | * | ..   | = | .. |
// | ..    an bn |   | outn |   | rn |
//
// For out0 .. outn in place in the given slice.
func TriDiagAt(as, bs, cs, rs, out []float64) {
	if len(as) != len(bs) || len(as) != len(cs) ||
		len(as) != len(out) || len(as) != len(rs) {

		log.Fatal("Length of arugments to TriDiagAt are unequal.")
	}

	tmp := make([]float64, len(as))

	beta := bs[0]
	if beta == 0 {
		log.Fatal("TriDiagAt cannot solve given system.")
	}
	out[0] = rs[0] / beta

	for i := 1; i < len(out); i++ {
		tmp[i] = cs[i-1] / beta
		beta = bs[i] - as[i]*tmp[i]
		if beta == 0 {
			log.Fatal("TriDiagAt cannot solve given system")
		}
		out[i] = (rs[i] - as[i]*out[i-1]) / beta

	}

	for i := len(out) - 2; i >= 0; i-- {
		out[i] -= tmp[i+1] * out[i+1]
	}
}

// TriTiag solves the system of equations
//
// | b0 c0 ..    |   | u0 |   | r0 |
// | a1 a2 c2 .. |   | u1 |   | r1 |
// | ..          | * | .. | = | .. |
// | ..    an bn |   | un |   | rn |
//
// For u0 .. un.
func TriDiag(as, bs, cs, rs []float64) []float64 {
	us := make([]float64, len(as))
	TriDiagAt(as, bs, cs, rs, us)
	return us
}
