package interpolate

import (
	"log"
)

// Spline represents a 1D cubic spline which can be used to interpolate between
// points
type Spline struct {
	xs, ys, y2s, sqrs []float64

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
	sp.y2s = make([]float64, len(xs))
	sp.sqrs = make([]float64, len(xs)-1)

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

	sp.secondDerivative()
	for i := range sp.sqrs {
		sp.sqrs[i] = (xs[i+1] - xs[i]) * (xs[i+1] - xs[i])
	}

	return sp
}

// Interpolate interpolates the table of x and y values given in NewSpline to
// the point x.
//
// x must be within the range of x values given to NewSpline().
func (sp *Spline) Interpolate(x float64) float64 {
	if x < sp.xs[0] == sp.incr || x > sp.xs[len(sp.xs)-1] == sp.incr {
		log.Fatalf("Point %g given to Spline.Interpolate() out of bounds.", x)
	}

	lo := sp.bsearch(x)
	if lo == -1 {
		lo = 0
	}
	hi := lo + 1
	if hi == len(sp.xs) {
		hi = len(sp.xs) - 1
	}

	A := (sp.xs[hi] - x) / (sp.xs[hi] - sp.xs[lo])
	B := 1 - A
	C := (A*A*A - A) * sp.sqrs[lo] / 6
	D := (B*B*B - B) * sp.sqrs[lo] / 6
	return A*sp.ys[lo] + B*sp.ys[hi] + C*sp.y2s[lo] + D*sp.y2s[hi]
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
	return lo
}

// secondDerivative computes the second derivative at every point in the table
// given in NewSpline.
func (sp *Spline) secondDerivative() {
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

	TriDiagAt(as, bs, cs, rs, sp.y2s)
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

// TriTiagAt solves the system of equations
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
