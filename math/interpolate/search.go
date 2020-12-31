package interpolate

import (
	"fmt"
)

type searcher struct {
	xs          []float64
	x0, dx, lim float64
	n           int
	unif, incr  bool
}

func (s *searcher) init(xs []float64) {
	s.xs = xs
	s.x0 = xs[0]
	s.lim = xs[len(xs)-1]
	s.dx = (s.lim - s.x0) / float64(len(xs)-1)
	s.n = len(xs)
	s.unif = false
	s.incr = s.dx > 0
}

func (s *searcher) unifInit(x0, dx float64, n int) {
	s.xs = nil
	s.x0 = x0
	s.lim = float64(n-1)*dx + x0
	s.dx = dx
	s.n = n
	s.unif = true
	s.incr = s.dx > 0
}

func (s *searcher) search(x float64) int {
	if x > s.lim || x < s.x0 {
		panic(fmt.Sprintf(
			"Value %g out of range bounds [%g, %g]", x, s.x0, s.lim,
		))
	}

	if s.unif {
		idx := int((x - s.x0) / s.dx)
		if idx == s.n-1 {
			idx--
		}
		return idx
	} else {

		// Guess under the assumption of uniform spacing.
		guess := int((x - s.xs[0]) / s.dx)
		if guess >= 0 && guess < len(s.xs)-1 &&
			(s.xs[guess] <= x == s.incr) &&
			(s.xs[guess+1] >= x == s.incr) {

			return guess
		}

		// Binary search.
		lo, hi := 0, s.n-1
		for hi-lo > 1 {
			mid := (lo + hi) / 2
			if s.incr == (x >= s.xs[mid]) {
				lo = mid
			} else {
				hi = mid
			}
		}

		return lo
	}
}

func (s *searcher) val(i int) float64 {
	if s.unif {
		return float64(i)*s.dx + s.x0
	}
	return s.xs[i]
}
