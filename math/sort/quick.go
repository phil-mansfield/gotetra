package sort

import (
	"math/rand"
)

const (
	manualLen = 25
)

func sort3(x, y, z float64) (max, mid, min float64) {
	if x > y {
		if x > z {
			if y > z {
				return x, y, z
			} else {
				return x, z, y
			}
		} else {
			return z, x, y
		}
	} else {
		if y > z {
			if x > z {
				return y, x, z
			} else {
				return y, z, x
			}
		} else {
			return z, y, x
		}
	}
}

// Quick sorts an array in place via quicksort (and returns the result for
// convenience.)
func Quick(xs []float64) []float64 {
	if len(xs) < manualLen {
		return Shell(xs)
	} else {
		pivIdx := partition(xs)
		Quick(xs[0: pivIdx])
		Quick(xs[pivIdx: len(xs)])
		return xs
	}
}

func pivIdx(xs []float64) int {
	return rand.Intn(len(xs))
}

func partition(xs []float64) int {
	n, n2 := len(xs), len(xs) / 2
	// Take three values. The median will be the pivot, the other two will
	// be sentinel values so that we cna avoid bounds checks.
	max, mid, min := sort3(xs[0], xs[n2], xs[n-1])
	xs[0], xs[n2], xs[n-1] = min, mid, max
	xs[1], xs[n2] = xs[n2], xs[1]

	lo, hi := 1, n-1
	for {
		lo++
		for xs[lo] < mid { lo++ }
		hi--
		for xs[hi] > mid { hi-- }
		if hi < lo { break }
		xs[lo], xs[hi] = xs[hi], xs[lo]
	}
	
	// Swap the pivot into the middle
	xs[1], xs[hi] = xs[hi], xs[1]

	return hi
}
