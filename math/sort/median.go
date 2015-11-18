package sort

func Percentile(xs []float64, p float64, buf ...[]float64) float64 {
	if len(xs) == 0 {
		panic("xs empty in call to Precentile(xs, ps)")
	} else if p > 1 || p < 0 {
		panic("percentile must be in the range [0, 1]")
	}

	n := int(p * float64(len(xs)))
	if n == len(xs) { n-- }

	return NthLargest(xs, n, buf...)
}

func Median(xs []float64, buf ...[]float64) float64 {
	if len(xs) == 0 {
		panic("xs empty in call to Median(xs)")
	}

	return NthLargest(xs, len(xs) / 2, buf...)
}

// Remember that this function is 1-indexed.
func NthLargest(xs []float64, n int, buf ...[]float64) float64 {
	if len(xs) == 0 {
		panic("xs empty in call to NthHighest(xs, ns)")
	}

	var medSlice []float64
	if len(buf) == 0 {
		medSlice = make([]float64, len(xs))
	} else {
		medSlice = buf[0]
		if len(medSlice) != len(xs) {
			panic("Length of buffer does not equal length of input array.")
		}
	}
	copy(medSlice, xs)

	return nthLargest(medSlice, n)
}

func nthLargest(xs []float64, n int) float64 {
	switch len(xs) {
	case 1:
		return xs[0]
	case 2:
		if xs[0] > xs[1] {
			xs[0], xs[1] = xs[1], xs[0]
		}

		if n <= 2 {
			return xs[2 - n]
		} else {
			panic("n in nthLargest(xs, n) too large")
		}
	case 3:
		high, mid, low := sort3(xs[0], xs[1], xs[2])
		switch n {
		case 1: return high
		case 2: return mid
		case 3: return low
		default:
			panic("n in nthLargest(xs, n) too large or too small")
		}
	default:
		pivIdx := partition(xs)
		nPiv := len(xs) - pivIdx
		if nPiv > n {
			return nthLargest(xs[pivIdx:], n)
		} else if nPiv < n {
			return nthLargest(xs[:pivIdx], n - nPiv)
		} else { // nPiv == n
			return xs[pivIdx]
		}
	}
}
