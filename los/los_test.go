package los

import (
	"testing"
)

func intSliceEq(xs, ys []int) bool {
	if len(xs) != len(ys) { return false }
	for i := range xs {
		if xs[i] != ys[i] { return false }
	}
	return true
}

func TestSplits(t *testing.T) {
	table := []struct {
		intr []bool
		workers int
		out []int
	} {
		{[]bool{true, true, true, true}, 1, []int{0, 4}},
		{[]bool{false, true, true, false}, 1, []int{0, 3}},
		{[]bool{true, true, true, true}, 2, []int{0, 2, 4}},
		{[]bool{true, true, true, false}, 2, []int{0, 1, 4}},
		{[]bool{false, true, true, false, true,
			false, true, true, true, false}, 3, []int{0, 3, 7, 10}},
	}

	for i, test := range table {
		res, _ := splits(test.intr, test.workers)
		if !intSliceEq(res, test.out) {
			t.Errorf("%d) Expected splits(%v, %v) = %v, but got %v",
				i+1, test.intr, test.workers, test.out, res)
		}
	}
}
