package analyze

import (
	"math"
	"testing"
)

func almostEq(x, y, eps float64) bool {
	return x + eps > y && x - eps < y
}

func sliceAlmostEq(xs, ys []float64, eps float64) bool {
	if len(xs) != len(ys) { return false }
	for i := range xs {
		if !almostEq(xs[i], ys[i], eps) { return false }
	}
	return true
}

func TestNaNSplit(t *testing.T) {
	nan := math.NaN()
	eps := 1e-6

	table := []struct{
		in []float64
		out [][]float64
	} {
		{[]float64{}, [][]float64{}},
		{[]float64{nan}, [][]float64{}},
		{[]float64{1}, [][]float64{{1}}},
		{[]float64{1, nan}, [][]float64{{1}}},
		{[]float64{nan, 1}, [][]float64{{1}}},
		{[]float64{nan, 1, nan}, [][]float64{{1}}},
		{[]float64{1, nan, 2}, [][]float64{{1}, {2}}},

		{[]float64{1, 2, 3, nan, 5, nan, nan, 8, 9, nan},
			[][]float64{{1, 2, 3}, {5}, {8, 9}}},
	}

	for i, test := range table {
		out, _ := NaNSplit(test.in)
		if len(out) != len(test.out) {
			t.Errorf(
				"%d) Expected NaNSplit(%v) -> %v, but got %v",
				i + 1, test.in, test.out, out,
			)
		}
		for j := range out {
			if !sliceAlmostEq(out[j], test.out[j], eps) {
				t.Errorf(
					"%d) Expected NaNSplit(%v) -> %v, but got %v",
					i + 1, test.in, test.out, out,
				)
			}
		}
	}
}

func TestNaNSplitOptions(t *testing.T) {
	nan := math.NaN()
	eps := 1e-6

	in := []float64{1, 2, 3, nan, 5, nan, nan, 8, 9, nan}
	aux := []float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	auxOut := [][][]float64{{{0, 0, 0}, {0}, {0, 0}}}

	_, auxRes := NaNSplit(in, Aux(aux), AuxSets([][][]float64{{{}, {}, {}}}))
	if len(auxRes) != len(auxOut) {
		t.Errorf(
			"Expected NaNSplit(%v, Aux(%v)) -> %v, but got %v",
			in, aux, auxOut, auxRes,
		)
	}

	for j := range auxOut[0] {
		if !sliceAlmostEq(auxOut[0][j], auxRes[0][j], eps) {
			t.Errorf(
				"Expected NaNSplit(%v, Aux(%v)) -> %v, but got %v",
				in, aux, auxOut, auxRes,
			)
		}
	}
}
