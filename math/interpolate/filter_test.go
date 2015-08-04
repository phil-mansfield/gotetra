package interpolate

import (
	"testing"
)

func BenchmarkConvolveArray200Filter21(b *testing.B) {
	out, xs := make([]float64, 200), make([]float64, 200)
	k := NewTophatKernel(21)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		k.ConvolveAt(xs, Extension, out)
	}
}

func BenchmarkNewSavGolKernel21(b *testing.B) {
	for i := 0; i < b.N; i++ {
		NewSavGolKernel(4, 21)
	}
}

func almostEq(xs, ys []float64) bool {
	if len(xs) != len(ys) { return false }
	eps := 1e-3
	for i := range xs {
		if !(xs[i] + eps > ys[i] && xs[i] - eps < ys[i]) {
			return false
		}
	}
	return true
}

func TestSavGolKernel(t *testing.T) {
	table := []struct{
		order, width int
		cs []float64
	} {
		{2, 5, []float64{-0.086, 0.343, 0.486, 0.343, -0.086}},
		{2, 11, []float64{-0.084, 0.021, 0.103, 0.161, 0.196,
			0.207, 0.196, 0.161, 0.103, 0.021, -0.084}},
		{4, 9, []float64{0.035, -0.128, 0.070, 0.315,
			0.417, 0.315, 0.070, -0.128, 0.035}},
		{4, 11, []float64{0.042, -0.105, -0.023, 0.140, 0.280,
			0.333, 0.280, 0.140, -0.023, -0.105, 0.042}},
	}
	for i, test := range table {
		k := NewSavGolKernel(test.order, test.width)
		if !almostEq(k.cs, test.cs) {
			t.Errorf("%d) Expected %.3f for coefficients. Got %.3f.",
				i+1, test.cs, k.cs)
		}
	}
}
