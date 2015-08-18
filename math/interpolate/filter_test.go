package interpolate

import (
	"math"
	"math/rand"
	"testing"

	plt "github.com/phil-mansfield/pyplot"
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

func linspace(low, high float64, n int) []float64 {
	xs := make([]float64, n)
	dx := (high - low) / float64(n - 1)
	for i := range xs { xs[i] = low + dx*float64(i) }
	xs[len(xs) - 1] = high
	return xs
}

func gaussian(x0, sigma, A, x float64) float64 {
	return A * math.Exp(-(x-x0)*(x-x0)/(2*sigma*sigma))
}

func rawSavGolFunc(x float64) float64 {
	return gaussian(2, 1, 1.5, x) + gaussian(4, 0.5, 1.5, x) + 
		gaussian(5.5, 0.125, 1.5, x) + gaussian(0.5, 0.125, 1.5, x)
}

func TestPyplotSavGol(t *testing.T) {

	xs := linspace(0, 6, 200)
	rawYs := make([]float64, 200)
	noiseYs := make([]float64, 200)

	rand.Seed(0)
	for i, x := range xs {
		rawYs[i] = rawSavGolFunc(x)
		noiseYs[i] = rawYs[i] + rand.Float64() - 0.5
	}

	window := 41
	windowSize := float64(window) / float64(len(xs)) * (xs[len(xs)-1]-xs[0])
	sigma := windowSize / 5

	tk := NewTophatKernel(window)
	gk:= NewGaussianKernel(window, sigma, xs[1]-xs[0])
	sgk := NewSavGolKernel(4, window)

	plt.Reset()

	plt.Plot(xs, rawYs, "m", plt.Label("Underlying Function"), plt.LW(3))
	plt.Plot(xs, noiseYs, "k", plt.Label("Noisy Function"), plt.LW(3))
	plt.Plot(xs, tk.Convolve(noiseYs, Extension), "r",
		plt.Label("Tophat"), plt.LW(3))
	plt.Plot(xs, gk.Convolve(noiseYs, Extension), "g",
		plt.Label("Gaussian"), plt.LW(3))
	plt.Plot(xs, sgk.Convolve(noiseYs, Extension), "b",
		plt.Label("Savitzky-Golay"), plt.LW(3))

	plt.Legend(plt.Loc("lower left"), plt.FrameOn(false))
	plt.Show()
}
