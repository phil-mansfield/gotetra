package analyze

import (
	"math"
	"math/rand"
	"testing"

	plt "github.com/phil-mansfield/pyplot"
)

func TestGaussianKDE(t *testing.T) {
	plt.Reset()

	xs := make([]float64, 20)
	for i := range xs { xs[i] = rand.Float64() }
	sp := GaussianKDE(xs, 0.1, 0, 1, 100)
	evalXs, evalYs := make([]float64, 200), make([]float64, 200)
	ptYs := make([]float64, len(xs))
	for i := 0; i < len(evalXs); i++ {
		evalXs[i] = float64(i) / float64(len(evalXs) - 1)
	}
	evalXs[len(evalXs) - 1] = 1
	for i, x := range evalXs { evalYs[i] = sp.Eval(x) }
	for i, x := range xs { ptYs[i] = sp.Eval(x) }

	plt.Plot(xs, ptYs, "ok")
	plt.Plot(evalXs, evalYs, "r", plt.LW(3))

	plt.Show()
}

func BenchmarkGaussianKDE(b *testing.B) {
	xs := make([]float64, 1028)
	for i := range xs { xs[i] = rand.Float64() }
	for i := 0; i < b.N; i++ { GaussianKDE(xs, 0.1, 0, 1, 100) }
}

func TestGaussianKDETree(t *testing.T) {
	n := 1000
	rs, ths := make([]float64, n), make([]float64, n)
	xs, ys  := make([]float64, n), make([]float64, n)

	for i := 0; i < n; i++ {
		rs[i] = rand.NormFloat64() / 2 + 4
		if rs[i] < 0 { rs[i] = rand.Float64() + 3.5 }
		ths[i] = rand.Float64() * 2 * math.Pi
	}
	for i := 0; i < n; i++ {
		sin, cos := math.Sincos(ths[i])
		xs[i], ys[i] = cos * rs[i], sin * rs[i]
	}
	

	kt := NewKDETree(rs, ths, 5)
	f := kt.GetRFunc(5, Cartesian)
	spXs, spYs := make([]float64, 200), make([]float64, 200)
	for i := 0; i < len(spXs) - 1; i++ {
		th := 2 * math.Pi * (float64(i) + 0.5) / 200
		sin, cos := math.Sincos(th)
		r := f(th)
		spXs[i], spYs[i] = r * cos, r * sin
	}
	spXs[len(spXs) - 1], spYs[len(spXs) - 1] = spXs[0], spYs[0]

	plt.Reset()
	plt.Plot(xs, ys, "ow")
	plt.Plot(spXs, spYs, "r", plt.LW(3))
	plt.Show()
}

func BenchmarkGaussianKDETree(b *testing.B) {
	n := 1000
	rs, ths := make([]float64, n), make([]float64, n)
	for i := 0; i < n; i++ {
		rs[i] = rand.NormFloat64() + 4
		if rs[i] < 0 { rs[i] = rand.Float64() + 3.5 }
		ths[i] = rand.Float64() * 2 * math.Pi
	}
	
	for i := 0; i < b.N; i++ { NewKDETree(rs, ths, 2) }
}
