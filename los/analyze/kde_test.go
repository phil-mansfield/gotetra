package analyze

import (
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
