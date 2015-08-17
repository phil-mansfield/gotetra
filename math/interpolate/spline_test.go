package interpolate

import (
	"math/rand"
	"sort"
	"testing"

	plt "github.com/phil-mansfield/pyplot"
)

type func1D func(float64) float64

func splinePlots(xs, ys []float64) {
	spXs := linspace(xs[0], xs[len(xs) - 1], 100)
	spYs := make([]float64, 100)
	sp := NewSpline(xs, ys)
	for i := 0; i < 100; i++ { spYs[i] = sp.Eval(spXs[i]) }

	plt.Plot(spXs, spYs, "b", plt.Label("Spline"), plt.Lw(3))
	plt.Plot(xs, ys, "ok", plt.Label("Input"), plt.Lw(3))
}

func (sp *Spline) terms() int { return len(sp.coeffs) }
func (sp *Spline) mapTerm(i int, xs []float64) []float64 {
	a, b, c, d := sp.coeffs[i].a, sp.coeffs[i].b, sp.coeffs[i].c, sp.coeffs[i].d
	ys := make([]float64, len(xs))
	for j, x := range xs {
		dx := x - sp.xs[i]
		ys[j] = d + c*dx + b*dx*dx + a*dx*dx*dx
	}
	return ys
}

func randSeq(n int, lo, hi float64) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = rand.Float64() * (hi - lo) + lo
	}
	return out
}

func TestPyplotSpline(t *testing.T) {
	plt.Reset()

	plt.Figure(plt.Num(0))
	plt.Title("Linear")
	splinePlots([]float64{0, 1, 2, 3, 4}, []float64{2, 3, 4, 5, 6})
	plt.Figure(plt.Num(1))
	plt.Title("Quardatic")
	splinePlots([]float64{0, 0.5, 1, 1.5, 2}, []float64{0, 0.25, 1, 2.25, 4})
	plt.Figure(plt.Num(2))
	rand.Seed(0)
	randXs := linspace(-1, 1, 10)
	sort.Float64Slice(randXs).Sort()
	randYs := randSeq(10, 0, 1)
	splinePlots(randXs, randYs)

	plt.Legend(plt.Loc("upper left"))
	plt.Show()
}


