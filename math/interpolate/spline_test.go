package interpolate

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"testing"
)

type func1D func(float64) float64

func splinePlots(xs, ys []float64) (spPlot, rawPlot string) {
	spXs := linspace(xs[0], xs[len(xs) - 1], 100)
	spYs := make([]float64, 100)
	sp := NewSpline(xs, ys)
	for i := 0; i < 100; i++ { spYs[i] = sp.Eval(spXs[i]) }

	spPlot = pyplotPlotString(spXs, spYs, "b", "Spline")
	rawPlot = pyplotPlotString(xs, ys, "ok", "Input")
	return spPlot, rawPlot
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

func TestPyplotSpline(t *testing.T) {
	f, err := ioutil.TempFile(".", "spline_pyplot")
	defer os.Remove(f.Name())
    if err != nil { t.Fatal(err.Error()) }

	fig0 := pyplotFigString(0)
	plot1, plot2 := splinePlots(
		[]float64{0, 1, 2, 3, 4},
		[]float64{2, 3, 4, 5, 6},
	)
	fig1 := pyplotFigString(0)
	plot3, plot4 := splinePlots(
		[]float64{0, 0.5, 1, 1.5, 2},
		[]float64{0, 0.25, 1, 2.25, 4},
	)
	fig2 := pyplotFigString(0)
	plot5, plot6 := splinePlots(
		[]float64{0, 0.5, 1, 1.5, 2},
		[]float64{0, 0.125, 1, 3.375, 8},
	)

	loc := "upper left"
	fileBody := fmt.Sprintf(`import numpy as np
import matplotlib.pyplot as plt
%s
%s
%s
plt.legend(fontsize=12, loc='%s')
%s
%s
%s
plt.legend(fontsize=12, loc='%s')
%s
%s
%s
plt.legend(fontsize=12, loc='%s')

plt.show()`, fig0, plot1, plot2, loc,
		fig1, plot3, plot4, loc,
		fig2, plot5, plot6, loc,
	)

	f.Write([]byte(fileBody))
    c := exec.Command("python", f.Name())
    var out bytes.Buffer
    c.Stdout = &out

    err = c.Run()
    if err != nil { t.Fatal(err.Error()) }
}


