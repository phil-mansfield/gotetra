package sphere_halo

import (
	"math"
	"testing"
	
	"github.com/phil-mansfield/gotetra/los/geom"
)

func Float64SliceEq(xs, ys []float64) bool {
	if len(xs) != len(ys) { return false }
	for i := range xs {
		if !almostEq(xs[i], ys[i]) { return false }
	}
	return true
}

func TestInsertToRing(t *testing.T) {
	edges := make([]float64, 9)
	for i := range edges { edges[i] = math.Pow(10, 1 - float64(i)/4) }

	tests := []struct {
		los, n int
		vec geom.Vec
		radius float64
		res []float64
	} {
		{0, 8, geom.Vec{0, 0, 0}, 1, []float64{1, 1, 1, 1, 0, 0, 0, 0}},
		{0, 8, geom.Vec{0.25, 0, 0}, 0.75, []float64{1, 1, 1, 1, 0, 0, 0, 0}},
		{4, 8, geom.Vec{0.25, 0, 0}, 0.75, []float64{1,1,0.79588,0,0,0,0,0}},
		{0, 8, geom.Vec{float32(edges[4] + edges[3])/2, 0, 0},
			(edges[4] - edges[3])/2, []float64{0, 0, 0, 0, 1, 0, 0, 0}},
		{0, 8, geom.Vec{float32(edges[4] + edges[3])/2, 0, 0},
			(edges[4] - edges[3])/2, []float64{0, 0, 0, 0, 1, 0, 0, 0}},
		{0, 8, geom.Vec{float32(edges[8] + edges[0])/2, 0, 0},
			(edges[8] - edges[0])/2, []float64{1, 1, 1, 1, 1, 1, 1, 1}},
	}

	buf := make([]float64, 8)
	for i, test := range tests {
		h := SphereHalo{}
		h.Init([]geom.Vec{{0, 0, 1}}, [3]float64{1, 1, 1}, 0.1, 10, 8, test.n)
		h.insertToRing(test.vec, test.radius, 1, 0)
		h.GetRhos(0, test.los, buf)
		if !Float64SliceEq(buf, test.res) {
			t.Errorf("%d) h.InsertToRing() -> %.4g, but expected %.4g",
				i, buf, test.res)
		}
	}
}
