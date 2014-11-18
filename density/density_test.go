package density

import (
	"testing"

	"github.com/phil-mansfield/gotetra/geom"
)

func BenchmarkNGP(b *testing.B) {
	g, bg := Bounds(100, 1, 0, 0, 0)
	rhos := make([]float64, 100*100*100)
	pts := make([]geom.Vec, 1000)
	intr := NewInterpolator(NearestGridPoint, g, bg, 1.0, rhos)

	b.ResetTimer()

	for i := 0; i < (b.N/len(pts))+1; i++ {
		intr.Interpolate(1, pts)
	}
}

func BenchmarkCIC(b *testing.B) {
	g, bg := Bounds(100, 1, 0, 0, 0)
	rhos := make([]float64, 100*100*100)
	pts := make([]geom.Vec, 1000)
	intr := NewInterpolator(CloudInCell, g, bg, 1.0, rhos)

	b.ResetTimer()

	for i := 0; i < (b.N/len(pts))+1; i++ {
		intr.Interpolate(1, pts)
	}
}

func TestBounds(t *testing.T) {
	table := []struct {
		cells, gridWidth int
		gx, gy, gz       int
		g, bg            *geom.Grid
	}{
		{100, 1, 0, 0, 0,
			geom.NewGrid(&[3]int{0, 0, 0}, 100),
			geom.NewGrid(&[3]int{0, 0, 0}, 100)},
		{10, 10, 0, 0, 0,
			geom.NewGrid(&[3]int{0, 0, 0}, 10),
			geom.NewGrid(&[3]int{0, 0, 0}, 100)},
		{20, 10, 1, 2, 3,
			geom.NewGrid(&[3]int{20, 40, 60}, 20),
			geom.NewGrid(&[3]int{0, 0, 0}, 200)},
	}

	for i, test := range table {
		g, bg := Bounds(test.cells, test.gridWidth, test.gx, test.gy, test.gz)

		if !gridEq(g, test.g) {
			t.Errorf("%d) Expected g to be %v, got %v \n", i, g, test.g)
		}
		if !gridEq(bg, test.bg) {
			t.Errorf("%d) Expected bg to be %v, got %v \n", i, bg, test.bg)
		}
	}
}

func TestNGPInterpolate(t *testing.T) {
	table := []struct {
		gx10, gy10, gz10 int
		mass, width      float64
		pts              []geom.Vec
		rhos             []float64
	}{
		// Tests for both NGP and CIC
		{0, 0, 0, 1.0, 2.0, []geom.Vec{},
			[]float64{0, 0, 0, 0, 0, 0, 0, 0}},
		{0, 0, 0, 1.0, 2.0, []geom.Vec{{0.5, 0.5, 0.5}},
			[]float64{1, 0, 0, 0, 0, 0, 0, 0}},
		{0, 0, 0, 1.0, 4.0, []geom.Vec{{0.5, 0.5, 0.5}},
			[]float64{0.125, 0, 0, 0, 0, 0, 0, 0}},
		{0, 0, 0, 1.0, 2.0, []geom.Vec{{1.5, 0.5, 0.5}},
			[]float64{0, 1, 0, 0, 0, 0, 0, 0}},

		{0, 0, 0, 1.0, 2.0, []geom.Vec{{1.5, 0.5, 0.5}},
			[]float64{0, 1, 0, 0, 0, 0, 0, 0}},
		{0, 0, 0, 1.0, 2.0, []geom.Vec{{0.5, 1.5, 0.5}},
			[]float64{0, 0, 1, 0, 0, 0, 0, 0}},
		{0, 0, 0, 1.0, 2.0, []geom.Vec{{0.5, 0.5, 1.5}},
			[]float64{0, 0, 0, 0, 1, 0, 0, 0}},

		{0, 0, 0, 1.0, 2.0, []geom.Vec{{9.5, 0.5, 0.5}},
			[]float64{0, 0, 0, 0, 0, 0, 0, 0}},
		{0, 0, 0, 1.0, 2.0, []geom.Vec{{0.5, 9.5, 0.5}},
			[]float64{0, 0, 0, 0, 0, 0, 0, 0}},
		{0, 0, 0, 1.0, 2.0, []geom.Vec{{0.5, 0.5, 9.5}},
			[]float64{0, 0, 0, 0, 0, 0, 0, 0}},

		{2, 4, 6, 1.0, 2.0, []geom.Vec{{2.5, 4.5, 6.5}},
			[]float64{1, 0, 0, 0, 0, 0, 0, 0}},
		{0, 0, 0, 1.0, 2.0, []geom.Vec{{0.5, 0.5, 0.5}, {1.5, 1.5, 1.5}},
			[]float64{1, 0, 0, 0, 0, 0, 0, 1}},
		// Tests for only NGP
		{0, 0, 0, 1.0, 2.0, []geom.Vec{{0.75, 0.75, 0.75}},
			[]float64{1, 0, 0, 0, 0, 0, 0, 0}},
	}

	for i, test := range table {
		g, bg := Bounds(2, 5, test.gx10/2, test.gy10/2, test.gz10/2)

		rhos := make([]float64, len(test.rhos))
		intr := NewInterpolator(NearestGridPoint, g, bg, test.width, rhos)
		intr.Interpolate(test.mass, test.pts)

		if !sliceEq(rhos, test.rhos) {
			t.Errorf("%d) Expected density grid %v, got %v.\n",
				i, test.rhos, rhos)
		}
	}
}

func gridEq(g1, g2 *geom.Grid) bool {
	for i := 0; i < 3; i++ {
		if g1.Origin[i] != g2.Origin[i] {
			return false
		}
	}

	return g1.Width == g2.Width
}

func sliceEq(xs1, xs2 []float64) bool {
	if len(xs1) != len(xs2) {
		return false
	}

	for i := range xs1 {
		if xs1[i] != xs2[i] {
			return false
		}
	}

	return true
}
