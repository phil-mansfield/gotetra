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
			&geom.Grid{[3]int{0, 0, 0}, 100, 100 * 100, 100 * 100 * 100},
			&geom.Grid{[3]int{0, 0, 0}, 100, 100 * 100, 100 * 100 * 100}},
		{10, 10, 0, 0, 0,
			&geom.Grid{[3]int{0, 0, 0}, 10, 10 * 10, 10 * 10 * 10},
			&geom.Grid{[3]int{0, 0, 0}, 100, 100 * 100, 100 * 100 * 100}},
		{20, 10, 1, 2, 3,
			&geom.Grid{[3]int{20, 40, 60}, 20, 20 * 20, 20 * 20 * 20},
			&geom.Grid{[3]int{0, 0, 0}, 200, 200 * 200, 200 * 200 * 200}},
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

func gridEq(g1, g2 *geom.Grid) bool {
	for i := 0; i < 3; i++ {
		if g1.Origin[i] != g2.Origin[i] {
			return false
		}
	}

	return g1.Width == g2.Width
}
