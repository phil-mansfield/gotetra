package sphere_halo

import (
	"math"
	"testing"

	"github.com/phil-mansfield/gotetra/los/geom"
)

func almostEq(x, y float64) bool {
	return x - y < 0.0001 && y - x < 0.0001
}

func TestHalfAngularWidth(t *testing.T) {
	tests := []struct { dist2, r2, res float64 }{
		{2, 1, math.Pi/4},
		{8, 4, math.Pi/4},
		{4, 4, math.Pi/2},
		{10 * 1000, 1, 0.01},
	}

	for i, test := range tests {
		res := halfAngularWidth(test.dist2, test.r2)
		if !almostEq(res, test.res) {
			t.Errorf("%d) expected halfAngularWidth(%g, %g) = %g, but got %g",
				i, test.dist2, test.r2, test.res, res)
		}
	}
}

func TestTwoValIntrDist(t *testing.T) {
	tests := []struct { dist, rad, b, res1, res2 float64 } {
		{2, 1, 0, 1, 3},
		{1, 1, 0, 0, 2},
		{5, 3, 3, 4, 4},
		{math.Sqrt(90), 5, 3, 5, 13},
	}

	for i, test := range tests {
		res1, res2 := twoValIntrDist(
			test.dist*test.dist, test.rad*test.rad, test.b)
		if !almostEq(res1, test.res1) || !almostEq(res2, test.res2) {
			t.Errorf("%d) expected twoValIntrDist(%g, %g, %g) = (%g, %g), " + 
				"but got (%g, %g)", i, test.dist*test.dist, test.rad*test.rad,
				test.b, test.res1, test.res2, res1, res2)
		}
	}
}

func TestOneValIntrDist(t *testing.T) {
	tests := []struct { dist, rad, b, dir, res float64 } {
		{0, 1, 0, +1, 1},
		{0, 1, 0, -1, 1},
		{1, 1, 0, +1, 2},
		{1, 1, 0, -1, 0},
		{5, 5, 3, +1, 8},
		{5, 5, 3, -1, 0},
		{math.Sqrt(2), math.Sqrt(5), 1, +1, 3},
		{math.Sqrt(2), math.Sqrt(5), 1, -1, 1},
	}
	for i, test := range tests {
		res := oneValIntrDist(
			test.dist*test.dist, test.rad*test.rad, test.b, test.dir,
		)
		if !almostEq(res, test.res) {
			t.Errorf("%d) expected oneValIntrDist(%g, %g, %g) = %g, " + 
				"but got %g", i, test.dist*test.dist, test.rad*test.rad,
				test.b, test.res, res)
		}
	}
}

func TestIdxRange(t *testing.T) {
	norms := []geom.Vec{{0, 0, 1}}
	origin := [3]float64{0, 0, 0}
	rMin, rMax := 1.0, 2.0
	bins, n := 10, 12

	h := SphereHalo{}
	h.Init(norms, origin, rMin, rMax, bins, n)

	tests := []struct {
		phiLo, phiHi float64
		iLo1, iHi1, iLo2, iHi2 int
	} {
		{math.Pi - 0.002, math.Pi - 0.001, 5, 6, 0, 0},
		{math.Pi - 0.001, math.Pi + 0.001, 5, 7, 0, 0},
		{2*math.Pi - 0.001, 2*math.Pi + 0.001, 11, 12, 0, 1},
		{-0.001, +0.001, 11, 12, 0, 1},
	}

	for i, test := range tests {
		iLo1, iHi1, iLo2, iHi2 := h.idxRange(test.phiLo, test.phiHi)
		if iLo1 != test.iLo1 || iHi1 != test.iHi1 ||
			iLo2 != test.iLo2 || iHi2 != test.iHi2 {
			
			t.Errorf("%d) expected idxRange(%g, %g) = (%d, %d, %d, %d), but " +
				"got (%d, %d, %d, %d)", i, test.phiLo, test.phiHi, test.iLo1,
				test.iHi1, test.iLo2, test.iHi2, iLo1, iHi1, iLo2, iHi2)
		}
	}
}
