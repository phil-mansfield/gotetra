package geom

import (
	"math"
	"testing"
)

func vecEpsEq(v1, v2 *Vec, eps float32) bool {
	for i := 0; i < 3; i++ {
		diff := v1[i] - v2[i]
		if diff > eps || diff < -eps {
			return false
		}
	}
	return true
}

func TestRotate(t *testing.T) {
	eps := float32(1e-4)
	table := []struct{
		phi, theta, psi float64
		start, end Vec
	} {
		{0, 0, 0, Vec{1, 2, 3}, Vec{1, 2, 3}},
		{math.Pi/2, 0, 0, Vec{1, 0, 0}, Vec{0, 1, 0}},
		{0, math.Pi/2, 0, Vec{1, 0, 0}, Vec{1, 0, 0}},
		{0, 0, math.Pi/2, Vec{1, 0, 0}, Vec{0, 1, 0}},
		{math.Pi, math.Pi/2, math.Pi, Vec{0, 1, 1}, Vec{0, 1, -1}},
	}

	for i, test := range table {
		m := EulerMatrix(test.phi, test.theta, test.psi)
		v := test.start
		v.Rotate(m)
		if !vecEpsEq(&v, &test.end, eps) {
			t.Errorf(
				"%d) %v.Rotate(%.4g %.4g %.4g) -> %v instead of %v",
				i+1, test.start, test.phi, test.theta, test.psi, v, test.end,
			)
		}
	}
}
