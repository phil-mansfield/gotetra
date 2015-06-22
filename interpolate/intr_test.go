package interpolate

import (
	"testing"
)

func TestSpline(t *testing.T) {
	xs := []float64{5, 4, 3, 2, 1, 0}
	xs = []float64{0, 1, 2, 3, 4, 5}
	ys := make([]float64, 6)

	sp := NewSpline(xs, ys)
	sp.Interpolate(1.5)
}
