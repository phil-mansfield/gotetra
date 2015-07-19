package interpolate

import (
	"fmt"
	"testing"
)

func TestSpline(t *testing.T) {
	xs := []float64{ 0, 1, 1.5, 2, 3, 4, 5 }
	ys := []float64{ 2, 1, 1, 0, 2, 3, 1 }

	sp := NewSpline(xs, ys)
	dx := 5.0 / 100.0
	for x := 0.0; x <= 5.0; x += dx {
		fmt.Println(x, sp.Interpolate(x))
	}
}
