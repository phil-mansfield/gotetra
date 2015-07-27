package geom

import (
	"testing"
)

func TestSolve(t *testing.T) {
	table := []struct{
		x1, y1, x2, y2, xintr, yintr float32
	} {
		{ 1, 1, 2, 2, 3, -3 },
		{ 1, 1, 2, 2, 1, 0 },
	}

	l1, l2 := new(Line), new(Line)

	for i, line := range table {
		l1.Init(line.x1, line.y1, line.xintr, line.yintr)
		l2.Init(line.x2, line.y2, line.xintr, line.yintr)
		
		x, y, ok := Solve(l1, l2)
		if !ok {
			t.Errorf("%d) Found that %g + %g * x intersects with %g + %g * " +
				"x are parallel.", i+1, l1.Y0, l1.M, l2.Y0, l2.M)
		}
		if !almostEq(x, line.xintr, 1e-5) || !almostEq(y, line.yintr, 1e-5) {
			t.Errorf("%d) Found that %g + %g * x intersects with %g + %g * " +
				"x at (%g, %g)\n", i+1, l1.Y0, l1.M, l2.Y0, l2.M, x, y)
		}
	}
}

func BenchmarkSolve(b *testing.B) {
	l1 := new(Line)
	l1.Init(5, 10, -1, 4)
	l2 := new(Line)
	l2.Init(-6, -5, 2, 2)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Solve(l1, l2)
	}
}
