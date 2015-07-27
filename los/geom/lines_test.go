package geom

import (
	"testing"
	"math/rand"
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
	lines := make([]Line, 1000)
	for i := range lines {
		lines[i].Init(
			rand.Float32(), rand.Float32(),
			rand.Float32(), rand.Float32(),
		)
	}

	j1, j2 := 0, 1

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Solve(&lines[j1], &lines[j2])
		if j2 + 1 == len(lines){ j2 -= len(lines) }
		if j1 + 1 == len(lines){ j1 -= len(lines) }
	}
}
