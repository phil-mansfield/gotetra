package geom

import (
	"math"
	"math/rand"
	"testing"
)

func randomTranslations(n int) []Vec {
	vs := make([]Vec, n)
	for i := range vs {
		for j := 0; j < 3; j++ {
			vs[i][j] = rand.Float32() - 0.5
		}
	}
	return vs
}

func TestPluckerTranslate(t *testing.T) {
	tet := Tetra{Vec{4, 0, 1}, Vec{0, 4, 1}, Vec{0, 0, 1}, Vec{0, 0, 2}}
	P, L := Vec{1, 1, 0}, Vec{0, 0, 1}
	p := new(AnchoredPluckerVec)
	pt := new(PluckerTetra)
	p.Init(&P, &L)
	pt.Init(&tet)
	targetEnter := float32(1.0)
	
	n := 1000
	dxs := append([]Vec{{0, 0, 0}}, randomTranslations(n - 1)...)
	w := new(IntersectionWorkspace)
	
	for i := range dxs {
		p.Translate(&dxs[i])
		pt.Translate(&dxs[i])
		tet.Translate(&dxs[i])
		
		enter, _, ok := w.IntersectionDistance(pt, &tet, p)

		if !ok {
			t.Errorf(
				"%d) No intersection with dx = %v", i+1, dxs[i],
			)
		} else if !almostEq(enter, targetEnter, 1e-4) {
			t.Errorf(
				"%d) Intersection distance of %g instead of %g with dx = %v",
				i + 1, enter, targetEnter, dxs[i],
			)
		}
	}
}

func almostEq64(x, y float64) bool { return x - 1e-4 < y && x + 1e-4 > y }

func TestTetraVolume(t *testing.T) {
	table := []struct {
		t Tetra
		vol float64
	} {
		{Tetra{{0, 2, 0}, {3, 0, 0}, {0, 0, 1}, {0, 0, 0}}, 1},
	}

	for i, test := range table {
		vol := test.t.Volume()
		if !almostEq64(vol, test.vol) {
			t.Errorf("%d) Expected %v.Volume() = %g, got %g.",
				i+1, test.t, test.vol, vol)
		}
	}
}

func TestSphereSphereIntersect(t *testing.T) {
	table := []struct {
		s1, s2 Sphere
		res bool
	} {
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 0, 0}, 1}, true},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 0, 0}, 2}, true},
		{Sphere{Vec{0, 0, 0}, 2}, Sphere{Vec{0, 0, 0}, 1}, true},
		
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 0, 1}, 1}, true},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 1, 0}, 1}, true},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{1, 0, 0}, 1}, true},

		{Sphere{Vec{0, 0, 0}, 3}, Sphere{Vec{0, 0, 1.5}, 1}, true},
		{Sphere{Vec{0, 0, 0}, 3}, Sphere{Vec{0, 1.5, 0}, 1}, true},
		{Sphere{Vec{0, 0, 0}, 3}, Sphere{Vec{1.5, 0, 0}, 1}, true},

		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 0, 3}, 1}, false},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 3, 0}, 1}, false},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{3, 0, 0}, 1}, false},
	}

	for i, test := range table {
		if test.s1.SphereIntersect(&test.s2) != test.res {
			t.Errorf("%d) %v.SphereIntersect(%v) -> %v",
				i+1, test.s1, test.s2, test.res)
		}
	}
}

func TestSphereSphereContain(t *testing.T) {
	table := []struct {
		s1, s2 Sphere
		res bool
	} {
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 0, 0}, 1}, false},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 0, 0}, 2}, false},
		{Sphere{Vec{0, 0, 0}, 2}, Sphere{Vec{0, 0, 0}, 1}, true},
		
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 0, 1}, 1}, false},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 1, 0}, 1}, false},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{1, 0, 0}, 1}, false},

		{Sphere{Vec{0, 0, 0}, 3}, Sphere{Vec{0, 0, 1.5}, 1}, true},
		{Sphere{Vec{0, 0, 0}, 3}, Sphere{Vec{0, 1.5, 0}, 1}, true},
		{Sphere{Vec{0, 0, 0}, 3}, Sphere{Vec{1.5, 0, 0}, 1}, true},

		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 0, 3}, 1}, false},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{0, 3, 0}, 1}, false},
		{Sphere{Vec{0, 0, 0}, 1}, Sphere{Vec{3, 0, 0}, 1}, false},
	}

	for i, test := range table {
		if test.s1.SphereContain(&test.s2) != test.res {
			t.Errorf("%d) %v.SphereContain(%v) -> %v",
				i+1, test.s1, test.s2, test.res)
		}
	}
}

func TestSphereLineSegmentIntersect(t *testing.T) {
	eps := float32(1e-4)
	rt2 := 1/float32(math.Sqrt(2))
	table := []struct {
		s Sphere
		ls LineSegment
		enter, exit float32
		enters, exits bool
	} {
		{Sphere{Vec{0, 0, 0}, 1},
			LineSegment{Vec{0, 0, 0}, Vec{1, 0, 0}, -2, +2},
			-1, +1, true, true},
		{Sphere{Vec{0, 0, 0}, 1},
			LineSegment{Vec{0, 0, 0}, Vec{rt2, rt2, 0}, -1.1, +2},
			-1, +1, true, true},
		{Sphere{Vec{0, 0, 10}, 1},
			LineSegment{Vec{0, 0, 0}, Vec{1, 0, 0}, -2, +2},
			0, 0, false, false},
		{Sphere{Vec{7, 0, 0}, 1},
			LineSegment{Vec{0, 0, 0}, Vec{1, 0, 0}, 1, 10},
			6, 8, true, true},
		{Sphere{Vec{10, 0, 0}, 1},
			LineSegment{Vec{0, 0, 0}, Vec{1, 0, 0}, 1, 10},
			9, 0, true, false},
		{Sphere{Vec{1, 0, 0}, 1},
			LineSegment{Vec{0, 0, 0}, Vec{1, 0, 0}, 1, 10},
			0, 2, false, true},
	}

	for i, test := range table {
		enter, exit, enters, exits := test.s.LineSegmentIntersect(&test.ls)
		if enters != test.enters || exits != test.exits ||
			(enters && !almostEq(enter, test.enter,eps)) || 
			(exits && !almostEq(exit, test.exit,eps)) {
			t.Errorf(
				"%d) expected intersect result of (%g %g %v %v), got " + 
					"(%g %g %v %v) ", i + 1, test.enter, test.exit, test.enters,
				test.exits, enter, exit, enters, exits,
			)
		}
	}
}

func BenchmarkVecTranslate(b *testing.B) {
	n := 1000
	dxs := randomTranslations(n)
	v := new(Vec)
	for i := 0; i < b.N; i++ {
		for j := 0; j < 3; j++ { v[j] += dxs[i%n][j] }
	}
}

func BenchmarkTetraTranslate(b *testing.B) {
	n := 1000
	dxs := randomTranslations(n)
	t := new(Tetra)
	for i := 0; i < b.N; i++ { t.Translate(&dxs[i % n]) }	
}

func BenchmarkPluckerVecTranslate(b *testing.B) {
	n := 1000
	dxs := randomTranslations(n)
	p := new(PluckerVec)
	for i := 0; i < b.N; i++ { p.Translate(&dxs[i % n]) }
}

func BenchmarkPluckerTetraTranslate(b *testing.B) {
	n := 1000
	dxs := randomTranslations(n)
	pt := new(PluckerTetra)
	for i := 0; i < b.N; i++ { pt.Translate(&dxs[i % n])}
}

func BenchmarkSphereIntersect(b *testing.B) {
	n := 1000
	ts := make([]Tetra, n)
	for i := range ts {
		for j := 0; j < 4; j++ {
			for k := 0; k < 3; k++ {
				ts[i][j][k] = rand.Float32()
			}
		}
	}

	ss := make([]Sphere, n)
	for i := range ss {
		ts[i].BoundingSphere(&ss[i])
	}

	b.ResetTimer()
	s := ss[0]
	idx := 0
	for i := 0; i < b.N; i++ {
		s.SphereIntersect(&ss[idx])
		idx++
		if idx == n { idx = 0 }
	}
}

func BenchmarkTetraBoundingSphere(b *testing.B) {
	n := 1000
	ts := make([]Tetra, n)
	for i := range ts {
		for j := 0; j < 4; j++ {
			for k := 0; k < 3; k++ {
				ts[i][j][k] = rand.Float32()
			}
		}
	}

	ss := make([]Sphere, n)

	b.ResetTimer()
	idx := 0
	for i := 0; i < b.N; i++ {
		ts[idx].BoundingSphere(&ss[idx])

		idx++
		if idx == n { idx = 0 }
	} 
}

func BenchmarkTetraVolume(b *testing.B) {
	t := new(Tetra)
	for i := 0; i < 4; i++ {
		for j := 0; j < 3; j++ {
			t[i][j] = rand.Float32()
		}
	}

	for i := 0; i < b.N; i++ { t.Volume() }
}
