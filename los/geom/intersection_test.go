package geom

import (
	"testing"
	"math"
	"math/rand"
)

func randomizeTetra(t *Tetra, low, high float32) {
	for v := 0; v < 4; v++ {
		for i := 0; i < 3; i++ {
			t[v][i] = (high - low) * rand.Float32() + low
		}
	}
}

func randomAnchoredPluckerVec() *AnchoredPluckerVec {
	x := rand.Float32()
	y := rand.Float32()
	z := rand.Float32()
	norm := float32(math.Sqrt(float64(x*x + y*y + z*z)))

	P := &Vec{0, 0, 0}
	L := &Vec{ x/norm, y/norm, z/norm }
	p := new(AnchoredPluckerVec)
	p.Init(P, L)
	return p
}


func BenchmarkPluckerTetraInit(b *testing.B) {
	ts := make([]Tetra, 1<<10)
	pts := make([]PluckerTetra, len(ts))
	for i := range ts {
		randomizeTetra(&ts[i], 0, 1)
	}

	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		i := n % len(ts)
		pts[i].Init(&ts[i])
	}
}

func BenchmarkIntersectionBary(b *testing.B) {
	ts := make([]Tetra, 1<<10)
	pts := make([]PluckerTetra, len(ts))
	for i := range ts {
		randomizeTetra(&ts[i], 0, 1)
		pts[i].Init(&ts[i])
	}

	w := new(IntersectionWorkspace)
	ap := randomAnchoredPluckerVec()

	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		i := n % len(ts)
		w.IntersectionBary(&pts[i], &ap.PluckerVec)
	}
}

func BenchmarkIntersectionDistance(b *testing.B) {
	ts := make([]Tetra, 1<<10)
	pts := make([]PluckerTetra, len(ts))
	for i := range ts {
		randomizeTetra(&ts[i], 0, 1)
		pts[i].Init(&ts[i])
	}

	w := new(IntersectionWorkspace)
	ap := randomAnchoredPluckerVec()

	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		i := n % len(ts)
		w.IntersectionDistance(&pts[i], &ts[i], ap)
	}
}

func almostEq(x1, x2, eps float32) bool {
	return x1 + eps > x2 && x1 - eps < x2
}

func TestIntersectionDistance(t *testing.T) {
	P, L := Vec{-0.5, 1, 0.5}, Vec{1, 0, 0}
	eps := float32(1e-4)

	ap := new(AnchoredPluckerVec)
	pt := new(PluckerTetra)
	w := new(IntersectionWorkspace)

	ap.Init(&P, &L)

	table := []struct{
		t Tetra
		enter, exit float32 
		ok bool
	}{	
		{Tetra{Vec{1,4,0},Vec{1,0,4},Vec{1,0,0},Vec{5,0,0}},1.5,2.5,true},
		{Tetra{Vec{1,0,4},Vec{1,4,0},Vec{1,0,0},Vec{5,0,0}},1.5,2.5,true},

		{Tetra{Vec{-1,4,0},Vec{-1,0,4},Vec{-1,0,0},Vec{3,0,0}},-0.5,0.5,true},
		{Tetra{Vec{9,4,0},Vec{9,0,4},Vec{9,0,0},Vec{13,0,0}},9.5,10.5,true},


		{Tetra{Vec{1,0,4},Vec{1,0,0},Vec{5,0,0},Vec{1,4,0}},1.5,2.5,true},
		{Tetra{Vec{1,0,0},Vec{5,0,0},Vec{1,4,0},Vec{1,0,4}},1.5,2.5,true},
		{Tetra{Vec{5,0,0},Vec{1,4,0},Vec{1,0,4},Vec{1,0,0},},1.5,2.5,true},
		
		{Tetra{Vec{1,6,0},Vec{1,2,4},Vec{1,2,0},Vec{5,2,0}},0,0,false},
	}

	for i, test := range table {
		test.t.Orient(+1)
		pt.Init(&test.t)
		enter, exit, ok := w.IntersectionDistance(pt, &test.t, ap)
		if ok != test.ok {
			t.Errorf("%d) Expected ok = %v, but got %v.", i+1, test.ok, ok)
		} else if !almostEq(enter, test.enter, eps) {
			t.Errorf(
				"%d) Expected enter = %g, but got %g", i+1, test.enter, enter,
			)
		} else if !almostEq(exit, test.exit, eps) {
			t.Errorf(
				"%d) Expected leave = %g, but got %g", i+1, test.exit, exit,
			)
		}
	}
}
