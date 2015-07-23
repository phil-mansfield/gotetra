package geom

import (
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

