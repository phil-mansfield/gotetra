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


func BenchmarkNaiveIntersection(b *testing.B) {
	targetIntersectionMode = naiveMode

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
		w.Intersection(&pts[i], &ap.PluckerVec)
	}
}
