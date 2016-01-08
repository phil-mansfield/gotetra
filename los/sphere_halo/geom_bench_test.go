package sphere_halo

import (
	"math/rand"
	"testing"

	"github.com/phil-mansfield/gotetra/los/geom"
)

func BenchmarkHalfAngularWidth(b *testing.B) {
	for i := 0; i < b.N; i++ {
		halfAngularWidth(10, 2)
	}
}

func BenchmarkTwoValIntrDist(b *testing.B) {
	for i := 0; i < b.N; i++ {
		twoValIntrDist(10, 2, 1)
	}
}

func BenchmarkOneValIntrDist(b *testing.B) {
	for i := 0; i < b.N; i++ {
		oneValIntrDist(10, 2, 1, -1)
	}
}

func randomNorms(n int) []geom.Vec {
	norms := make([]geom.Vec, n)
	for i := range norms {
		x := float32(rand.Float64()*2 - 1)
		y := float32(rand.Float64()*2 - 1)
		z := float32(rand.Float64()*2 - 1)
		sum := x + y + z
		if sum == 0 {
			norms[i] = geom.Vec{0, 0, 1}
		} else {
			norms[i] = geom.Vec{x, y, z}
		}
	}
	return norms
}

func BenchmarkSphereIntersectHalo100(b *testing.B) {
	rings := 100
	norms := randomNorms(rings)
	h := SphereHalo{}
	h.Init(norms, [3]float64{0, 0, 0}, 1, 2, 1, 1)

	b.ResetTimer()
	v := geom.Vec{0, 0, 0.5}
	for i := 0; i < b.N; i++ {
		for ring := 0; ring < rings; ring++ {
			h.sphereIntersectRing(v, 0.1, ring)
		}
	}
}
