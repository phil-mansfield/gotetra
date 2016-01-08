package sphere_halo

import (
	"math/rand"
	"testing"

	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/los/geom"
//	"github.com/phil-mansfield/gotetra/los"
)

func BenchmarkTransform10000000(b *testing.B) {
	vecs := make([]rgeom.Vec, 10 * 1000 * 1000)
	for i := range vecs {
		vecs[i][0] = float32(rand.Float64())
		vecs[i][1] = float32(rand.Float64())
		vecs[i][2] = float32(rand.Float64())
	}

	h := SphereHalo{}
	h.Init(nil, [3]float64{0.25, 0.25, 0.25}, 0, 0, 0, 0)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		h.Transform(vecs, 0.5)
	}
}

func BenchmarkIntersect10000000(b *testing.B) {
	vecs := make([]rgeom.Vec, 10 * 1000 * 1000)
	intr := make([]bool, 10 * 1000 * 1000)
	for i := range vecs {
		vecs[i][0] = float32(rand.Float64())
		vecs[i][1] = float32(rand.Float64())
		vecs[i][2] = float32(rand.Float64())
	}

	h := SphereHalo{}
	h.Init(nil, [3]float64{0.5, 0.5, 0.5}, 0.25, 0.5, 0, 0)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		h.Intersect(vecs, 0.1, intr)
	}
}

func BenchmarkSplitJoin16(b *testing.B) {
	norms := make([]geom.Vec, 100)
	for i := range norms { norms[i] = geom.Vec{0, 0, 1} }
	h := SphereHalo{}
	h.Init(norms, [3]float64{1, 1, 1}, 0.5, 5.0, 200, 256)
	hs := make([]SphereHalo, 15)
	h.Split(hs)
	
	for i := 0; i < b.N; i++ {
		h.Split(hs)
		h.Join(hs)
	}
}

func BenchmarkSplit16(b *testing.B) {
	norms := make([]geom.Vec, 100)
	for i := range norms { norms[i] = geom.Vec{0, 0, 1} }
	h := SphereHalo{}
	h.Init(norms, [3]float64{1, 1, 1}, 0.5, 5.0, 200, 256)
	hs := make([]SphereHalo, 15)
	h.Split(hs)

	for i := 0; i < b.N; i++ {
		h.Split(hs)
	}
}

func BenchmarkJoin16(b *testing.B) {
	norms := make([]geom.Vec, 100)
	for i := range norms { norms[i] = geom.Vec{0, 0, 1} }
	h := SphereHalo{}
	h.Init(norms, [3]float64{1, 1, 1}, 0.5, 5.0, 200, 256)
	hs := make([]SphereHalo, 15)
	h.Split(hs)

	for i := 0; i < b.N; i++ {
		h.Join(hs)
	}
}
