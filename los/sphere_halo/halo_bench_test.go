package sphere_halo

import (
	"math"
	"math/rand"
	"testing"

	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/los/geom"
)

func BenchmarkTransform10000000(b *testing.B) {
	vecs := make([]rgeom.Vec, 10 * 1000 * 1000)
	for i := range vecs {
		vecs[i][0] = float32(rand.Float64())
		vecs[i][1] = float32(rand.Float64())
		vecs[i][2] = float32(rand.Float64())
	}

	h := SphereHalo{}
	h.Init(nil, [3]float64{0.25, 0.25, 0.25}, 0, 0, 0, 0, 0)
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
	h.Init(nil, [3]float64{0.5, 0.5, 0.5}, 0.25, 0.5, 0, 0, 0)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		h.Intersect(vecs, 0.1, intr)
	}
}

func BenchmarkSplitJoin16(b *testing.B) {
	norms := make([]geom.Vec, 100)
	for i := range norms { norms[i] = geom.Vec{0, 0, 1} }
	h := SphereHalo{}
	h.Init(norms, [3]float64{1, 1, 1}, 0.5, 5.0, 200, 256, 0)
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
	h.Init(norms, [3]float64{1, 1, 1}, 0.5, 5.0, 200, 256, 0)
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
	h.Init(norms, [3]float64{1, 1, 1}, 0.5, 5.0, 200, 256, 0)
	hs := make([]SphereHalo, 15)
	h.Split(hs)

	for i := 0; i < b.N; i++ {
		h.Join(hs)
	}
}

func BenchmarkGetRhos(b *testing.B) {
	h := SphereHalo{}
	norms := make([]geom.Vec, 100)
	for i := range norms { norms[i] = geom.Vec{0, 0, 1} }
	h.Init(norms, [3]float64{0, 0, 0}, 0.5, 5.0, 200, 256, 0)
	buf := make([]float64, 200)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		h.GetRhos(50, 100, buf)
	}
}

func BenchmarkGetRhosFull(b *testing.B) {
	h := SphereHalo{}
	norms := make([]geom.Vec, 100)
	for i := range norms { norms[i] = geom.Vec{0, 0, 1} }
	h.Init(norms, [3]float64{0, 0, 0}, 0.5, 5.0, 200, 256, 0)
	bufs := make([][][]float64, 100)
	for i := range bufs {
		bufs[i] = make([][]float64, 256)
		for j := range bufs[i] {
			bufs[i][j] = make([]float64, 200)
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for ring := 0; ring < 100; ring++ {
			for los := 0; los < 256; los++ {
				h.GetRhos(ring, los, bufs[ring][los])
			}
		}
	}
}

func randomDirs(n int) []geom.Vec {
	vecs := make([]geom.Vec, n)
	for i := range vecs {
		for {
			x := rand.Float64()*2 - 1
			y := rand.Float64()*2 - 1
			z := rand.Float64()*2 - 1
			r := math.Sqrt(x*x + y*y + z*z)
			if r > 1 { continue }
			vecs[i] = geom.Vec{
				float32(x / r),
				float32(y / r),
				float32(z / r),
			}
			break
		}
	}
	return vecs
}

func BenchmarkInsert1(b *testing.B) {
	h := SphereHalo{}
	h.Init(randomDirs(100), [3]float64{0, 0, 0}, 0.3, 3, 200, 256, 0)

	var vecR float32 = 1.0
	vecs := randomDirs(10000)
	for i := range vecs {
		vecs[i][0] *= vecR
		vecs[i][1] *= vecR
		vecs[i][2] *= vecR
	}

	sphR := 0.1
	idx := 0

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		h.Insert(vecs[idx], sphR, 1)

		idx++
		if idx == 10000 { idx = 0 }
	}
}

func BenchmarkInsert3(b *testing.B) {
	h := SphereHalo{}
	h.Init(randomDirs(100), [3]float64{0, 0, 0}, 0.3, 3, 200, 256, 0)

	var vecR float32 = 3.0
	vecs := randomDirs(10000)
	for i := range vecs {
		vecs[i][0] *= vecR
		vecs[i][1] *= vecR
		vecs[i][2] *= vecR
	}

	sphR := 0.1
	idx := 0

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		h.Insert(vecs[idx], sphR, 1)

		idx++
		if idx == 10000 { idx = 0 }
	}
}

func BenchmarkInsert0_3(b *testing.B) {
	h := SphereHalo{}
	h.Init(randomDirs(100), [3]float64{0, 0, 0}, 0.3, 3, 200, 256, 0)

	var vecR float32 = 0.3
	vecs := randomDirs(10000)
	for i := range vecs {
		vecs[i][0] *= vecR
		vecs[i][1] *= vecR
		vecs[i][2] *= vecR
	}

	sphR := 0.1
	idx := 0

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		h.Insert(vecs[idx], sphR, 1)

		idx++
		if idx == 10000 { idx = 0 }
	}
}
