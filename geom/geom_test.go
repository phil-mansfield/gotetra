package geom

import (
	"github.com/phil-mansfield/num/rand"
	"testing"
)

const (
	genType = rand.Tausworthe
)

func (t *Tetra) random(gen *rand.Generator, width float64) {
	c1, c2, c3, c4 := &Vec{}, &Vec{}, &Vec{}, &Vec{}
	c1.random(gen, width)
	c2.random(gen, width)
	c3.random(gen, width)
	c4.random(gen, width)

	t.Init(c1, c2, c3, c4, width)
}

func (v *Vec) random(gen *rand.Generator, width float64) {
	v[0] = float32(gen.Uniform(0, 1))
	v[1] = float32(gen.Uniform(0, 1))
	v[2] = float32(gen.Uniform(0, 1))
}

func BenchmarkVolume(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		t := &ts[i % len(ts)]
		t.volumeValid = false
		t.Volume()
	}
}

func BenchmarkBarycenter(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		t := &ts[i % len(ts)]
		t.baryValid = false
		t.Barycenter()
	}
}

func BenchmarkCellBounds(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}

	g := NewGrid(&[3]int{500, 500, 500}, 100)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		t := &ts[i % len(ts)]
		t.baryValid = false
		t.CellBounds(g)
	}
}

func BenchmarkCellBoundsAt(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}

	g := NewGrid(&[3]int{500, 500, 500}, 100)
	cb := &CellBounds{}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		t := &ts[i % len(ts)]
		t.baryValid = false
		t.CellBoundsAt(g, cb)
	}
}

func BenchmarkContains(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}

	vs := make([]Vec, len(ts) + 1)
	for i := range vs {
		vs[i].random(gen, 1.0)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		t := &ts[i % len(ts)]
		t.Contains(&vs[i % len(vs)])
	}
}

func BenchmarkSample1(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}
	
	randBuf := make([]float64, 1 * 3)
	vecBuf := make([]Vec, 1)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		t := &ts[i % len(ts)]
		t.baryValid = false
		t.Sample(gen, randBuf, vecBuf)
	}
}

func BenchmarkSample10(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}
	
	randBuf := make([]float64, 10 * 3)
	vecBuf := make([]Vec, 10)

	b.ResetTimer()
	for i := 0; i < b.N / 10 + 1; i++ {
		t := &ts[i % len(ts)]
		t.baryValid = false
		t.Sample(gen, randBuf, vecBuf)
	}
}

func BenchmarkSample100(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}
	
	randBuf := make([]float64, 100 * 3)
	vecBuf := make([]Vec, 100)

	b.ResetTimer()
	for i := 0; i < b.N / 100 + 1; i++ {
		t := &ts[i % len(ts)]
		t.baryValid = false
		t.Sample(gen, randBuf, vecBuf)
	}
}

func BenchmarkSample1000(b *testing.B) {
	ts := make([]Tetra, b.N / 20 + 1)
	gen := rand.NewTimeSeed(genType)
	for i := range ts {
		ts[i].random(gen, 1.0)
	}
	
	randBuf := make([]float64, 1000 * 3)
	vecBuf := make([]Vec, 1000)

	b.ResetTimer()
	for i := 0; i < b.N / 1000 + 1; i++ {
		t := &ts[i % len(ts)]
		t.baryValid = false
		t.Sample(gen, randBuf, vecBuf)
	}
}
