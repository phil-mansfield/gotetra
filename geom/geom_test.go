package geom

import (
	"testing"

	"github.com/phil-mansfield/num/rand"
)

const (
	genType = rand.Tausworthe
	testEps = 1e-6
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

func TestVolume(t *testing.T) {
	table := []struct{
		c1, c2, c3, c4 *Vec
		width, vol float64
		valid bool
	} {
		{&Vec{0, 0, 0}, &Vec{1, 0, 0}, &Vec{0, 1, 0}, &Vec{0, 0, 1},
			100, 1.0/6.0, true},
		{&Vec{0, 0, 0}, &Vec{1, 0, 0}, nil, &Vec{0, 0, 1},
			100, 1.0/6.0, false},
		{&Vec{0, 0, 0}, &Vec{1, 0, 0}, &Vec{0.5, 0, 0}, &Vec{0, 0, 1},
			100, 0, true},
	}

	for i, test := range table {
		tet, valid := NewTetra(test.c1, test.c2, test.c3, test.c4, test.width)
		
		if valid != test.valid {
			t.Errorf("%d) Tetra {%v %v %v %v} gives validity %b, not %b\n",
				i, test.c1, test.c2, test.c3, test.c4, valid, test.valid)
		} else if valid && !epsEq(test.vol, tet.Volume(), testEps) {
			t.Errorf("%d) Tetra {%v %v %v %v} gives volume %g, not %g\n",
				i, test.c1, test.c2, test.c3, test.c4, tet.Volume(), test.vol)
		}
	}
}

func TestContains(t *testing.T) {
	table := []struct {
		c1, c2, c3, c4, pt *Vec
		width float64
		res bool
	} {
		{&Vec{1, 1, 1}, &Vec{2, 1, 1}, &Vec{1, 2, 1}, &Vec{1, 1, 2},
			&Vec{1.1, 1.1, 1.1}, 3, true},
		{&Vec{1, 1, 1}, &Vec{2, 1, 1}, &Vec{1, 2, 1}, &Vec{1, 1, 2},
			&Vec{0.9, 1.1, 1.1}, 3, false},
		{&Vec{1, 1, 1}, &Vec{2, 1, 1}, &Vec{1, 2, 1}, &Vec{1, 1, 2},
			&Vec{1.1, 0.9, 1.1}, 3, false},
		{&Vec{1, 1, 1}, &Vec{2, 1, 1}, &Vec{1, 2, 1}, &Vec{1, 1, 2},
			&Vec{1.1, 1.1, 0.9}, 3, false},
		{&Vec{1, 1, 1}, &Vec{2, 1, 1}, &Vec{1, 2, 1}, &Vec{1, 1, 2},
			&Vec{1.5, 1.5, 1.5}, 3, false},	

		{&Vec{-.25, -.25, -.25}, &Vec{.75, -.25, -.25},
			&Vec{-.25, .75, -.25}, &Vec{-.25, -.25, .75},
			&Vec{-.2, -.2, -.2}, 3, true},	
		{&Vec{-.25, -.25, -.25}, &Vec{.75, -.25, -.25},
			&Vec{-.25, .75, -.25}, &Vec{-.25, -.25, .75},
			&Vec{2.8, 2.8, 2.8}, 3, true},	
	}

	for i, test := range table {
		tet, _ := NewTetra(test.c1, test.c2, test.c3, test.c4, test.width)
		
		if tet.Contains(test.pt) != test.res {
			t.Errorf("%d) [%v %v %v %v, %g].Contains(%v) = %v, not %v \n",
				i, test.c1, test.c2, test.c3, test.c4, test.width,
				test.pt, tet.Contains(test.pt), test.res,
			)
		}
	}
}

func TestContainsMC(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping TestContainsMC")
	}
	
	xs := make([]float64, 10000)
	ys := make([]float64, 10000)
	zs := make([]float64, 10000)
	
	gen := rand.NewTimeSeed(genType)
	tet, _ := NewTetra(
		&Vec{2.8, 2.8, 2.8}, &Vec{0.8, 2.8, 2.8}, 
		&Vec{2.8, 0.8, 2.8}, &Vec{2.8, 2.8, 0.8}, 3,
	)
	v := &Vec{}

	inside, total := 0, 0
	for i := 0; i < 100; i++ {
		gen.UniformAt(0.0, 3.0, xs)
		gen.UniformAt(0.0, 3.0, ys)
		gen.UniformAt(0.0, 3.0, zs)

		for j := range xs {
			v[0], v[1], v[2] = float32(xs[j]), float32(ys[j]), float32(zs[j])

			if tet.Contains(v) {
				inside++
			}
			total++
		}
	}

	exp := float64(total) / (27 * 6)

	if !epsEq(exp, float64(inside), 0.05) {
		normFrac := float64(inside) / float64(total) * (27 * 6)
		t.Errorf("%d million point MC integration of tetrahedron volume " +
			"gives %f/27 instead of 1/27.", total / 1000000, normFrac)
	}
}
