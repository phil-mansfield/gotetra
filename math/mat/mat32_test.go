package mat

import (
	"testing"
	"math/rand"
)

func randomMatrices32(w, h, n int) []*Matrix32 {
	ms := make([]*Matrix32, n)
	for i := range ms {
		ms[i] = NewMatrix32(make([]float32, w*h), w, h)
		for j := range ms[i].Vals { ms[i].Vals[j] = rand.Float32() }
	}
	return ms
}

func Benchmark32Mult3(b *testing.B) {
	n := 10
	m1s := randomMatrices32(3, 3, n)
	m2s := randomMatrices32(3, 3, n)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		m1s[idx].Mult(m2s[idx])
	}
}

func Benchmark32MultAt3(b *testing.B) {
	n := 10
	m1s := randomMatrices32(3, 3, n)
	m2s := randomMatrices32(3, 3, n)
	out := NewMatrix32(make([]float32, 9), 3, 3)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		m1s[idx].MultAt(m2s[idx], out)
	}
}

func Benchmark32Mult1024(b *testing.B) {
	n := 10
	m1s := randomMatrices32(1<<10, 1<<10, n)
	m2s := randomMatrices32(1<<10, 1<<10, n)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		m1s[idx].Mult(m2s[idx])
	}
}

func Benchmark32LU3(b *testing.B) {
	n := 10
	ms := randomMatrices32(3, 3, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		ms[idx].LU()
	}
}

func Benchmark32LUFactorsAt3(b *testing.B) {
	n := 10
	ms := randomMatrices32(3, 3, n)
	lu := NewLUFactors32(3)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		ms[idx].LUFactorsAt(lu)
	}
}

func Benchmark32Determinant3(b *testing.B) {
	n := 10
	ms := randomMatrices32(3, 3, n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		ms[idx].Determinant()
	}
}

func Benchmark32LUDeterminant3(b *testing.B) {
	n := 10
	ms := randomMatrices32(3, 3, n)
	lus := make([]*LUFactors32, n)
	for i := range lus {
		lus[i] = NewLUFactors32(3)
		ms[i].LUFactorsAt(lus[i])
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		lus[idx].Determinant()
	}
}

func Benchmark32Invert3(b *testing.B) {
	n := 10
	ms := randomMatrices32(3, 3, n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		ms[idx].Invert()
	}
}

func Benchmark32LUInvert3(b *testing.B) {
	n := 10
	ms := randomMatrices32(3, 3, n)
	lus := make([]*LUFactors32, n)
	for i := range lus {
		lus[i] = NewLUFactors32(3)
		ms[i].LUFactorsAt(lus[i])
	}
	out := NewMatrix32(make([]float32, 9), 3, 3)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		lus[idx].InvertAt(out)
	}
}
