package mat

import (
	"testing"
	"math/rand"
)

func randomMatrices(w, h, n int) []*Matrix {
	ms := make([]*Matrix, n)
	for i := range ms {
		ms[i] = NewMatrix(make([]float64, w*h), w, h)
		for j := range ms[i].Vals { ms[i].Vals[j] = rand.Float64() }
	}
	return ms
}

func BenchmarkMult3(b *testing.B) {
	n := 10
	m1s := randomMatrices(3, 3, n)
	m2s := randomMatrices(3, 3, n)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		m1s[idx].Mult(m2s[idx])
	}
}

func BenchmarkMultAt3(b *testing.B) {
	n := 10
	m1s := randomMatrices(3, 3, n)
	m2s := randomMatrices(3, 3, n)
	out := NewMatrix(make([]float64, 9), 3, 3)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		m1s[idx].MultAt(m2s[idx], out)
	}
}

func BenchmarkMult1024(b *testing.B) {
	n := 10
	m1s := randomMatrices(1<<10, 1<<10, n)
	m2s := randomMatrices(1<<10, 1<<10, n)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		m1s[idx].Mult(m2s[idx])
	}
}

func BenchmarkLU3(b *testing.B) {
	n := 10
	ms := randomMatrices(3, 3, n)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		ms[idx].LU()
	}
}

func BenchmarkLUFactorsAt3(b *testing.B) {
	n := 10
	ms := randomMatrices(3, 3, n)
	lu := NewLUFactors(3)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		ms[idx].LUFactorsAt(lu)
	}
}

func BenchmarkDeterminant3(b *testing.B) {
	n := 10
	ms := randomMatrices(3, 3, n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		ms[idx].Determinant()
	}
}

func BenchmarkLUDeterminant3(b *testing.B) {
	n := 10
	ms := randomMatrices(3, 3, n)
	lus := make([]*LUFactors, n)
	for i := range lus {
		lus[i] = NewLUFactors(3)
		ms[i].LUFactorsAt(lus[i])
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		lus[idx].Determinant()
	}
}

func BenchmarkInvert3(b *testing.B) {
	n := 10
	ms := randomMatrices(3, 3, n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		ms[idx].Invert()
	}
}

func BenchmarkLUInvert3(b *testing.B) {
	n := 10
	ms := randomMatrices(3, 3, n)
	lus := make([]*LUFactors, n)
	for i := range lus {
		lus[i] = NewLUFactors(3)
		ms[i].LUFactorsAt(lus[i])
	}
	out := NewMatrix(make([]float64, 9), 3, 3)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		idx := i % n
		lus[idx].InvertAt(out)
	}
}
