package sort

import (
	"sort"
	"testing"
	"math/rand"
)

func sliceEq(xs, ys []float64) bool {
	if len(xs) != len(ys) {
		return false
	}
	for i := range xs {
		if xs[i] != ys[i] {
			return false
		}
	}

	return true
}

func randSlice(n int) []float64 {
	xs := make([]float64, n)
	for i := range xs {
		xs[i] = rand.Float64()
	}
	return xs
}

func TestReverse(t *testing.T) {
	if !sliceEq([]float64{1, 2, 3, 4, 5}, Reverse([]float64{5, 4, 3, 2, 1})) ||
		!sliceEq([]float64{2, 3, 4, 5}, Reverse([]float64{5, 4, 3, 2})) {
		t.Errorf("Welp, I hope you're proud of yourself.")
	}
}

func BenchmarkReverse10(b *testing.B) {
	xs := make([]float64, 10)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Reverse(xs)
	}
}

func BenchmarkReverse1000(b *testing.B) {
	xs := make([]float64, 1000)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Reverse(xs)
	}
}

func BenchmarkReverse1000000(b *testing.B) {
	xs := make([]float64, 1000000)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Reverse(xs)
	}
}

func BenchmarkShell10(b *testing.B) {
	xs := randSlice(10)
	buf := make([]float64, 10)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		Shell(buf)
	}
}

func BenchmarkShell100(b *testing.B) {
	xs := randSlice(100)
	buf := make([]float64, 100)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		Shell(buf)
	}
}

func BenchmarkShell1000(b *testing.B) {
	xs := randSlice(1000)
	buf := make([]float64, 1000)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		Shell(buf)
	}
}

func BenchmarkShell10000(b *testing.B) {
	xs := randSlice(10000)
	buf := make([]float64, 10000)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		Shell(buf)
	}
}

func BenchmarkQuick10(b *testing.B) {
	xs := randSlice(10)
	buf := make([]float64, 10)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		Quick(buf)
	}
}

func BenchmarkQuick100(b *testing.B) {
	xs := randSlice(100)
	buf := make([]float64, 100)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		Quick(buf)
	}
}

func BenchmarkQuick1000(b *testing.B) {
	xs := randSlice(1000)
	buf := make([]float64, 1000)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		Quick(buf)
	}
}

func BenchmarkQuick10000(b *testing.B) {
	xs := randSlice(10000)
	buf := make([]float64, 10000)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		Quick(buf)
	}
}

func BenchmarkGo10(b *testing.B) {
	xs := randSlice(10)
	buf := make([]float64, 10)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		sort.Float64s(buf)
	}
}

func BenchmarkGo100(b *testing.B) {
	xs := randSlice(100)
	buf := make([]float64, 100)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		sort.Float64s(buf)
	}
}

func BenchmarkGo1000(b *testing.B) {
	xs := randSlice(1000)
	buf := make([]float64, 1000)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		sort.Float64s(buf)
	}
}

func BenchmarkGo10000(b *testing.B) {
	xs := randSlice(10000)
	buf := make([]float64, 10000)
	for i := 0; i < b.N; i++ {
		copy(buf, xs)
		sort.Float64s(buf)
	}
}

func BenchmarkMedian10(b *testing.B) {
	xs := randSlice(10)
	buf := make([]float64, len(xs))
	n := len(xs) / 2
	for i := 0; i < b.N; i++ {
		NthLargest(xs, n, buf)
	}
}


func BenchmarkMedian100(b *testing.B) {
	xs := randSlice(100)
	buf := make([]float64, len(xs))
	n := len(xs) / 2
	for i := 0; i < b.N; i++ {
		NthLargest(xs, n, buf)
	}
}


func BenchmarkMedian1000(b *testing.B) {
	xs := randSlice(1000)
	buf := make([]float64, len(xs))
	n := len(xs) / 2
	for i := 0; i < b.N; i++ {
		NthLargest(xs, n, buf)
	}
}


func BenchmarkMedian10000(b *testing.B) {
	xs := randSlice(10000)
	buf := make([]float64, len(xs))
	n := len(xs) / 2
	for i := 0; i < b.N; i++ {
		NthLargest(xs, n, buf)
	}
}

// Tests

func TestShell(t *testing.T) {
	for i := 0; i < 10; i++ {
		xs := randSlice(1000)
		Shell(xs)
		if !sort.Float64sAreSorted(xs) {
			t.Errorf("Failed to sort.")
		}
	}
}

func TestQuick(t *testing.T) {
	for i := 0; i < 10; i++ {
		xs := randSlice(1000)
		Quick(xs)
		if !sort.Float64sAreSorted(xs) {
			t.Errorf("Failed to sort.")
		}
	}
}


func TestMedian(t *testing.T) {
	buf := make([]float64, 1000)
	for i := 0; i < 10; i++ {
		xs := randSlice(len(buf))
		Quick(xs)

		perm := rand.Perm(len(buf))
		mixed := make([]float64, len(buf))
		for j := range mixed {
			mixed[j] = xs[perm[j]]
		}
		
		for j := 1; j <= len(buf); j++ {
			val := NthLargest(mixed, j, buf)
			if val != xs[len(xs) - j] {
				t.Errorf("Failed to find NthLargest.")
			}
		}
	}
}
