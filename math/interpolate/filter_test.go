package interpolate

import (
	"testing"
)

func BenchmarkConvolveArray200Filter21(b *testing.B) {
	out, xs := make([]float64, 200), make([]float64, 200)
	k := NewTophatKernel(21)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		k.ConvolveAt(xs, Extension, out)
	}
}
