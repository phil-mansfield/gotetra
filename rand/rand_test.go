package rand

import (
	"testing"
)

func benchmarkUniform(gt GeneratorType, b *testing.B) {
	gen := NewTimeSeed(gt)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		_ = gen.Uniform(0, 13)
	}
}

func benchmarkUniformAt(gt GeneratorType, tLen int, b *testing.B) {
	gen := NewTimeSeed(gt)
	b.ResetTimer()

	target := make([]float64, tLen)

	n := 0
	for n < b.N {
		if n + tLen > b.N { target = target[0: b.N - n] }
		gen.UniformAt(0, 13, target)
		n += tLen
	}
}



func BenchmarkUniformGolang(b *testing.B) { benchmarkUniform(Golang, b) }
func BenchmarkUniformXorshift(b *testing.B) { benchmarkUniform(Xorshift, b) }
func BenchmarkUniformTausworthe(b *testing.B) { benchmarkUniform(Tausworthe, b) }

func BenchmarkUniformAtGolang(b *testing.B) { benchmarkUniformAt(Golang, DefaultBufSize, b) }
func BenchmarkUniformAtXorshift(b *testing.B) { benchmarkUniformAt(Xorshift, DefaultBufSize, b) }
func BenchmarkUniformAtTausworthe(b *testing.B) { benchmarkUniformAt(Tausworthe, DefaultBufSize, b) }

func BenchmarkNewSobolSequence(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = NewSobolSequence()
	}
}

func BenchmarkNextSobol6Sequence(b *testing.B) {
	seq := NewSobolSequence()
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		_ = seq.Next(6)
	}
}

func BenchmarkNextAtSobol6Sequence(b *testing.B) {
	seq := NewSobolSequence()
	vec := make([]float64, 6)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		seq.NextAt(vec)
	}
}
