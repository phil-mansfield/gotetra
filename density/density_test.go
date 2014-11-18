package density

import (
	"testing"

	"github.com/phil-mansfield/gotetra/geom"
)

func BenchmarkNGP(b *testing.B) {
	g, bg := Bounds(100, 1, 0, 0, 0)
	rhos := make([]float64, 100 * 100 * 100)
	pts := make([]geom.Vec, 1000)
	intr := NewInterpolator(NearestGridPoint, g, bg, 1.0, rhos)

	b.ResetTimer()

	for i:= 0; i < (b.N / len(pts)) + 1; i++ {
		intr.Interpolate(1, pts)
	}
}

func BenchmarkCIC(b *testing.B) {
	g, bg := Bounds(100, 1, 0, 0, 0)
	rhos := make([]float64, 100 * 100 * 100)
	pts := make([]geom.Vec, 1000)
	intr := NewInterpolator(CloudInCell, g, bg, 1.0, rhos)

	b.ResetTimer()

	for i:= 0; i < (b.N / len(pts)) + 1; i++ {
		intr.Interpolate(1, pts)
	}
}
