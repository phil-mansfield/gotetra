package analyze

import (
//	"math"
	"math/rand"
	"time"
	"fmt"
	"testing"
)

func sphere(r float64) Shell {
	return func(phi, theta float64) float64 {
		return r
	}
}

func TestEverything(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	s := sphere(1.0)
	samples := 10 * 1000
	fmt.Printf("Volume:  %8.3g\n", s.Volume(samples))
	Ix, Iy, Iz := s.Moments(samples)
	fmt.Printf("Moments: %8.3g %8.3g %8.3g\n", Ix, Iy, Iz)
	fmt.Printf("Area:    %8.3g\n", s.SurfaceArea(samples))
}
