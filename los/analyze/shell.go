package analyze

import (
//	"fmt"
	"math"
	"math/rand"
)

func randomAngle() (phi, theta float64) {
	u, v := rand.Float64(), rand.Float64()
	return 2 * math.Pi * u, math.Acos(2 * v - 1)
}

func cartesian(phi, theta, r float64) (x, y, z float64) {
	sinP, cosP := math.Sincos(phi)
	sinT, cosT := math.Sincos(theta)
	return r * sinT * cosP, r * sinT * sinP, r * cosT
}

type Shell func(phi, theta float64) float64

func (s Shell) Volume(samples int) float64 {
	sum := 0.0
	for i := 0; i < samples; i++ {
		phi, theta := randomAngle()
		r := s(phi, theta)
		sum += r*r*r
	}
	r3 := sum / float64(samples)
	return r3 * 4 * (math.Pi / 3)
}

func (s Shell) Moments(samples int) (Ix, Iy, Iz float64) {
	xSum, ySum, zSum, rSum := 0.0, 0.0, 0.0, 0.0
	for i := 0; i < samples; i++ {
		phi, theta := randomAngle()
		r := s(phi, theta)
		x, y, z := cartesian(phi, theta, r)
		xSum += (y*y + z*z) * r*r
		ySum += (x*x + z*z) * r*r
		zSum += (x*x + y*y) * r*r
		rSum += r*r
	}
	return xSum / rSum, ySum / rSum, zSum / rSum
}

func (s Shell) SurfaceArea(samples int) float64 {
	sum := 0.0
	for i := 0; i < samples; i++ {
		phi, theta := randomAngle()
		r := s(phi, theta)
		sum += r*r
	}
	return sum / float64(samples) * 4 * math.Pi
}

func (s1 Shell) DiffVolume(s2 Shell, samples int) float64 {
	sum := 0.0
	for i := 0; i < samples; i++ {
		phi, theta := randomAngle()
		r1, r2 := s1(phi, theta), s2(phi, theta)
		r := (r1 + r2) / 2
		dr := math.Abs(r1 - r2)
		sum += dr*r*r
	}
	return sum / float64(samples) * (4 * math.Pi) / 3
}

func (s1 Shell) MaxDiff(s2 Shell, samples int) float64 {
	max := 0.0
	for i := 0; i < samples; i++ {
		phi, theta := randomAngle()
		r1, r2 := s1(phi, theta), s2(phi, theta)
		dr := math.Abs(r1 - r2)
		if dr > max { max = dr }
	}
	return max
}
