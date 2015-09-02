package analyze

import (
	"math"
	"math/rand"
	
	"github.com/phil-mansfield/gotetra/los"
)

type Shell func(phi, theta float64) float64
type ProjectedShell func(h *los.HaloProfiles, ring int, phi float64) float64

func randomAngle() (phi, theta float64) {
	u, v := rand.Float64(), rand.Float64()
	return 2 * math.Pi * u, math.Acos(2 * v - 1)
}

func cartesian(phi, theta, r float64) (x, y, z float64) {
	sinP, cosP := math.Sincos(phi)
	sinT, cosT := math.Sincos(theta)
	return r * sinT * cosP, r * sinT * sinP, r * cosT
}

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

func CumulativeShells(
	xs, ys [][]float64, h *los.HaloProfiles,
	I, J, start, stop, step int,
) (ringCounts []int, shells []Shell) {
	n := 0
	for i := range xs { n += len(xs[i]) }
	fXs, fYs, fZs := make([]float64, n), make([]float64, n), make([]float64, n)

	idx := 0
	for i := range xs {
		for j := range xs[i] {
			fXs[idx], fYs[idx], fZs[idx] =
				h.PlaneToVolume(i, xs[i][j], ys[i][j])
			idx++
		}
	}

	shells, ringCounts = []Shell{}, []int{}
	for rings := start; rings < stop; rings++ {
		end := 0
		for _, x := range xs[:rings] { end += len(x) }
		cs := PennaCoeffs(fXs[:end], fYs[:end], fZs[:end], I, J, 2)
		shells = append(shells, PennaFunc(cs, I, J, 2))
		ringCounts = append(ringCounts, rings)
	}

	return ringCounts, shells
}

type Tracers struct {
	Vol, Sa, Ix, Iy, Iz float64
}

func CumulativeTracers(
	shells [][]Shell, samples int,
) (means, stds []Tracers) {
	rings := len(shells[0])
	sums, sqrs := make([]Tracers, rings), make([]Tracers, rings)

	for ir := range shells[0] {
		for ih := range shells {
			shell := shells[ih][ir]

			vol := shell.Volume(samples)
			sa  := shell.SurfaceArea(samples)
			ix, iy, iz := shell.Moments(samples)

			sums[ir].Vol += vol
			sums[ir].Sa += sa
			sums[ir].Ix += ix
			sums[ir].Iy += iy
			sums[ir].Iz += iz

			sqrs[ir].Vol += vol*vol
			sqrs[ir].Sa += sa*sa
			sqrs[ir].Ix += ix*ix
			sqrs[ir].Iy += iy*iy
			sqrs[ir].Iz += iz*iz
		}
	}
	
	means, stds = make([]Tracers, rings), make([]Tracers, rings)
	n := float64(len(shells))
	for i := range means {
		means[i].Vol = sums[i].Vol / n
		means[i].Sa =  sums[i].Sa  / n
		means[i].Ix =  sums[i].Ix  / n
		means[i].Iy =  sums[i].Iy  / n
		means[i].Iz =  sums[i].Iz  / n

		stds[i].Vol = math.Sqrt(sqrs[i].Vol / n - means[i].Vol*means[i].Vol)
		stds[i].Sa =  math.Sqrt(sqrs[i].Sa  / n - means[i].Sa*means[i].Sa)
		stds[i].Ix =  math.Sqrt(sqrs[i].Ix  / n - means[i].Ix*means[i].Ix)
		stds[i].Iy =  math.Sqrt(sqrs[i].Iy  / n - means[i].Iy*means[i].Iy)
		stds[i].Iz =  math.Sqrt(sqrs[i].Iz  / n - means[i].Iz*means[i].Iz)
	}

	return means, stds
}
