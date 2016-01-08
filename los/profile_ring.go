/*package los computes line of sight densities oriented in rings around a point.
*/
package los

import (
	"math"

	"github.com/phil-mansfield/gotetra/los/geom"
)

// ProfileRing is a ring of uniformly spaced profiles. The reported profiles are
// the average densities at each radial bin in the profile.
//
// The saved format of the profiles is not neccessarily the same as the output
// format of the profiles, so they must be accessed with the Retreive method.
type ProfileRing struct {
	derivs []float64 // Contiguous block of pofile data. Column-major.
	Lines []geom.Line
	bins int // Length of an individual profile.
	n int // Number of profiles.
	lowR, highR, dr float64
}

func (p *ProfileRing) Reuse(lowR, highR float64) {
	p.lowR = lowR
	p.highR = highR
	p.dr = (highR - lowR) / float64(p.bins)	
}

func (p1 *ProfileRing) Join(p2 *ProfileRing) {
	if p1.n != p2.n || p1.bins != p2.bins {
		panic("ProfileRing sizes do not match.")
	} else {
		for i, val := range p2.derivs { p1.derivs[i] += val }
	}
}

func (p1 *ProfileRing) Split(p2 *ProfileRing) {
	if p1.n != p2.n || p1.bins != p2.bins {
		p2.Init(p1.lowR, p1.highR, p1.bins, p1.n)
	} else {
		p2.bins = p1.bins
		p2.lowR = p1.lowR
		p2.highR = p1.highR
		p2.dr = p2.dr

		for i := range p2.derivs { p2.derivs[i] = 0 }
	}
}

// Init initializes a profile ring made up of n profiles each of which consist
// of the given number of radial bins and extend between the two specified 
// radii.
func (p *ProfileRing) Init(lowR, highR float64, bins, n int) {
	p.derivs = make([]float64, bins * n)
	p.bins = bins
	p.n = n
	p.lowR = lowR
	p.highR = highR
	p.dr = (highR - lowR) / float64(bins)
	
	p.Lines = make([]geom.Line, n)
	for i := 0; i < n; i++ {
		sin, cos := math.Sincos(p.Angle(i))
		p.Lines[i].Init(0, 0, float32(cos), float32(sin))
	}
}

// Insert inserts a plateau with the given radial extent and density to the
// profile.
func (p *ProfileRing) Insert(start, end, rho float64, i int) {
	if end <= p.lowR || start >= p.highR { return }

	// You could be a bit more careful with floating point ops here.
	if start > p.lowR {
		fidx, rem := math.Modf((start - p.lowR) / p.dr)
		idx := int(fidx)
		p.derivs[i*p.bins + idx] += rho * (1 - rem)
		if idx < p.bins - 1 {
			p.derivs[i*p.bins + idx+1] += rho * rem
		}
	} else {
		p.derivs[i*p.bins] += rho
	}

	if end < p.highR {
		fidx, rem := math.Modf((end - p.lowR) / p.dr)
		idx := int(fidx)
		p.derivs[i*p.bins + idx] -= rho * (1 - rem)
		if idx < p.bins - 1 {
			p.derivs[i*p.bins + idx+1] -= rho * rem
		}
	}
}

// Retrieve does any neccessary post-processing on the specified profile and
// writes in to an out buffer.
func (p *ProfileRing) Retrieve(i int, out []float64) {
	sum := float64(0)
	for j := 0; j < p.bins; j++ {
		sum += p.derivs[j + p.bins*i]
		out[j] = sum
	}
}

// Angle returns the angle that the line with the specified index points in.
func (p *ProfileRing) Angle(i int) float64 {
	return math.Pi * 2 * float64(i) / float64(p.n)
}
