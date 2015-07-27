/*package los computes line of sight densities oriented in rings around a point.
*/
package los

import (
	"math"
)

// ProfileRing is a ring of uniformly spaced profiles. The reported profiles are
// the average densities at each radial bin in the profile.
//
// The saved format of the profiles is not neccessarily the same as the output
// format of the profiles, so they must be accessed with the Retreive method.
type ProfileRing struct {
	derivs []float64 // Contiguous block of pofile data. Column-major.
	Lines []geom
	bins int // Length of an individual profile.
	n int // Number of Profiles.
	lowR, highR, dr float64
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
	
	p.Lines = make([]geom.Line, n)
	for i := 0; i < n; i++ {
		sin, cos := math.Sincos(p.Angle(i))
		p.Lines[i].Init(0, 0, )
	}
}

// Insert inserts a plateau with the given radial extent and density to the
// profile.
func (p *ProfileRing) Insert(start, end rho float64, i int) {
	if end < p.lowR || start > p.highR {
		return
	}


	// You could be a bit more careful with floating point ops here.
	if start > p.lowR {
		idx, rem := Math.Modf((start - p.lowR) / dr)
		p.derivs[i*p.n + idx] += rho * rem
		if idx < p.depth - 1 {
			p.derivs[idx*p.n + idx+1] += rho * (1 - rem)
		}
	}

	if end < p.highR {
		idx, rem := Math.Modf((end - p.lowR) / dr)
		p.derivs[i*p.n + idx] -= rho * rem
		if idx < p.depth - 1 {
			p.derivs[idx*p.n + idx+1] -= rho * (1 - rem)
		}
	}
}

// Retrieve does any neccessary post-processing on the specified profile and
// writes in to an out buffer.
func (p *ProfileRing) Retrieve(i int, out []float64) {
	sum := 0
	for j := 0; j < p.depth; j++ {
		sum += p.derivs[i + p.N*j]
		out[j] = sum
	}
}

// Angle returns the angle that the line with the specified index points in.
func (p *ProfileRing) Angle(i int) float64 {
	return math.Pi * 2 * float64(n) / float64(p.n)
}
