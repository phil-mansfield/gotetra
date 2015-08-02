package los

import (
	"fmt"
	"math"

	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/math/mat"
)

type haloRing struct {
	ProfileRing
	phis []float32
	rot, irot mat.Matrix32
	rMax float64

	// Workspace objects.
	dr geom.Vec
	t geom.Tetra
	pt geom.PluckerTetra
	poly geom.TetraSlice
}

// Init initialized a haloRing.
func (hr *haloRing) Init(
	norm, origin *geom.Vec, rMin, rMax float64, bins, n int,
) {
	hr.ProfileRing.Init(rMin, rMax, bins, n)
	zAxis := &geom.Vec{0, 0, 1}
	
	hr.rot.Init(make([]float32, 9), 3, 3)
	hr.irot.Init(make([]float32, 9), 3, 3)
	geom.EulerMatrixBetweenAt(norm, zAxis, &hr.rot)
	geom.EulerMatrixBetweenAt(zAxis, norm, &hr.irot)
	hr.phis = make([]float32, n)
	for i := 0; i < n; i++ {
		hr.phis[i] = float32(i) / float32(n) * (2 * math.Pi)
	}

	hr.rMax = rMax
	hr.dr = *origin
	for i := 0; i < 3; i++ { hr.dr[i] *= -1 }
}

func (hr *haloRing) insert(t *geom.Tetra, rho float64) {
	hr.t = *t
	hr.t.Translate(&hr.dr)
	hr.t.Rotate(&hr.rot)
	hr.pt.Init(&hr.t) // This is slower than it has to be by a lot!
	
	rSqrMin, rSqrMax := float32(hr.lowR*hr.lowR), float32(hr.highR*hr.highR)
	if hr.t.ZPlaneSlice(&hr.pt, 0, &hr.poly) {
		// Stop early if possible.
		rSqrTMin, rSqrTMax := hr.poly.RSqrMinMax()
		if rSqrTMin > rSqrMax || rSqrTMax < rSqrMin { return }

		// Find the intersected LoS range and check each line in it for
		// intersection distance.
		lowPhi, phiWidth := hr.poly.AngleRange()
		lowIdx, idxWidth := geom.AngleBinRange(lowPhi, phiWidth, hr.n)
		
		idx := lowIdx

		for i := 0; i < idxWidth; i++ {
			if idx >= hr.n { idx -= hr.n }
			l := &hr.Lines[idx]
			l1, l2 := hr.poly.IntersectingLines(hr.phis[idx])

			// No intersections. Happens sometimes due to floating points
			/// fuzziness.
			if l1 == nil { continue }

			var rEnter, rExit float64
			if l2 != nil {
				// The polygon does not enclose the origin.
				enterX, enterY, _ := geom.Solve(l, l1)
				exitX, exitY, _ := geom.Solve(l, l2)
				
				rSqrEnter := enterX*enterX + enterY*enterY
				rSqrExit := exitX*exitX + exitY*exitY
				
				rEnter = math.Sqrt(float64(rSqrEnter))
				rExit = math.Sqrt(float64(rSqrExit))
				if rExit < rEnter { rEnter, rExit = rExit, rEnter }
			} else {
				// The polygon encloses the origin.
				exitX, exitY, _ := geom.Solve(l, l1)
				rSqrExit := exitX*exitX + exitY*exitY
				rExit = math.Sqrt(float64(rSqrExit))
				rEnter = 0.0
			}

			hr.Insert(rEnter, rExit, rho, idx)

			idx++
		}
	}
}

// Count inserts a tetrahedron into the haloRing so that its profiles represent
// overlap counts.
func (hr *haloRing) Count(t *geom.Tetra) { hr.insert(t, 1) }

// Density inserts a tetrahedorn into the haloRing so that its profiles
// represent densties.
func (hr *haloRing) Density(t *geom.Tetra, rho float64) { hr.insert(t, rho) }

// Add adds the contents of hr2 to hr 1.
func (hr1 *haloRing) Add(hr2 *haloRing) {
	for i, x := range hr2.derivs { hr1.derivs[i] += x }
}

// Clear resets the contents of the haloRing.
func (hr *haloRing) Clear() {
	for i := range hr.derivs { hr.derivs[i] = 0 }
}

// Add adds the contents of hp2 to hp1.
func (hp1 *HaloProfiles) Add(hp2 *HaloProfiles) {
	for i := range hp1.rs { hp1.rs[i].Add(&hp2.rs[i]) }
}

// Clear resets the conents of the HaloProfiles.
func (hp *HaloProfiles) Clear() {
	for i := range hp.rs { hp.rs[i].Clear() }
}

// HaloProfiles is a terribly-named struct which represents a halo and all its
// LoS profiles.
type HaloProfiles struct {
	geom.Sphere
	cCopy geom.Vec
	minSphere geom.Sphere

	rs []haloRing
	rMin, rMax float64
	id, bins, n int
}

// Init initializes a HaloProfiles struct with the given parameters.
func (hp *HaloProfiles) Init(
	id, rings int, origin *geom.Vec, rMin, rMax float64,
	bins, n int,
) *HaloProfiles {
	// We might be able to do better than this.
	var norms []geom.Vec
	if rings >= 3 {
		solid, ok := geom.NewUniquePlatonicSolid(rings)
		norms = solid.UniqueNormals()
		if !ok {
			panic(fmt.Sprintf("Cannot uniformly space %d rings.", rings))
		}
	} else if rings == 1 {
		norms = []geom.Vec{{0, 0, 1}}
	} else if rings == 2 {
		norms = []geom.Vec{{0, 0, 1}, {0, 1, 0}}
	}

	hp.rs = make([]haloRing, rings)
	hp.C, hp.R = *origin, float32(rMax)
	hp.cCopy = hp.C
	hp.minSphere.C, hp.minSphere.R = *origin, float32(rMin)
	hp.rMin, hp.rMax = rMin, rMax
	hp.id, hp.bins, hp.n = id, bins, n

	for i := 0; i < rings; i++ {
		hp.rs[i].Init(&norms[i], origin, rMin, rMax, bins, n)
	}

	return hp
}

// Count inserts the given tetrahedron such that the resulting profiles give
// tetrahedron overlap counts.
func (hp *HaloProfiles) Count(t *geom.Tetra) {
	for i := range hp.rs { hp.rs[i].Count(t) }
}

// Density inserts a tetrahedron such that the resulting profiles give
// densities.
func (hp *HaloProfiles) Density(t *geom.Tetra, rho float64) {
	for i := range hp.rs { hp.rs[i].Density(t, rho) }
}

func wrapDist(x1, x2, width float32) float32 {
	dist := x1 - x2
	if dist > width / 2 {
		return dist - width
	} else if dist < width / -2 {
		return dist + width
	} else {
		return dist
	}
}

func inRange(x, r, low, width, tw float32) bool {
	return wrapDist(x, low, tw) > -r && wrapDist(x, low + width, tw) < r
}

// SheetIntersect returns true if the given halo and sheet intersect one another
// and false otherwise.
func (hp *HaloProfiles) SheetIntersect(hd *io.SheetHeader) bool {
	tw := float32(hd.TotalWidth)
	return inRange(hp.C[0], hp.R, hd.Origin[0], hd.Width[0], tw) &&
		inRange(hp.C[1], hp.R, hd.Origin[1], hd.Width[1], tw) &&
		inRange(hp.C[2], hp.R, hd.Origin[2], hd.Width[2], tw)
}

// SphereIntersect returns true if the given halo and sphere intersect and false
// otherwise.
func (hp *HaloProfiles) SphereIntersect(s *geom.Sphere) bool {
	return hp.Sphere.SphereIntersect(s) && !hp.minSphere.SphereContain(s)
}

// VecIntersect returns true if the given vector is contained in the given halo
// and false otherwise.
func (hp *HaloProfiles) VecIntersect(v *geom.Vec) bool {
	return hp.Sphere.VecIntersect(v) && !hp.minSphere.VecIntersect(v)
}

// TetraIntersect returns true if the given vector and tetrahedron overlap.
func (hp *HaloProfiles) TetraIntersect(t *geom.Tetra) bool {
	return hp.Sphere.TetraIntersect(t) && !hp.minSphere.TetraContain(t)
}

// ChangeCenter updates the center of the halo to a new position. This includes
// updating several pieces of internal state.
func (hp *HaloProfiles) ChangeCenter(v *geom.Vec) {
	for i := 0; i < 3; i++ {
		hp.C[i] = v[i]
		hp.minSphere.C[i] = v[i]
		for r := range hp.rs {
			hp.rs[r].dr[i] = -v[i]
		}
	}
}

// Mass returns the mass of the halo as estimated by averaging the halo's
// profiles.
func (hp *HaloProfiles) Mass(rhoM float64) float64 {
	vol := (4 * math.Pi / 3) * hp.R*hp.R*hp.R
	return float64(vol) * hp.Rho() * rhoM
}

// rho returns the total enclosed density of the halo as estimate by averaging
// the halo's profiles.
func (hp *HaloProfiles) Rho() float64 {
	rBuf, rSum := make([]float64, hp.bins), make([]float64, hp.bins)
	profs := hp.n * len(hp.rs)

	// Find the spherically averaged rho profile
	for r := 0; r < len(hp.rs); r++ {
		for i := 0; i < hp.n; i++ {
			hp.rs[r].Retrieve(i, rBuf)
			for j := range rBuf { rSum[j] += rBuf[j] }
		}
	}
	for j := range rSum { rSum[j] /= float64(profs) }

	// Integrate
	dr := float64(hp.R - hp.minSphere.R) / float64(len(rSum))
	sum := 0.0
	for i, rho := range rSum {
		r := (float64(i) + 0.5)*dr + float64(hp.minSphere.R)
		sum += r*r*dr*rho
	}

	sum *= 4*math.Pi	
	vol := (4 * math.Pi / 3) * hp.R*hp.R*hp.R
	return sum / float64(vol)
}
