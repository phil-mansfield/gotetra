package los

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"

	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/math/mat"
	"github.com/phil-mansfield/gotetra/math/sort"
)

type haloRing struct {
	ProfileRing
	phis []float32
	rot, irot mat.Matrix32
	norm geom.Vec
	rMax float64

	// Options
	log bool

	// Workspace objects.
	dr geom.Vec
	t geom.Tetra
	pt geom.PluckerTetra
	poly geom.TetraSlice
}

func (hr *haloRing) Reuse(origin *geom.Vec, rMin, rMax float64) {
	hr.dr = *origin
	for i := 0; i < 3; i++ { hr.dr[i] *= -1 }
	hr.rMax = rMax
	hr.ProfileRing.Reuse(rMin, rMax)
	hr.Clear()
}

type internalOption func(*haloRing)
type Option internalOption

func Log(log bool) Option {
	return func(hr *haloRing) { hr.log = log }
}

func Rotate(phi, theta, psi float32) Option {
	return func(hr *haloRing) {
		rot := geom.EulerMatrix(phi, theta, psi)
		hr.norm.Rotate(rot)
	}
}

// Init initialized a haloRing.
func (hr *haloRing) Init(
	norm, origin *geom.Vec, rMin, rMax float64, bins, n int, opts ...Option,
) {
	hr.log = false
	hr.norm = *norm
	for _, opt := range opts { opt(hr) }

	if hr.log {
		if rMax <= 0 || rMin <= 0 {
			panic("Non-positive bounding radius given to logarithmic haloRing.")
		}

		rMin = math.Log(rMin)
		rMax = math.Log(rMax)
	}

	hr.ProfileRing.Init(rMin, rMax, bins, n)
	zAxis := &geom.Vec{0, 0, 1}
	
	hr.rot.Init(make([]float32, 9), 3, 3)
	hr.irot.Init(make([]float32, 9), 3, 3)
	geom.EulerMatrixBetweenAt(&hr.norm, zAxis, &hr.rot)
	geom.EulerMatrixBetweenAt(zAxis, &hr.norm, &hr.irot)
	hr.phis = make([]float32, n)
	for i := 0; i < n; i++ {
		hr.phis[i] = float32(i) / float32(n) * (2 * math.Pi)
	}

	hr.rMax = rMax
	hr.dr = *origin
	for i := 0; i < 3; i++ { hr.dr[i] *= -1 }
}

// This is it. This is the critical innermost for loop.
func (hr *haloRing) insert(t *geom.Tetra, rho float64) {
	hr.t = *t
	hr.t.Translate(&hr.dr)
	hr.t.Rotate(&hr.rot)
	hr.pt.Init(&hr.t) // This is slower than it has to be by a lot!

	var rSqrMin, rSqrMax float32
	if hr.log {
		rMin, rMax := float32(math.Exp(hr.lowR)), float32(math.Exp(hr.highR))
		rSqrMin, rSqrMax = rMin*rMin, rMax*rMax
	} else {
		rSqrMin, rSqrMax = float32(hr.lowR*hr.lowR), float32(hr.highR*hr.highR)
	}
	if hr.t.ZPlaneSlice(&hr.pt, 0, &hr.poly) {
		// Stop early if possible.
		rSqrTMin, rSqrTMax := hr.poly.RSqrMinMax()
		if rSqrTMin > rSqrMax || rSqrTMax < rSqrMin {
			return
		}

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

				if hr.log {
					rEnter = math.Log(float64(rSqrEnter)) / 2
					rExit = math.Log(float64(rSqrExit)) / 2
				} else {
					rEnter = math.Sqrt(float64(rSqrEnter))
					rExit = math.Sqrt(float64(rSqrExit))
				}
				if rExit < rEnter { rEnter, rExit = rExit, rEnter }
			} else {
				// The polygon encloses the origin.
				exitX, exitY, _ := geom.Solve(l, l1)
				rSqrExit := exitX*exitX + exitY*exitY
				if hr.log {
					rExit = math.Log(float64(rSqrExit))
					rEnter = math.Inf(-1)
				} else {
					rExit = math.Sqrt(float64(rSqrExit))
					rEnter = 0.0
				}
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

// BinIdx computes the radial bin containing the given radius.
func (hr *haloRing) BinIdx(r float64) int {
	if hr.log {
		return int((math.Log(r) - hr.lowR) / hr.ProfileRing.dr)
	} else {
		return int((r - hr.lowR) / hr.ProfileRing.dr)
	}
}

// LineSegment write a line segment corresponding to the given profile 
func (hr *haloRing) LineSegment(prof int, out *geom.LineSegment) {
	vec := geom.Vec{}
	sin, cos := math.Sincos(float64(hr.phis[prof]))
	vec[0], vec[1] = float32(cos), float32(sin)
	vec.Rotate(&hr.irot)
	
	if hr.log {
		*out = geom.LineSegment{Origin: hr.dr, Dir: vec,
			StartR: float32(math.Exp(hr.lowR)),
			EndR: float32(math.Exp(hr.highR)) }
		for i := 0; i < 3; i++ { out.Origin[i] *= -1 }
	} else {
		*out = geom.LineSegment{Origin: hr.dr, Dir: vec,
			StartR: float32(hr.lowR), EndR: float32(hr.highR) }
		for i := 0; i < 3; i++ { out.Origin[i] *= -1 }
	}
}

// Add adds the contents of hp2 to hp1.
func (hp1 *HaloProfiles) Add(hp2 *HaloProfiles) {
	for i := range hp1.rs { hp1.rs[i].Add(&hp2.rs[i]) }
}

// Clear resets the conents of the HaloProfiles.
func (hp *HaloProfiles) Clear() {
	for i := range hp.rs { hp.rs[i].Clear() }
}

func (hp *HaloProfiles) Reuse(id int, origin *geom.Vec, rMin, rMax float64) {
	hp.C, hp.R = *origin, float32(rMax)
	hp.cCopy = hp.C
	hp.minSphere.C, hp.minSphere.R = *origin, float32(rMin)
	hp.rMin, hp.rMax = rMin, rMax
	hp.id = id
	for i := range hp.rs { hp.rs[i].Reuse(origin, rMin, rMax) }
}

func ParallelClearHaloProfiles(hs []HaloProfiles) {	
	workers := len(hs)
	runtime.GOMAXPROCS(workers)

	out := make(chan int, workers)
	for i := 0; i < workers -1; i++ { go hs[i].chanClear(out) }
	hs[workers - 1].chanClear(out)

	for i := 0; i < workers; i++ { <-out }
}

func (hp *HaloProfiles) chanClear(out chan<- int) {
	hp.Clear()
	out <- 1
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
	boxWidth float32

	log bool
	IsValid bool
}

// Init initializes a HaloProfiles struct with the given parameters.
func (hp *HaloProfiles) Init(
	id, rings int, origin *geom.Vec, rMin, rMax float64,
	bins, n int, boxWidth float64, opts ...Option,
) *HaloProfiles {
	// We might be able to do better than this.
	var norms []geom.Vec
	switch {
	case rings > 10:
		solid, _ := geom.NewUniquePlatonicSolid(10)
		norms = solid.UniqueNormals()
		for i := 10; i < rings; i++ {
			v := geom.Vec{
				float32(rand.Float64() - 0.5),
				float32(rand.Float64() - 0.5),
				float32(rand.Float64() - 0.5),
			}
			sum := 0.0
			for _, x := range v { sum += float64(x*x) }
			for i := range v { v[i] /= float32(math.Sqrt(sum)) }
			norms = append(norms, v)

		}

	case rings >= 3:
		solid, ok := geom.NewUniquePlatonicSolid(rings)
		norms = solid.UniqueNormals()
		if !ok {
			panic(fmt.Sprintf("Cannot uniformly space %d rings.", rings))
		}
	case rings == 2:
		norms = []geom.Vec{{0, 0, 1}, {0, 1, 0}}
	case rings == 1:
		norms = []geom.Vec{{0, 0, 1}}
	default:
		panic("Invalid ring number.")
	}

	hp.rs = make([]haloRing, rings)
	hp.C, hp.R = *origin, float32(rMax)
	hp.cCopy = hp.C
	hp.minSphere.C, hp.minSphere.R = *origin, float32(rMin)
	hp.rMin, hp.rMax = rMin, rMax
	hp.id, hp.bins, hp.n = id, bins, n
	hp.boxWidth = float32(boxWidth)

	for i := 0; i < rings; i++ {
		hp.rs[i].Init(&norms[i], origin, rMin, rMax, bins, n, opts...)
	}

	hp.log = hp.rs[0].log
	hp.IsValid = true
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
	rBuf, rAvg := make([]float64, hp.bins), make([]float64, hp.bins)
	count := make([]int, hp.bins)

	// Find the spherically averaged rho profile
	for r := 0; r < len(hp.rs); r++ {
		for i := 0; i < hp.n; i++ {
			hp.GetRhos(r, i, rBuf)
			hp.rs[r].Retrieve(i, rBuf)
			for j := range rBuf {
				if !math.IsNaN(rBuf[j]) {
					rAvg[j] += rBuf[j]
					count[j]++
				}
			}
		}
	}

	for j := range rAvg {
		if count[j] == 0 { 
			rAvg[j] = 0 
		} else {
			rAvg[j] /= float64(count[j])
		}
	}

	// Integrate
	sum := 0.0
	if hp.log {
		minlr := hp.rs[0].lowR
		dlr := hp.rs[0].ProfileRing.dr
		for i, rho := range rAvg {
			lr := (float64(i) + 0.5)*dlr + minlr
			r := math.Exp(lr)
			sum += r*r*r*dlr*rho
		}
	} else {
		dr := float64(hp.R - hp.minSphere.R) / float64(len(rAvg))
		for i, rho := range rAvg {
			r := (float64(i) + 0.5)*dr + float64(hp.minSphere.R)
			sum += r*r*dr*rho
		}
	}

	sum *= 4*math.Pi	
	vol := (4 * math.Pi / 3) * hp.R*hp.R*hp.R
	return sum / float64(vol)
}

func (hp *HaloProfiles) ID() int { return hp.id }
func (hp *HaloProfiles) Rings() int { return len(hp.rs) }
func (hp *HaloProfiles) Bins() int { return hp.bins }
func (hp *HaloProfiles) Profiles() int { return hp.n }

func (hp *HaloProfiles) Phi(prof int) float32 { return hp.rs[0].phis[prof] }

func (hr *HaloProfiles) LineSegment(ring, prof int, out *geom.LineSegment) {
	hr.rs[ring].LineSegment(prof, out)
}

func (hp *HaloProfiles) GetRs(out []float64) {
	if hp.bins != len(out) { panic("Length of out array != hp.Bins().") }

	rMin := hp.rs[0].lowR
	dr := hp.rs[0].ProfileRing.dr
	for i := range out { out[i] = (float64(i) + 0.5) * dr + rMin }
	if hp.log {
		for i := range out { out[i] = math.Exp(out[i]) }
	}
}

func (hp *HaloProfiles) PlaneToVolume(
	ring int, px, py float64,
) (x, y, z float64) {
	v := &geom.Vec{float32(px), float32(py), 0}
	v.Rotate(&hp.rs[ring].irot)
	return float64(v[0]), float64(v[1]), float64(v[2])
}

func (hp *HaloProfiles) VolumeToPlane(
	ring int, x, y, z float64,
) (px, py float64) {
	v := &geom.Vec{float32(x), float32(y), float32(z)}
	v.Rotate(&hp.rs[ring].rot)
	return float64(v[0]), float64(v[1])
}

func dist(x, y float32) float32 {
	diff := x - y
	if diff < 0 { return -diff }
	return diff
}

func change(x, anchor, bw float32) float32 {
	normDist := dist(x, anchor)
	if dist(x + bw, anchor) < normDist { return +bw }
	if dist(x - bw, anchor) < normDist { return -bw }
	return 0
}

func (hp *HaloProfiles) GetRhos(ring, prof int, out []float64) {
	if hp.bins != len(out) { panic("Length of out array != hp.Bins().") }
	r := &hp.rs[ring]
	r.Retrieve(prof, out)
}

func LoadDensities(
	hs []HaloProfiles, hds []io.SheetHeader,
	files []string, buf *Buffers,
) {
	ptrs := make([]*HaloProfiles, len(hs))
	for i := range ptrs { ptrs[i] = &hs[i] }
	LoadPtrDensities(ptrs, hds, files, buf)
}

func LoadPtrDensities(
	hs []*HaloProfiles, hds []io.SheetHeader,
	files []string, buf *Buffers,
) {
	for i, file := range files {
		hd := &hds[i]
		WrapHalo(hs, hd)
		buf.ParallelRead(file, hd)
		for j := range hs {
			buf.ParallelDensity(hs[j])
		}
	}
}

func (h *HaloProfiles) MedianProfile() []float64 {
	// Read Densities
	rhoBufs := make([][]float64, h.n * len(h.rs))
	for i := range rhoBufs {
		rhoBufs[i] = make([]float64, h.bins)
	}
	
	idx := 0
	for ring := range h.rs {
		for prof := 0; prof < h.n; prof++ {
			h.GetRhos(ring, prof, rhoBufs[idx])
			idx++
		}
	}
	
	// Find median of each radial bin.
	medBuf := make([]float64, h.n * len(h.rs))
	out := make([]float64, h.bins)
	for j := 0; j < h.bins; j++ {
		for i := range medBuf {
			medBuf[i] = rhoBufs[i][j]
		}
		out[j] = sort.Median(medBuf, medBuf)
	}
	
	return out
}

func (hp *HaloProfiles) MeanProfile() []float64 {
	mean := make([]float64, hp.bins)
	buf := make([]float64, hp.bins)
	
	// Find the spherically averaged rho profile
	for r := 0; r < len(hp.rs); r++ {
		for i := 0; i < hp.n; i++ {
			hp.GetRhos(r, i, buf)
			for j := range buf { mean[j] += buf[j] }
		}
	}

	n := float64(len(hp.rs) * hp.n)
	for j := range mean { mean[j] /= n }
	return mean
}
