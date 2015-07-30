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

	// Workspace objects.
	dr geom.Vec
	t geom.Tetra
	pt geom.PluckerTetra
	poly geom.TetraSlice
}

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

	hr.dr = *origin
	for i := 0; i < 3; i++ { hr.dr[i] *= -1 }
}

func (hr *haloRing) insert(t *geom.Tetra, rho float64) {
	hr.t = *t
	hr.t.Translate(&hr.dr)
	hr.t.Rotate(&hr.rot)
	hr.pt.Init(&hr.t) // This is slower than it has to be by a lot!
	
	if hr.t.ZPlaneSlice(&hr.pt, 0, &hr.poly) {
		lowPhi, phiWidth := hr.poly.AngleRange()
		lowIdx, idxWidth := geom.AngleBinRange(lowPhi, phiWidth, hr.n)
		
		idx := lowIdx
		var rEnter, rExit float64
		for i := 0; i < idxWidth; i++ {
			if idx >= hr.n { idx -= hr.n }
			l := &hr.Lines[idx]
			l1, l2 := hr.poly.IntersectingLines(hr.phis[idx])


			if l1 == nil { continue }

			if l2 != nil {
				enterX, enterY, _ := geom.Solve(l, l1)
				exitX, exitY, _ := geom.Solve(l, l2)
				
				rEnterSqr := enterX*enterX + enterY*enterY
				rExitSqr := exitX*exitX + exitY*exitY
				
				rEnter = math.Sqrt(float64(rEnterSqr))
				rExit = math.Sqrt(float64(rExitSqr))
				if rExit < rEnter { rEnter, rExit = rExit, rEnter }
			} else {
				exitX, exitY, _ := geom.Solve(l, l1)
				rExitSqr := exitX*exitX + exitY*exitY
				rExit = math.Sqrt(float64(rExitSqr))
				rEnter = 0.0
			}

			hr.Insert(rEnter, rExit, rho, idx)

			idx++
		}
	}
}

func (hr *haloRing) Count(t *geom.Tetra) { hr.insert(t, 1) }
func (hr *haloRing) Density(t *geom.Tetra, rho float64) { hr.insert(t, rho) }

type HaloProfiles struct {
	geom.Sphere

	rs []haloRing
	rMin, rMax float64
	id, bins, n int
}

func (hp *HaloProfiles) Init(
	id, rings int, origin *geom.Vec, rMin, rMax float64,
	bins, n int,
) *HaloProfiles {
	// We might be able to do better than this.
	solid, ok := geom.NewUniquePlatonicSolid(rings)
	if !ok {
		panic(fmt.Sprintf("Cannot uniformly space %d rings.", rings))
	}

	hp.rs = make([]haloRing, rings)
	hp.C, hp.R = *origin, float32(rMax)
	hp.rMin, hp.rMax = rMin, rMax
	hp.id, hp.bins, hp.n = id, bins, n

	norms := solid.UniqueNormals()
	println(rings, len(norms))
	for i := 0; i < rings; i++ {
		hp.rs[i].Init(&norms[i], origin, rMin, rMax, bins, n)
	}

	return hp
}

func (hp *HaloProfiles) Count(t *geom.Tetra) {
	for i := range hp.rs { hp.rs[i].Count(t) }
}

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

func (hp *HaloProfiles) SheetIntersect(hd *io.SheetHeader) bool {
	tw := float32(hd.TotalWidth)
	return inRange(hp.C[0], hp.R, hd.Origin[0], hd.Width[0], tw) &&
		inRange(hp.C[1], hp.R, hd.Origin[1], hd.Width[1], tw) &&
		inRange(hp.C[2], hp.R, hd.Origin[2], hd.Width[2], tw)
}

func (hp *HaloProfiles) ChangeCenter(v *geom.Vec) {
	for i := 0; i < 3; i++ {
		hp.C[i] = v[i]
		for r := range hp.rs {
			hp.rs[r].dr[i] = -v[i]
		}
	}
}

func (hp *HaloProfiles) Mass(rhoM float64) float64 {
	vol := (4 * math.Pi / 3) * hp.R*hp.R*hp.R
	return float64(vol) * hp.Rho() * rhoM
}

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
	dr := (hp.rMax - hp.rMin) / float64(len(rSum))
	sum := 0.0
	rs := make([]float64, len(rSum))
	for i, rho := range rSum {
		r := (float64(i) + 0.5)*dr + hp.rMin
		rs[i] = r
		sum += r*r*dr*rho
	}

	sum *= 4*math.Pi	
	vol := (4 * math.Pi / 3) * hp.R*hp.R*hp.R
	return sum / float64(vol)
}
