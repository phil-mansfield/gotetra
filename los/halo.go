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
	geom.EulerMatrixBetweenAt(norm, zAxis, &hr.rot)
	geom.EulerMatrixBetweenAt(zAxis, norm, &hr.irot)
	hr.phis = make([]float32, n)
	for i := 0; i < n; i++ {
		hr.phis[i] = float32(i) / float32(n) * (2 * math.Pi)
	}

	hr.dr = *origin
	for i := 0; i < 3; i++ { hr.dr[i] *= -1 }
}

func (hr *haloRing) Count(t *geom.Tetra) {
	hr.t = *t
	hr.t.Translate(&hr.dr)
	hr.t.Rotate(&hr.rot)
	hr.pt.Init(&hr.t) // This is slower than it has to be by a lot!
	
	if t.ZPlaneSlice(&hr.pt, 0, &hr.poly) {
		lowPhi, phiWidth := hr.poly.AngleRange()
		lowIdx, idxWidth := geom.AngleBinRange(lowPhi, phiWidth, hr.n)
		
		idx := lowIdx
		var rEnter, rExit float64
		for i := 0; i < idxWidth; i++ {
			l := &hr.Lines[idx]
			l1, l2 := hr.poly.IntersectingLines(hr.phis[idx])
			if l2 != nil {
				enterX, enterY, _ := geom.Solve(l, l1)
				exitX, exitY, _ := geom.Solve(l, l2)

				rEnterSqr := enterX*enterX + enterY*enterY
				rExitSqr := exitX*exitX + exitY*exitY
				
				rEnter = math.Sqrt(float64(rEnterSqr))
				rExit = math.Sqrt(float64(rExitSqr))
			} else {
				exitX, exitY, _ := geom.Solve(l, l1)
				rExitSqr := exitX*exitX + exitY*exitY
				rExit = math.Sqrt(float64(rExitSqr))
				rEnter = 0.0
			}

			hr.Insert(rEnter, rExit, 1, idx)

			idx++
			if idx == hr.n { idx = 0 }
		}
	}
}

type HaloProfiles struct {
	rs []haloRing
	sph geom.Sphere
	rMin, rMax float64
	id int
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
	hp.sph.X, hp.sph.Y = origin[0], origin[1]
	hp.sph.Z, hp.sph.R = origin[2], float32(rMax)
	hp.rMin, hp.rMax = rMin, rMax
	hp.id = id

	norms := solid.UniqueNormals()
	for i := 0; i < rings; i++ {
		hp.rs[i].Init(&norms[i], origin, rMin, rMax, bins, n)
	}

	return hp
}

func (hp *HaloProfiles) Count(t *geom.Tetra) {
	for i := range hp.rs {
		hp.rs[i].Count(t)
	}
}

func (hp *HaloProfiles) SphereIntersect(sph *geom.Sphere) bool {
	return hp.sph.Intersect(sph)
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
	s := &hp.sph
	x, y, z := s.X, s.Y, s.Z
	tw := float32(hd.TotalWidth)
	return inRange(x, s.R, hd.Origin[0], hd.Width[0], tw) &&
		inRange(y, s.R, hd.Origin[1], hd.Width[1], tw) &&
		inRange(z, s.R, hd.Origin[2], hd.Width[2], tw)
}
