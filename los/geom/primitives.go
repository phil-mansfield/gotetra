/*package geom contains routines for computing geometric quantities.

Contains implementations of algorithms described in Platis & Theoharis, 2015
as well as Schneider & Eberly.

The calling signatures might be more convoluted than they have to be because
I was too worried about memory consumption when I started it. I should include
examples.
*/
package geom

import (
	"math"
)

// Vec is a three dimensional vector. (Duh!)
type Vec [3]float32

// PluckerVec represents a ray. If P and L are the position of the ray's 
// origin and the unit vector representing its direction, respectively, then
// U = L and V = L cross P.
type PluckerVec struct {
	U, V Vec
}

// AnchoredPluckerVec is a Plucker vector which also saves the position of
// the ray's origin.
type AnchoredPluckerVec struct {
	PluckerVec
	P Vec
}

// Init initializes a Plucker vector given a ray origin, P, and a unit
// direction vector, L.
func (p *PluckerVec) Init(P, L *Vec) {
	p.U = *L
	
	p.V[0] = -P[1]*L[2] + P[2]*L[1]
    p.V[1] = -P[2]*L[0] + P[0]*L[2]
    p.V[2] = -P[0]*L[1] + P[1]*L[0]
}

// InitFromSegment initialized a Plucker vector which corresponds to a ray
// pointing from the position vector P1 to the position vector P2.
func (p *PluckerVec) InitFromSegment(P1, P2 *Vec) {
	var sum float32
	for i := 0; i < 3; i++ {
		p.U[i] = P2[i] - P1[i]
		sum += p.U[i]*p.U[i]
	}
	sum = float32(math.Sqrt(float64(sum)))
	//for i := 0; i < 3; i++ { p.U[i] /= sum }

	p.V[0] = -P1[1]*p.U[2] + P1[2]*p.U[1]
    p.V[1] = -P1[2]*p.U[0] + P1[0]*p.U[2]
    p.V[2] = -P1[0]*p.U[1] + P1[1]*p.U[0]
}

// Translate translates a Plucker vector along the given vector.
func (p *PluckerVec) Translate(dx *Vec) {
	p.V[0] += -dx[1]*p.U[2] + dx[2]*p.U[1]
    p.V[1] += -dx[2]*p.U[0] + dx[0]*p.U[2]
    p.V[2] += -dx[0]*p.U[1] + dx[1]*p.U[0]
}

// Dot computes the permuted inner product of p1 and p2, i.e.
// p1.U*p2.V + p1.V*p2.U.
func (p1 *PluckerVec) Dot(p2 *PluckerVec, flip bool) float32 {
	var sum float32
	for i := 0; i < 3; i++ {
		sum += p1.U[i]*p2.V[i] + p1.V[i]*p2.U[i]
	}
	if flip {
		return sum
	} else {
		return -sum
	}
}

// Dot computes the permuted inner product of p1 and p2, i.e.
// p1.U*p2.V + p1.V*p2.U and also returns a sign flag of -1, 0, or +1 if
// that product is negative, zero, or positive, respectively.
func (p1 *PluckerVec) SignDot(p2 *PluckerVec, flip bool) (float32, int) {
	dot := p1.Dot(p2, flip)
	if dot == 0 {
		return dot, 0
	} else if dot > 0 {
		return dot, +1
	} else {
		return dot, -1
	}
}

// Init initializes an anchored Plucker vector given a ray origin, P, and a
// unit direction vector, L.
func (ap *AnchoredPluckerVec) Init(P, L *Vec) {
	ap.PluckerVec.Init(P, L)
	ap.P = *P
}

// InitFromSegment initialized a Plucker vector which corresponds to a ray
// pointing from the position vector P1 to the position vector P2.
func (ap *AnchoredPluckerVec) InitFromSegment(P1, P2 *Vec) {
	ap.PluckerVec.InitFromSegment(P1, P2)
	ap.P = *P1
}

// Translate translates a Plucker vector along the given vector.
func (ap *AnchoredPluckerVec) Translate(dx *Vec) {
	ap.PluckerVec.Translate(dx)
	for i := 0; i < 3; i++ { ap.P[i] += dx[i] }
}

// Tetra is a tetrahedron. (Duh!)
//
// Face ordering is:
// F0(V3, V2, V1)
// F1(V2, V3, V0)
// F2(V1, V0, V3)
// F3(V0, V1, V2)
type Tetra [4]Vec

var tetraIdxs = [4][3]int {
	[3]int{ 3, 2, 1 },
	[3]int{ 2, 3, 0 },
	[3]int{ 1, 0, 3 },
	[3]int{ 0, 1, 2 },
}

// VertexIdx returns the index into the given tetrahedron corresponding to
// the specified face and vertex.
func (_ *Tetra) VertexIdx(face, vertex int) int {
	return tetraIdxs[face][vertex]
}

// Orient arranges tetrahedron points so that all faces point outward for
// dir = +1 and inward for dir = -1.
func (t *Tetra) Orient(dir int) {
	v, w, n := Vec{}, Vec{}, Vec{}
	for i := 0; i < 3; i++ {
		v[i] = t[1][i] - t[0][i]
		w[i] = t[2][i] - t[0][i]
	}
	n[0] = -v[1]*w[2] + v[2]*w[1]
    n[1] = -v[2]*w[0] + v[0]*w[2]
    n[2] = -v[0]*w[1] + v[1]*w[0]

	var dot float32
	for i := 0; i < 3; i++ {
		dot += n[i] * (t[3][i] - t[0][i])
	}

	if (dot < 0 && dir == -1) || (dot > 0 && dir == +1) {
		t[0], t[1] = t[1], t[0]
	}
}

// Translate translates a tetrahedron by the given vector.
func (t *Tetra) Translate(dx *Vec) {
	for i := 0; i < 4; i++ {
		for j := 0; j < 3; j++ {
			t[i][j] += dx[j]
		}
	}
}

// Volume computes the volume of a tetrahedron.
func (t *Tetra) Volume() float64 {
	leg10, leg11, leg12 := t[1][0]-t[0][0], t[1][1]-t[0][1], t[1][2]-t[0][2]
	leg20, leg21, leg22 := t[2][0]-t[0][0], t[2][1]-t[0][1], t[2][2]-t[0][2]
	leg30, leg31, leg32 := t[3][0]-t[0][0], t[3][1]-t[0][1], t[3][2]-t[0][2]

	n0 := leg10*(leg21*leg32 - leg22*leg31)
    n1 := leg11*(leg22*leg30 - leg20*leg32)
    n2 := leg12*(leg20*leg31 - leg21*leg30)
	vol := float64((n0 + n1 + n2) / 6)
	if vol < 0 { return -vol }
	return vol
}

// Sphere is exactly what you think it is.
type Sphere struct {
	C Vec
	R float32
}

// Intersect returns true if the two spheres intersect and false otherwise.
func (s1 *Sphere) SphereIntersect(s2 *Sphere) bool {
	sum := float32(0)
	dr := s1.R + s2.R
	dr2 := dr*dr
	for i := 0; i < 3; i++ {
		dx := s1.C[0] - s2.C[0]
		sum += dx*dx
		if sum > dr2 { return false }
	}
	return true
}

func (s *Sphere) VecIntersect(v *Vec) bool {
	sum := float32(0)
	dr2 := s.R*s.R
	for i := 0; i < 3; i++ {
		dx := s.C[0] - v[0]
		sum += dx*dx
		if sum > dr2 { return false }
	}
	return true

}

// BoundingSphere draws a bounding sphere aorund the given tetrahedron.
func (t *Tetra) BoundingSphere(sph *Sphere) {
	bx := (t[0][0] + t[1][0] + t[2][0] + t[3][0]) / 4
	by := (t[0][1] + t[1][1] + t[2][1] + t[3][1]) / 4
	bz := (t[0][2] + t[1][2] + t[2][2] + t[3][2]) / 4

	dx, dy, dz := bx-t[0][0], by-t[0][1], bz-t[0][2]
	maxRSqr := dx*dx + dy*dy + dz*dz
	for i := 1; i < 4; i++ {
		dx, dy, dz = bx-t[i][0], by-t[i][1], bz-t[i][2]
		rSqr := dx*dx + dy*dy + dz*dz
		if rSqr > maxRSqr { maxRSqr = rSqr }
	}

	sph.C[0], sph.C[1], sph.C[2] = bx, by, bz
	sph.R = float32(math.Sqrt(float64(maxRSqr)))
}

// TetraFaceBary contains information specifying the barycentric coordinates
// of a point on a face of a tetrahedron.
type TetraFaceBary struct {
	w [3]float32
	face int
}

// Distance calculates the distance from an anchored Plucker vector to a point
// in a tetrahedron described by the given unscaled barycentric coordinates.
func (t *Tetra) Distance(ap *AnchoredPluckerVec, bary *TetraFaceBary) float32 {
	// Computes one coordinate of the intersection point from the barycentric
	// coordinates of the intersection, then solves P_intr = P + t * L for t.
	var sum float32
	for i := 0; i < 3; i++ { sum += bary.w[i] }
	u0, u1, u2 := bary.w[0] / sum, bary.w[1] / sum, bary.w[2] / sum
	var dim int
	for dim = 0; dim < 3; dim++ {
		if ap.U[dim] > 1e-6 || ap.U[dim] < -1e-6 { break }
	}

	p0 := t[t.VertexIdx(bary.face, 0)][dim]
	p1 := t[t.VertexIdx(bary.face, 1)][dim]
	p2 := t[t.VertexIdx(bary.face, 2)][dim]

	d := ((u0*p0 + u1*p1 + u2*p2) - ap.P[dim]) / ap.U[dim]

	return d
}

// PluckerTetra is a tetrahedron represented by the Plucker vectors that make
// up its edges. It is used for Platis & Theoharis's interseciton detection
// algorithm.
//
// The raw ordering of edges is
// F0(V3, V2, V1)
// F1(V2, V3, V0)
// F2(V1, V0, V3)
// F3(V0, V1, V2)
// {0-1, 0-2, 0-3, 1-2, 1-3, 2-3}
type PluckerTetra [6]PluckerVec

var pluckerTetraEdges = [4][3]int{
	[3]int{ 3, 4, 5 }, // 2-1, 1-3, 3-2
	[3]int{ 2, 1, 5 }, // 3-0, 0-2, 2-3
	[3]int{ 2, 4, 0 }, // 0-3, 3-1, 1-0
	[3]int{ 3, 1, 0 }, // 1-2, 2-0, 0-1
}

var pluckerTetraFlips = [4][3]bool{
	[3]bool{true, false, true},
	[3]bool{true, false, false},
	[3]bool{false, true, true},
	[3]bool{false, true, false},
}

var pluckerTetraFaceShare = [6][6]bool {
	[6]bool{ false, true,  true,  true,  true,  false },
	[6]bool{ true,  false, true,  true,  false, true  },
	[6]bool{ true,  true,  false, false, true,  true  },
	[6]bool{ true,  true,  false, false, true,  true  },
	[6]bool{ true,  false, true,  true,  false, true  },
	[6]bool{ false, true,  true,  true,  true,  false },
}

var tetraEdgeStarts = [6]int{ 0, 0, 0, 1, 1, 2 }
var tetraEdgeEnds = [6]int{ 1, 2, 3, 2, 3, 3 }

// Init initializes a Plucker Tetrahedron from a normal Tetrahedron.
func (pt *PluckerTetra) Init(t *Tetra) {
	pt[0].InitFromSegment(&t[0], &t[1])
	pt[1].InitFromSegment(&t[0], &t[2])
	pt[2].InitFromSegment(&t[0], &t[3])
	pt[3].InitFromSegment(&t[1], &t[2])
	pt[4].InitFromSegment(&t[1], &t[3])
	pt[5].InitFromSegment(&t[2], &t[3])
}

// Translate translates a Plucker tetrahedron along the given vector.
func (pt *PluckerTetra) Translate(dx *Vec) {
	for i := 0; i < 6; i++ { pt[i].Translate(dx) }
}

// EdgeIdx returns the index into pt which corresponds to the requested
// face and edge. A flag is also returned indicating whether the vector stored
// in pt needs to be flipped when doing operations on that face.
func (_ *PluckerTetra) EdgeIdx(face, edge int) (idx int, flip bool) {
	idx = pluckerTetraEdges[face][edge]
	flip = pluckerTetraFlips[face][edge]
	return idx, flip
}

// TetraVertexIdx returns the indices of the vertices in a Tetra object which
// correspond to end points of a given PluckerVector within a PluckerTetra.
func (_ *PluckerTetra) TetraVertices(i int) (start, end int) {
	return tetraEdgeStarts[i], tetraEdgeEnds[i]
}
