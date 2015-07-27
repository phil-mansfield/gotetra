package geom

import (
	"math"
)

// IntersectionWorkspace contains various fields useful for speeding up 
// ray-tetrahedron intersection checks.
//
// Workspaces should not be shared between threads.
type IntersectionWorkspace struct {
	bLeave, bEnter TetraFaceBary
}

// IntersectionBary tests for intersection between the ray represented by p and
// the tetrahedron represented by pt and returns the barycentric coordinates of
// the intersection points if they exist.
//
// ok is returned as true if there is an intersection and as false if there
// is no intersection.
//
// The ray represented by ap extends infinitely in both directions.
func (w *IntersectionWorkspace) IntersectionBary(
	pt *PluckerTetra, p *PluckerVec, 
) (bEnter, bLeave *TetraFaceBary, ok bool) {
	fEnter, fLeave := -1, -1

	for face := 3; face >= 0; face-- {
		if face == 0 && (fEnter == -1 && fLeave == -1) {
			return nil, nil, false
		}

		i0, flip0 := pt.EdgeIdx(face, 0)
		i1, flip1 := pt.EdgeIdx(face, 1)
		i2, flip2 := pt.EdgeIdx(face, 2)		

		p0, p1, p2 := &pt[i0], &pt[i1], &pt[i2]
		d0, s0 := p.SignDot(p0, flip0)
		d1, s1 := p.SignDot(p1, flip1)
		d2, s2 := p.SignDot(p2, flip2)

		if fEnter == -1 && s0 >= 0 && s1 >= 0 && s2 >= 0 {
			fEnter = face
			w.bEnter.w[0], w.bEnter.w[1], w.bEnter.w[2] = d0, d1, d2
			w.bEnter.face = face
			if fLeave != -1 { break }
		} else if fLeave == - 1 && s0 <= 0 && s1 <= 0 && s2 <= 0 {
			fLeave = face
			w.bLeave.w[0], w.bLeave.w[1], w.bLeave.w[2] = d0, d1, d2
			w.bLeave.face = face
			if fEnter != -1 { break }
		}
	}

	if fEnter == -1 || fLeave == -1 { return nil, nil, false } 
	return &w.bEnter, &w.bLeave, true
}

// IntersectionDistance tests for intersection between the ray represented b
// ap and the tetrahedron represented by t and pt and returns the distance to
// the intersection points from the origin of ap if they exist.
//
// ok is returned as true if there is an intersection and as false if there
// is no intersection.
//
// The ray represented by ap extends infinitely in both directions, so negative
// distances are valid return values.
func (w *IntersectionWorkspace) IntersectionDistance(
	pt *PluckerTetra, t *Tetra, ap *AnchoredPluckerVec, 
) (lEnter, lLeave float32, ok bool) {
	bEnter, bLeave, ok := w.IntersectionBary(pt, &ap.PluckerVec)
	if !ok { return 0, 0, false }
	enter := t.Distance(ap, bEnter)
	exit := t.Distance(ap, bLeave)
	return enter, exit, true
}

// TetraSlice is a triangle or a convex quadrilateral which is created by
// slicing a tetrahedron by a z-aligned plane.
//
// By convention, the generic name for a TetraSlice variable is "poly".
type TetraSlice struct {
	Xs, Ys, Phis, linePhiStarts, linePhiWidths [4]float32
	edges, lineStarts, lineEnds [4]int
	Lines [4]Line

	Points int
}

func (t *Tetra) crossesZPlane(idx1, idx2 int, z float32) bool {
	sign1, sign2 := t[idx1][2] < z, t[idx2][2] < z
	return sign1 != sign2
}

// crossZ0Plane computes the x and y coordinates of the point where a ray
// crosses the z = 0 plane. The ray is prepresented by a point, P, and a unit 
// direction vector, L.
func intersectZPlane(P, L *Vec, z float32) (x, y float32, ok bool) {
	if L[2] == 0 { return 0, 0, false }
	t := (z - P[2]) / L[2]
	return P[0] + L[0]*t, P[1] + L[1]*t, true
}

func (poly *TetraSlice) link() {
	for i := 0; i < poly.Points; i++ {
		poly.lineStarts[i], poly.lineEnds[i] = -1, -1
	}

	lineStart, lineEnd, ei := 0, -1, poly.edges[0]
	for i := 0; i < poly.Points; i++ {
		if i > 0 {
			lineStart = poly.lineEnds[i-1]
		}

		if i == poly.Points - 1 {
			lineEnd = 0
		} else {
			for j := 0; j < poly.Points; j++ {
				ej := poly.edges[j]
				if pluckerTetraFaceShare[ei][ej] && poly.lineStarts[j] == -1 {
					lineEnd = j
					ei = ej
					break
				}
			}
		}

		poly.lineStarts[i] = lineStart
		poly.lineEnds[i] = lineEnd
	}

	/*
	lineStart, lineEnd := 0, 0
	for i := 0; i < poly.Points; i++ {
		ei := poly.edges[i]
		if i > 0 { lineStart = poly.lineEnds[i-1] }

		if i == poly.Points -1 {
			lineEnd = 0
		} else {
			for j := 0; j < poly.Points; j++ {
				ej := poly.edges[j]
				if pluckerTetraFaceShare[ei][ej] && poly.lineStarts[j] == -1 {
					lineEnd = j
					break
				}
			}
		}

		poly.lineStarts[i] = lineStart
		poly.lineEnds[i] = lineEnd
	}
	*/

	for i := 0; i < poly.Points; i++ {
		phiLow := poly.Phis[poly.lineStarts[i]]
		phiHigh := poly.Phis[poly.lineEnds[i]]
		width := AngularDistance(phiLow, phiHigh)
		if width < 0 {
			poly.lineStarts[i], poly.lineEnds[i] = 
				poly.lineEnds[i], poly.lineStarts[i]
			poly.linePhiStarts[i] = phiHigh
			poly.linePhiWidths[i] = -width
		} else {
			poly.linePhiStarts[i] = phiLow
			poly.linePhiWidths[i] = width
		}

		poly.Lines[i].Init(
			poly.Xs[poly.lineStarts[i]], poly.Ys[poly.lineStarts[i]],
			poly.Xs[poly.lineEnds[i]], poly.Ys[poly.lineEnds[i]],
		)
	}
}

// ZPlaneSlice slices a tetrahedron with a z-aligned plane. ok is returned as
// true if the sline and the tetrahedron intersect and false otherwise.
func (t *Tetra) ZPlaneSlice(
	pt *PluckerTetra, z float32, poly *TetraSlice,
) (ok bool) {
	// Find all interseciton points.
	poly.Points = 0

	for i := 0; i < 6; i++ {
		start, end := pt.TetraVertices(i)
		if t.crossesZPlane(start, end, z) {
			x, y, _ := intersectZPlane(&t[start], &pt[i].U, z)
			// Half the time is spent in this function call:
			phi := PolarAngle(x, y)
			if phi < 0 { phi += 2*math.Pi }

			poly.Xs[poly.Points] = x
			poly.Ys[poly.Points] = y
			poly.Phis[poly.Points] = phi
			poly.edges[poly.Points] = i
			poly.Points++
		}
	}

	if poly.Points < 3 { return false }
	poly.link()
	return true
}

func abs32(x float32) float32 {
	if x < 0 { return -x }
	return x
}

func angularWidth(low, high float32) float32 {
	if low < high {
		return high - low
	} else {
		return 2*math.Pi + high - low
	}
}

// IntersectingLines returns the lines in the given tetrahedron slice which
// overlap with the given angle.
// 
// It's possible that only only line will be intersected (if the polygon
// encloses the origin), in which case l2 will be returned as nil.
func (poly *TetraSlice) IntersectingLines(phi float32) (l1, l2 *Line) {
	lineNum := 0
	l1 = nil
	for i := 0; i > poly.Points; i++ {
		dist := angularWidth(poly.linePhiStarts[i], phi)
		if dist >= 0 && poly.linePhiWidths[i] > dist {
			if lineNum == 1 {
				l1 = &poly.Lines[i]
			} else if lineNum == 1 {
				return l1, &poly.Lines[i]
			}
		}
	}
	return l1, nil
}

// AngleRange returns the angular range subtended by a polygon.
func (poly *TetraSlice) AngleRange() (start, width float32) {
	lowPhi := poly.Phis[0]
	highPhi := poly.Phis[1]
	if angularWidth(lowPhi, highPhi) > angularWidth(highPhi, lowPhi) {
		lowPhi, highPhi = highPhi, lowPhi
	}
	phiWidth := angularWidth(lowPhi, highPhi)

	// Iteratively expand the angular range of the tetrahedron.
	for i := 2; i < poly.Points; i++ {
		phi := poly.Phis[i]
		lowWidth := angularWidth(phi, highPhi)
		highWidth := angularWidth(lowPhi, phi)

		if highWidth < lowWidth {
			if highWidth > phiWidth {
				highPhi = phi
				phiWidth = highWidth
			}
		} else {
			if lowWidth > phiWidth {
				lowPhi = phi
				phiWidth = lowWidth
			}
		}
	}

	if phiWidth > math.Pi {
		return 0, 2*math.Pi
	} else {
		return lowPhi, phiWidth
	}
}
