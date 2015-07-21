package geom

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
	return t.Distance(ap, bEnter), t.Distance(ap, bLeave), true
}
