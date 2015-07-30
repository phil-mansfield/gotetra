package geom

// LineEps represents the maximum distance between two two floating point
// numbers at which they are still considered as equal by this module when
// calulating properties of lines. Specifically, when detecting whether or not
// a line is vertical.
//
// You can change it if you want, but since I use float32s, making it much
// smaller than this will probably be a mistake on your part.
var LineEps float32 = 1e-5

// Line is a possibly vertical 2D line.
type Line struct {
	Y0, M float32
	Vertical bool
}

func lineEpsEq(x, y float32) bool {
	return (x + LineEps > y) && (x - LineEps < y)
}

// Init initializes a line so that it passes though both the supplied points.
// Init returns false if a line cannot be unambiguously drawn between the given
// points, false is returned. Otherwise true is returned.
func (l *Line) Init(x1, y1, x2, y2 float32) (ok bool) {
	if lineEpsEq(x1, x2) {
		if lineEpsEq(y1, y2) {
			return false
		}
		l.Y0 = x1
		l.Vertical = true
	} else {
		l.M = (y1 - y2) / (x1 - x2)
		l.Y0 = y1 - l.M * x1
		l.Vertical = false
	}
	return true
}

// InitFromPlucker initializes a line so that it is aligned with the given
// Plucker vector.
func (l *Line) InitFromPlucker(ap *AnchoredPluckerVec) {
	l.Init(ap.P[0], ap.P[1], ap.U[0], ap.U[1])
}

// AreParallel returns true if both the given lines are parallel.
func AreParallel(l1, l2 *Line) bool {
	return (l1.Vertical && l2.Vertical) || lineEpsEq(l1.M, l2.M)
}

// Solve solves for the intersection point between l1 and l2 if it exists. If
// no intersection point exists, ok is returned as false. Otherwise it is
// returned as true.
func Solve(l1, l2 *Line) (x, y float32, ok bool) {
	if AreParallel(l1, l2) { return 0, 0, false }

	if l1.Vertical {
		return l1.Y0, l2.Y0 + l2.M * l1.Y0, true
	} else if l2.Vertical {
		return l2.Y0, l1.Y0 + l1.M * l2.Y0, true
	}

	x = (l2.Y0 - l1.Y0) / (l1.M - l2.M)
	y = l1.Y0 + l1.M * x
	return x, y, true
}
