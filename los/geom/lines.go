package geom

var LineEps float32 = 1e-5

type Line struct {
	Y0, M float32
	Vertical bool
}

func lineEpsEq(x, y float32) bool {
	return (x + LineEps > y) && (x - LineEps < y)
}

func (l *Line) Init(x1, y1, x2, y2 float32) {
	if lineEpsEq(x1, x2) {
		if lineEpsEq(y1, y2) {
			panic("Cannot make line between a point and itself.")
		}
		l.Y0 = x1
		l.Vertical = true
	} else {
		l.M = (y1 - y2) / (x1 - x2)
		l.Y0 = y1 - l.M * x1
	}
}

func (l *Line) InitFromPlucker(ap *AnchoredPluckerVec) {
	l.Init(ap.P[0], ap.P[1], ap.U[0], ap.U[1])
}

func AreParallel(l1, l2 *Line) bool {
	return (l1.Vertical && l2.Vertical) || lineEpsEq(l1.M, l2.M)
}

func Solve(l1, l2 *Line) (x, y float32) {
	if l1.Vertical {
		return l1.Y0, l2.Y0 + l2.M * l1.Y0
	} else if l2.Vertical {
		return l2.Y0, l1.Y0 + l1.M * l2.Y0
	}

	x = (l2.Y0 - l1.Y0) / (l1.M - l1.M)
	y = l1.Y0 + l1.M * x
	return x, y		
}
