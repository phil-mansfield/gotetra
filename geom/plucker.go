package geom

type PluckerVec struct {
	Us, Vs Vec
}

type pluckerTri [3]PluckerVec

type PluckerTetra [4]pluckerTri

func NewPluckerVec(ls, xs *Vec) *PluckerVec {
	p := &PluckerVec{}
	p.Init(ls, xs)
	return p
}

func (p *PluckerVec) Init(ls, xs *Vec) {
	// Don't reverse this order since p.Us is used as a buffer
	ls.CrossAt(xs, &p.Vs)
	p.Us = *ls
}

// Pdot computes the permuted inner product of two Plucker vectors.
func (p1 *PluckerVec) PDot(p2 *PluckerVec) float64 {
	return p1.Us.Dot(&p2.Vs) + p2.Us.Dot(&p1.Vs)
}

func NewPluckerTetra(r1, r2, r3, r4 *Vec) *PluckerTetra {
	tet := &PluckerTetra{}
	tet.Init(r1, r2, r3, r4)
	return tet
}

func (tet *PluckerTetra) Init(r1, r2, r3, r4 *Vec) {
	
}

func (tet *PluckerTetra) initLine(r1, r2, l12 *Vec) {
}
