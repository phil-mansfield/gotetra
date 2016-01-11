package analyze

import (
	"math"

	"github.com/phil-mansfield/gotetra/los"
	"github.com/phil-mansfield/gotetra/los/geom"
)

type RingBuffer struct {
	PlaneXs, PlaneYs, Xs, Ys, Zs, Rs, Phis []float64
	Oks []bool
	profRs, profRhos, smoothRhos, smoothDerivs []float64
	N, Bins int
}

func (r *RingBuffer) Init(n, bins int) {
	r.N, r.Bins = n, bins

	r.PlaneXs, r.PlaneYs = make([]float64, n), make([]float64, n)
	r.Xs, r.Ys, r.Zs = make([]float64, n), make([]float64, n), make([]float64,n)
	r.Phis, r.Rs = make([]float64, n), make([]float64, n)
	r.Oks = make([]bool, n)

	r.smoothRhos, r.smoothDerivs = make([]float64, bins), make([]float64, bins)
	r.profRs, r.profRhos = make([]float64, bins), make([]float64, bins)
}

func (r *RingBuffer) Clear() {
	for i := 0; i < r.N; i++ {
		r.PlaneXs[i], r.PlaneYs[i] = 0, 0 
		r.Xs[i], r.Ys[i], r.Zs[i] = 0, 0, 0
		r.Phis[i], r.Rs[i] = 0, 0
		r.Oks[i] = false
	}

	for i := 0; i < r.Bins; i++ {
		r.smoothRhos[i], r.smoothDerivs[i] = 0, 0
		r.profRs[i], r.profRhos[i] = 0, 0
	}
}

func (r *RingBuffer) Splashback(
	h los.Halo, ring int, window int, dLim float64,
) {
	h.GetRs(r.profRs)
	ls := new(geom.LineSegment)
	for i := 0; i < r.N; i++ {
		h.GetRhos(ring, i, r.profRhos)
		_, _, r.Oks[i] = Smooth(
			r.profRs, r.profRhos, window,
			Vals(r.smoothRhos), Derivs(r.smoothDerivs),
		)
		if !r.Oks[i] { continue }
		r.Rs[i], r.Oks[i] = SplashbackRadius(
			r.profRs, r.smoothRhos, r.smoothDerivs, DLim(dLim),
		)
		if !r.Oks[i] { continue }

		r.Phis[i] = float64(h.Phi(i))
		if r.Phis[i] < 0 { r.Phis[i] += math.Pi }
		sin, cos := math.Sincos(r.Phis[i])
		r.PlaneXs[i], r.PlaneYs[i] = cos * r.Rs[i], sin * r.Rs[i]

		h.LineSegment(ring, i, ls)
		r.Xs[i] = r.Rs[i] * float64(ls.Dir[0])
		r.Ys[i] = r.Rs[i] * float64(ls.Dir[1])
		r.Zs[i] = r.Rs[i] * float64(ls.Dir[2])
	}
}

func (r *RingBuffer) OkPlaneCoords(xs, ys []float64) (okXs, okYs []float64) {
	xs, ys = xs[:0], ys[:0]
	for i, ok := range r.Oks {
		if ok {
			xs = append(xs, r.PlaneXs[i])
			ys = append(ys, r.PlaneYs[i])
		}
	}
	return xs, ys
}

func (r *RingBuffer) OkPolarCoords(rs, phis []float64) (okRs, okPhis []float64) {
	rs, phis = rs[:0], phis[:0]
	for i, ok := range r.Oks {
		if ok {
			rs = append(rs, r.Rs[i])
			phis = append(phis, r.Phis[i])
		}
	}
	return rs, phis
}
