package analyze

import (
	"math"
	
	"github.com/phil-mansfield/gotetra/math/mat"
	"github.com/phil-mansfield/gotetra/los"

	"github.com/gonum/matrix/mat64"
)

func pinv(m, t *mat.Matrix) *mat.Matrix {
	// I HATE THIS
	gm := mat64.NewDense(m.Height, m.Width, m.Vals)
	gmt := mat64.NewDense(m.Width, m.Height, t.Vals)

	out1 := mat64.NewDense(m.Height, m.Height,
		make([]float64, m.Height * m.Height))
	out2 := mat64.NewDense(m.Width, m.Height,
		make([]float64, m.Height * m.Width))
	out1.Mul(gm, gmt)
	inv, err := mat64.Inverse(out1)
	if err != nil { panic(err.Error()) }
	out2.Mul(gmt, inv)

	vals := make([]float64, m.Width*m.Height)
	for y := 0; y < m.Width; y++ {
		for x := 0; x < m.Height; x++ {
			vals[y*m.Height + x] = out2.At(y, x)
		}
	}
	return mat.NewMatrix(vals, m.Height, m.Width)
}

func PennaCoeffs(xs, ys, zs []float64, I, J, K int) []float64 {
	N := len(xs)
	// TODO: Pass buffers to the function.
	rs := make([]float64, N)
	cosths := make([]float64, N)
	sinths := make([]float64, N)
	cosphis := make([]float64, N)
	sinphis := make([]float64, N)
	cs := make([]float64, I*J*K)

	// Precompute trig functions.
	for i := range rs {
		rs[i] = math.Sqrt(xs[i]*xs[i] + ys[i]*ys[i] + zs[i]*zs[i])
		cosths[i] = zs[i] / rs[i]
		sinths[i] = math.Sqrt(1 - cosths[i]*cosths[i])
		cosphis[i] = xs[i] / rs[i] / sinths[i]
		sinphis[i] = ys[i] / rs[i] / sinths[i]
	}

	MVals := make([]float64, I*J*K * len(xs))
	M := mat.NewMatrix(MVals, len(rs), I*J*K)

	// Populate matrix.
	for n := 0; n < N; n++ {
		m := 0
		for k := 0; k < K; k++ {
			costh := math.Pow(cosths[n], float64(k))
			for j := 0; j < J; j++ {
				sinphi := math.Pow(sinphis[n], float64(j))
				cosphi := 1.0
				for i := 0; i < I; i++ {
					MVals[m*M.Width + n] =
						math.Pow(sinths[n], float64(i+j)) *
						cosphi * costh * sinphi
					m++
					cosphi *= cosphis[n]
				}
			}
		}
	}

	// Solve.
	mat.VecMult(rs, pinv(M, M.Transpose()), cs)
	return cs
}

func PennaFunc(cs []float64, I, J, K int) Shell {
	return func(phi, th float64) float64 {
		idx, sum := 0, 0.0
		sinPhi, cosPhi := math.Sincos(phi)
		sinTh, cosTh := math.Sincos(th)

		for k := 0; k < K; k++ {
			cosK := math.Pow(cosTh, float64(k))
			for j := 0; j < J; j++ {
				sinJ := math.Pow(sinPhi, float64(j))
				for i := 0; i < I; i++ {
					cosI := math.Pow(cosPhi, float64(i))
					sinIJ := math.Pow(sinTh, float64(i+j))
					sum += cs[idx] * sinIJ * cosK * sinJ * cosI
					idx++
				}
			}
		}
		return sum
	}
}


func PennaVolumeFit(
	xs, ys [][]float64, h los.Halo, I, J int,
) (cs []float64, shell Shell) {
	n := 0
	for i := range xs { n += len(xs[i]) }
	fXs, fYs, fZs := make([]float64, n), make([]float64, n), make([]float64, n)
	
	idx := 0
	for i := range xs {
		for j := range xs[i] {
			fXs[idx], fYs[idx], fZs[idx] =
				h.PlaneToVolume(i, xs[i][j], ys[i][j])
			idx++
		}
	}

	cs = PennaCoeffs(fXs, fYs, fZs, I, J, 2)
	return cs, PennaFunc(cs, I, J, 2)
}

func PennaPlaneFit(
	xs, ys [][]float64, hRef *los.HaloProfiles, I, J int,
) (cs []float64, shell ProjectedShell) {
	n := 0
	for i := range xs { n += len(xs[i]) }
	fXs, fYs, fZs := make([]float64, n), make([]float64, n), make([]float64, n)
	
	idx := 0
	for i := range xs {
		for j := range xs[i] {
			fXs[idx], fYs[idx], fZs[idx] =
				hRef.PlaneToVolume(i, xs[i][j], ys[i][j])
			idx++
		}
	}

	cs = PennaCoeffs(fXs, fYs, fZs, I, J, 2)
	pf := PennaFunc(cs, I, J, 2)
	return cs, func (h *los.HaloProfiles, ring int, phi float64) float64 {
		sin, cos := math.Sincos(phi)
		x, y, z := h.PlaneToVolume(ring, cos, sin)
		pi2 := 2 * math.Pi
		return pf(math.Mod(math.Atan2(y, x) + pi2, pi2), math.Acos(z))
	}
}

func FilterPoints(
	rs []RingBuffer, levels int, hFactor float64,
) (pxs, pys [][]float64, ok bool) {
	pxs, pys = [][]float64{}, [][]float64{}
	
	for ri := range rs {
		r := &rs[ri]
		validXs := make([]float64, 0, r.N)
		validYs := make([]float64, 0, r.N)
		
		for i := 0; i < r.N; i++ {
			if r.Oks[i] {
				validXs = append(validXs, r.PlaneXs[i])
				validYs = append(validYs, r.PlaneYs[i])
			}
		}
		
		validRs, validPhis := []float64{}, []float64{}
		for i := range r.Rs {
			if r.Oks[i] {
				validRs = append(validRs, r.Rs[i])
				validPhis = append(validPhis, r.Phis[i])
			}
		}
		
		kt, ok := NewKDETree(validRs, validPhis, levels, hFactor)
		if !ok { return nil, nil, false }
		fRs, fThs, _ := kt.FilterNearby(validRs, validPhis, levels, kt.H() / 2)
		fXs, fYs := make([]float64, len(fRs)), make([]float64, len(fRs))
		for i := range fRs {
			sin, cos := math.Sincos(fThs[i])
			fXs[i], fYs[i] = fRs[i] * cos, fRs[i] * sin
		}

		pxs, pys = append(pxs, fXs), append(pys, fYs)
	}
	return pxs, pys, true
}
