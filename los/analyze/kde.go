package analyze

import (
	"math"

	intr "github.com/phil-mansfield/gotetra/math/interpolate"
)

func GaussianKDE(xs []float64, h, low, high float64, n int) *intr.Spline {
	dx := (high - low) / float64(n - 1)
	spXs, spYs := make([]float64, n), make([]float64, n)
	for i := 0; i < n-1; i++ { spXs[i] = low + dx*float64(i) }
	spXs[n - 1] = high

	maxDist := h * 3

	for _, x:= range xs {
		lowIdx := int((x - maxDist - low) / dx)
		highIdx := int((x + maxDist -low) / dx) + 1
		if lowIdx < 0 { lowIdx = 0 }
		if highIdx >= n { highIdx = n - 1 }
		for i := lowIdx; i <= highIdx; i++ {
			udx := (spXs[i] - x) / h
			spYs[i] += math.Exp(-udx*udx)
		}
	}

	return intr.NewSpline(spXs, spYs)
}

type KDETree struct {
	h, low, high float64
	spTree [][]*intr.Spline
	maxesTree [][][]float64
	thTree, connMaxes [][]float64
	spRs []float64
}

func NewKDETree(rs, phis []float64, splits int) *KDETree {
	kt := new(KDETree)

	hFactor := 10.0
	rn := 100

	kt.low, kt.high = 0, rs[0]
	for _, r := range rs {
		if r > kt.high { kt.high = r }
	}

	kt.spRs = make([]float64, rn)
	dr := (kt.high - kt.low) / float64(rn)
	for i := range kt.spRs {
		kt.spRs[i] = kt.low + dr*float64(i)
	}
	kt.spRs[rn - 1] = kt.high

	kt.h = (kt.high - kt.low)/ hFactor
	kt.spTree = [][]*intr.Spline{{GaussianKDE(rs, kt.h, kt.low, kt.high, 100)}}
	kt.thTree = [][]float64{{math.Pi}}

	kt.growTrees(rs, phis, splits)
	kt.findMaxes()
	kt.connectMaxes()

	return kt
}

func (kt *KDETree) growTrees(rs, phis []float64, splits int) {
	for split := 0; split < splits; split++ {
		bins := int(1 << uint((1 + split)))
		rBins, thBins := binByTheta(rs, phis, bins)
		sps := make([]*intr.Spline, bins)

		for i, rBin := range rBins {
			sps[i] = GaussianKDE(rBin, kt.h, kt.low, kt.high, 100)
		}
		kt.thTree = append(kt.thTree, thBins)
		kt.spTree = append(kt.spTree, sps)
	}
}

func binByTheta(
	rs, ths []float64, bins int,
) (rBins [][]float64, thBins []float64) {
	rBins = make([][]float64, bins)
	for i := range rBins { rBins[i] = []float64{} }
	dth := (2 * math.Pi) / float64(bins)
	for i := range rs {
		idx := int(ths[i] / dth)
		rBins[i] = append(rBins[idx], rs[i])
	}

	thBins = make([]float64, bins)
	for i := range thBins {
		thBins[i] = 2 * math.Pi * (float64(i) + 0.5) / float64(bins)
	}

	return rBins, thBins
}

func (kt *KDETree) findMaxes() {
	kt.maxesTree = [][][]float64{}
	for _, sps := range kt.spTree {
		levelMaxes := make([][]float64, len(sps))
		for j, sp := range sps {
			maxes := localSplineMaxes(kt.spRs, sp)
			levelMaxes[j] = maxes
		}

		kt.maxesTree = append(kt.maxesTree, levelMaxes)
	}
}

func localSplineMaxes(xs []float64, sp *intr.Spline) []float64 {
	prev, curr, next := sp.Eval(xs[0]), sp.Eval(xs[1]), sp.Eval(xs[2])
	maxes := []float64{}
	if curr > next && curr > prev { maxes = append(maxes, xs[1]) }
	for i := 2; i < len(xs); i++ {
		prev, curr, next = curr, next, xs[i]
		if curr > next && curr > prev { maxes = append(maxes, xs[1]) }
	}
	return maxes
}

func (kt *KDETree) connectMaxes() {
	kt.connMaxes = [][]float64{{kt.maxesTree[0][0][0]}}

	for split, maxes := range kt.maxesTree[1:] {
		prevMaxes := kt.connMaxes[len(kt.connMaxes) - 1]
		currMaxes := make([]float64, 2 * len(prevMaxes))
		for i := range currMaxes { currMaxes[i] = math.NaN() }

		for node := range maxes {
			nodeMaxes := maxes[node]
			if len(maxes) == 0 { continue }
			nodePrevMax := prevMaxes[node / 2]
			
			var connMax float64
			if math.IsNaN(nodePrevMax) {
				connMax = math.NaN()
			} else {
				connIdx, maxDist := -1, -1.0
				for i := range nodeMaxes {
					dist := math.Abs(nodePrevMax - nodeMaxes[connIdx])
					if dist > maxDist { connIdx, maxDist = i, dist }
				}
				connMax = nodeMaxes[connIdx]
			}
			
			if !math.IsNaN(connMax) && (split == 0 ||
				math.Abs(connMax - nodePrevMax) < kt.h) {
				currMaxes[node] = connMax
			} else {
				for _, max := range nodeMaxes {
					rFunc := kt.GetRFunc(split)
					spR := rFunc(kt.thTree[split+1][node])
					if math.Abs(max - spR) < kt.h { currMaxes[node] = max }
				}
			}

			kt.connMaxes = append(kt.connMaxes, currMaxes)
		}
	}
}

func (kt *KDETree) getFinestMax(idx, level int) float64 {
	for i := 0; i < level; i++ {
		r := kt.connMaxes[level - i][idx / (1 << uint(i))]
		if !math.IsNaN(r) { return r }
	}
	panic(":3")
}

func (kt *KDETree) GetConnMaxes(level int) (rs, ths []float64) {
	ths = kt.thTree[level]
	maxes := make([]float64, len(kt.connMaxes[level]))
	for i := range maxes {
		maxes[i] = kt.getFinestMax(i, level)
	}
	return maxes, ths
}

func extendAngularRange(maxes, ths []float64) (spMaxes, spThs []float64) {
	n := len(maxes)
	buf := 5
	spThs, spMaxes = make([]float64, 2*buf + n), make([]float64, 2*buf + n)

	j := n - buf
	for i := 0; i < buf; i++ {
		spThs[i], spMaxes[i] = ths[j] + 2*math.Pi, maxes[j]
		j++
	}
	j = 0
	for i := buf; i < n + buf; i++ {
		spThs[i], spMaxes[i] = ths[j], maxes[j]
		j++
	}
	j = 0
	for i := n + buf; i < 2*n + buf; i++ {
		spThs[i], spMaxes[i] = ths[j] - 2*math.Pi, maxes[j]
		j++
	}
	return spMaxes, spThs
}

func (kt *KDETree) GetRFunc(level int) (func(float64) float64) {
	maxes, ths := kt.GetConnMaxes(level)
	spMaxes, spThs := extendAngularRange(maxes, ths)
	spXs, spYs := make([]float64, len(spThs)), make([]float64, len(spThs))
	for i, th := range spThs {
		sin, cos := math.Sincos(th)
		spXs[i], spYs[i] = spMaxes[i] * cos, spMaxes[i] * sin
	}
	xSp, ySp := intr.NewSpline(spXs, spThs), intr.NewSpline(spYs, spThs)
	return func(th float64) float64 {
		x, y := xSp.Eval(th), ySp.Eval(th)
		return math.Sqrt(x*x + y*y)
	}
}

func (kt *KDETree) FilterNearby(
	rs, ths []float64, level int, dr float64,
) (fRs, fThs []float64, idxs []int) {
	rFunc := kt.GetRFunc(level)
	fRs, fThs, idxs = []float64{}, []float64{}, []int{}
	for i := range rs {
		if math.Abs(rFunc(ths[i]) - rs[i]) < dr {
			fRs = append(fRs, rs[i])
			fThs = append(fThs, ths[i])
			idxs = append(idxs, i)
		}
	}
	return fRs, fThs, idxs
}
