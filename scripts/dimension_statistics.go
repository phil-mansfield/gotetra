package main
import (
	"fmt"
	"math"
	"strconv"
	"os"

	tetra "github.com/phil-mansfield/gotetra"
	"github.com/phil-mansfield/gotetra/scripts/helper"
)

const (
	midX, midY, midZ = 0, 0, 1
	layers = 1
)

func main() {
	dir := os.Args[1]
	bins, err := strconv.Atoi(os.Args[2])
	if err != nil { panic(err) }

	h0, man, centerPs := helper.ReadCatalogs(dir, midX, midY, midZ, layers)

	lRs, lRFreq, lVs, lVFreq, f := dimensionBins(h0, centerPs, man, bins)

	fmt.Printf("# Center location: (%d, %d, %d)\n", midX, midY, midZ)
	fmt.Printf("# Miss Frac: %g\n", f)
	fmt.Printf("# %12s %12s %12s %12s\n", "R",
		"N/Ntot/dlR", "V", "N/Ntot/dlV")
	for i := range lRs {
		fmt.Printf("  %12.6g %12.6g %12.6g %12.6g\n",
			lRs[i], lRFreq[i], lVs[i], lVFreq[i])
	}
}

// returns binned dimensional statistics
func dimensionBins(
	h0 *tetra.Header,
	ps []tetra.Particle,
	man *tetra.ParticleManager,
	bins int,
) (lRs, lRFreqs, lVs, lVFreqs []float64, missFrac float64) {
	Vs := make([]float64, 6 * len(ps))
	Rs := make([]float64, 6 * 3 * len(ps))

	
	cornerBuf := []int64{ 0, 0, 0 }
	partBuf := []*tetra.Particle{ nil, nil, nil, nil }
	vb := tetra.NewVolumeBuffer()

	rIdx, vIdx := 0, 0

	for i := range ps {
		for dir := 0; dir < 6; dir++ {
			h0.TetraCorners(ps[i].Id, dir, cornerBuf)
			partBuf[3] = &ps[i]

			for j, cornerIdx := range cornerBuf {
				c := man.Get(cornerIdx)
				partBuf[j] = c
				
				if c != nil {
					Rs[rIdx] = h0.Distance(&ps[i], c)
					rIdx++
				}
			}
			
			if partBuf[0] != nil && partBuf[1] != nil && partBuf[2] != nil {
				vol := h0.Volume(partBuf, vb)
				Vs[vIdx] = vol
				vIdx++
			}
		}
	}


	lRs, lRFreqs = logBin(Rs[0: rIdx], bins)
	lVs, lVFreqs = logBin(Vs[0: vIdx], bins)
	missFrac = float64(vIdx) / float64(len(Vs))

	return lRs, lRFreqs, lVs, lVFreqs, missFrac
}

// bins up the input slice logarithmically
func logBin(xs []float64, bins int) (lXs, lXFreqs []float64) {
	minX, maxX := xs[0], xs[0]
	for _, x := range xs {
		if x < minX { minX = x }
		if x > maxX { maxX = x }
	}

	counts := make([]int64, bins)
	binWidth := math.Log10(maxX / minX) / float64(bins)

	for _, x := range xs {
		idx := int(math.Floor(math.Log10(x / minX) / binWidth))
		if idx == len(counts) { idx = len(counts) - 1 }

		if idx > len(counts) || idx < 0 {
			fmt.Println(x, math.Log10(x), idx)
		}

		counts[idx]++
	}

	lXs = make([]float64, bins)
	lXFreqs = make([]float64, bins)

	for i := range counts {
		lXs[i] = math.Pow(10.0, (float64(i) + 0.5) * binWidth) * minX
		lXFreqs[i] = float64(counts[i]) / float64(len(xs)) / binWidth
	}

	return lXs, lXFreqs
}
