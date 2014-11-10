package main
import (
	"fmt"
	"math"
	"strconv"
	"os"
	"path"
	"runtime"

	tetra "github.com/phil-mansfield/gotetra"
	"github.com/phil-mansfield/gotetra/catalog"
)

const (
	midX, midY, midZ = 0, 6, 5
	layers = 1
)

func main() {
	dir := os.Args[1]
	bins, err := strconv.Atoi(os.Args[2])
	if err != nil { panic(err) }

	ms := runtime.MemStats{}
	man := tetra.NewParticleManager()
	h0 := catalog.ReadHeader(path.Join(dir, "gridcell_0000.dat"))

	runtime.ReadMemStats(&ms)

	for x := midX - layers; x <= midX + layers; x++ {
		for y := midY - layers; y <= midY + layers; y++ {
			for z := midZ - layers; z <= midZ + layers; z++ {
				xIdx := (x + int(h0.GridWidth)) % int(h0.GridWidth)
				yIdx := (y + int(h0.GridWidth)) % int(h0.GridWidth)
				zIdx := (z + int(h0.GridWidth)) % int(h0.GridWidth)
					
				_, ps := readParticles(int(h0.GridWidth), xIdx, yIdx, zIdx, dir)
				man.Add(ps)
				runtime.ReadMemStats(&ms)
			}
		}
	}

	radii, fracs, counts := BinFraction(h0, man, bins)
	fmt.Printf("# Center Cell: (%d, %d, %d)\n", midX, midY, midZ)
	fmt.Printf("# %15s %15s %15s\n", "Fractions", "Radii", "Count")
	for i := range radii {
		fmt.Printf("  %15.6g %15.6g %15d\n",
			fracs[i], radii[i], counts[i])
	}
}

func readParticles(
	gridWidth, x, y, z int,
	dir string,
) (*tetra.Header, []tetra.Particle) {
	idx := x + y * gridWidth + z * gridWidth * gridWidth
	name := fmt.Sprintf("gridcell_%04d.dat", idx)
	path := path.Join(dir, name)
	return catalog.Read(path)
}

func BinFraction(
	h0 *tetra.Header,
	man *tetra.ParticleManager,
	bins int,
) (radii []float64, fracs []float64, counts []int) {

	counts = make([]int, bins)
	contained := make([]int, bins)

	binWidth := (h0.Width * math.Sqrt(3) *
		(float64(layers) + 0.5) / float64(bins))

	nBuf := make([]int64, 3)
	cp := &tetra.Particle{}
	cp.Xs[0] = float32((float64(midX) + 0.5) * h0.Width)
	cp.Xs[1] = float32((float64(midY) + 0.5) * h0.Width)
	cp.Xs[2] = float32((float64(midZ) + 0.5) * h0.Width)

	for id := range man.Locs {
 		p := man.Get(id)
		if p == nil { panic("Oh god.") }
		h0.TetraCorners(id, nBuf)

		binIdx := int(h0.Distance(p, cp) / binWidth)

		for _, n := range nBuf {
			np := man.Get(n)
			counts[binIdx]++
			if np != nil { contained[binIdx]++ }
		}
	}

	radii = make([]float64, bins)
	fracs = make([]float64, bins)

	for i := range fracs {
		radii[i] = (float64(i) + 0.5) * binWidth
		fracs[i] = float64(contained[i]) / float64(counts[i])
	}

	return radii, fracs, counts
}
