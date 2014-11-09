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
	lowX, highX = 0, 2
	lowY, highY = 0, 2
	lowZ, highZ = 0, 2
)

type cmpIdx struct {
	slice, p int32
}

type ParticleManager struct {
	ps [][]tetra.Particle
	Locs map[int64]cmpIdx
	Size int64
}

func NewParticleManager() *ParticleManager {
	man := &ParticleManager{
		[][]tetra.Particle{},
		make(map[int64]cmpIdx),
		0,
	}

	return man
}

func (man *ParticleManager) Add(ps []tetra.Particle) {
	man.ps = append(man.ps, ps)

	for i := range ps {
		man.Locs[ps[i].Id] = cmpIdx{int32(len(man.ps) - 1), int32(i)}
	}

	man.Size += int64(len(ps))
}

func (man *ParticleManager) Get(id int64) *tetra.Particle {
	idx, ok := man.Locs[id]

	if !ok {
		return nil
	}
		
	return &man.ps[idx.slice][idx.p]
}

func main() {
	dir := os.Args[1]
	bins, err := strconv.Atoi(os.Args[2])
	if err != nil { panic(err) }

	ms := runtime.MemStats{}
	man := NewParticleManager()
	h0 := catalog.ReadHeader(path.Join(dir, "gridcell_0000.dat"))

	runtime.ReadMemStats(&ms)
	numRead := 0

	for x := lowX; x <= highX; x++ {
		for y := lowY; y <= highY; y++ {
			for z := lowZ; z <= highZ; z++ {
				numRead += 1
				_, ps := readParticles(int(h0.GridWidth), x, y, z, dir)
				fmt.Printf("%2d files read (%d part)\n", numRead, len(ps))
				runtime.ReadMemStats(&ms)
				man.Add(ps)
				fmt.Printf("    Alloc = %5d MB,   TotalAlloc = %5d MB\n",
					ms.Alloc >> 20, ms.TotalAlloc >> 20)
				runtime.ReadMemStats(&ms)
				
			}
		}
	}

	radii, fracs := BinFraction(h0, man, bins)
	fmt.Printf("# %15s %15s\n", "Fractions", "Radii")
	for i := range radii {
		fmt.Printf("  %15g %15g\n", fracs[i], radii[i])
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
	man *ParticleManager,
	bins int,
) (radii []float64, fracs []float64) {

	counts := make([]int, bins)
	contained := make([]int, bins)

	// Todo: make this wrap.
	xDiff := float64(highX - lowX + 1)
	yDiff := float64(highY - lowY + 1)
	zDiff := float64(highZ - lowZ + 1)

	maxDiff := triMax(xDiff, yDiff, zDiff)
	binWidth := h0.Width * (maxDiff * math.Sqrt(3.0) / 2.0) / float64(bins)

	nBuf := make([]int64, 3)
	cp := &tetra.Particle{}
	cp.Xs[0] = float32((xDiff / 2.0 + lowX) * h0.Width)
	cp.Xs[1] = float32((yDiff / 2.0 + lowY) * h0.Width)
	cp.Xs[2] = float32((zDiff / 2.0 + lowZ) * h0.Width)

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

	fmt.Println("#",counts)

	for i := range fracs {
		radii[i] = (float64(i) + 0.5) * binWidth
		fracs[i] = float64(contained[i]) / float64(counts[i])
	}

	return radii, fracs
}

func min(x, y float64) float64 {
	if x < y {
		return x
	}
	return y
}

func triMax(x, y, z float64) float64 {
	if x > y && x > z {
		return x
	} else if  y > z {
		return y
	}
	return z
}
