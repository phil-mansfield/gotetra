package main

import (
	"fmt"
	"path"
	"os"

	"runtime/pprof"
	"unsafe"
	"flag"

	tetra "github.com/phil-mansfield/gotetra"
	"github.com/phil-mansfield/gotetra/catalog"
)

const (
	lowX, highX = 0, 0
	lowY, highY = 0, 0
	lowZ, highZ = 0, 0
)

type managerLocation struct {
	idx, slice int32
}

type ParticleManager struct {
	ps [][]tetra.Particle
	locs map[int64]managerLocation
	Size int64
}

func NewParticleManager() *ParticleManager {
	man := &ParticleManager{
		[][]tetra.Particle{},
		make(map[int64]managerLocation),
		0,
	}

	return man
}

func (man *ParticleManager) Add(ps []tetra.Particle) {
	man.ps = append(man.ps, ps)
	for i := range ps {
		man.locs[ps[i].Id] = managerLocation{int32(i), int32(len(man.ps) - 1)}
	}

	man.Size += int64(len(ps))
}

func (man *ParticleManager) Get(id int64) *tetra.Particle {
	loc, ok := man.locs[id]

	if !ok {
		return nil
	}

	return &man.ps[loc.slice][loc.idx]
}

var mprof = flag.String("mprof", "", "write memory profile to this file")

func main() {
	flag.Parse()

	N := 1<<15
	psSlice := make([][]int64, N)
	_ = psSlice
	size := 0

	for i := 0; i < N; i++ {
		ps := make([]int64, 1<<12)
		psSlice[i] = ps
		size += int(unsafe.Sizeof(ps[0])) * len(ps)
	}


	if *mprof != "" {
		f, err := os.Create(*mprof)
		if err != nil { panic(err) }
		pprof.WriteHeapProfile(f)
		f.Close()
	}

	fmt.Printf("total allocated: %d MB\n", size >> 20)

	// man := NewParticleManager()

	/*
	for x := lowX; x <= highX; x++ {
		for y := lowY; y <= highY; y++ {
			for z := lowZ; z <= highZ; z++ {
				_, ps := readParticles(int(h0.GridWidth), x, y, z, dir)
				// man.Add(ps)
				fmt.Printf("%d -> %d Mb\n", len(ps),
					(len(ps) * (8*3 + 8*3 + 8 + 8)) >> 20)
				fmt.Println(ps[0].Id)
			}
		}
	}
    */
}

func readParticles(gridWidth, x, y, z int, dir string) (*tetra.Header, []tetra.Particle) {
	idx := x + y * gridWidth + z * gridWidth * gridWidth
	name := fmt.Sprintf("gridcell_%04d.dat", idx)
	fmt.Println(name)
	path := path.Join(dir, name)
	return catalog.Read(path)
}
