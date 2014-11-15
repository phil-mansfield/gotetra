package main

import (
	"fmt"
	"os"
	"strconv"

	tetra "github.com/phil-mansfield/gotetra"
	"github.com/phil-mansfield/gotetra/scripts/helper"
)

type Buffers struct {
	Id []int64
	T tetra.Tetra
	Bounds tetra.Bounds
	Vol *tetra.VolumeBuffer
}

func NewBuffers() *Buffers {
	return &Buffers{
		[]int64{0, 0, 0, 0},
		[]*tetra.Particle{nil, nil, nil, nil},
		tetra.Bounds{},
		tetra.NewVolumeBuffer(),
	}
}

type Grid3D struct {
	Xs []float64 // Maybe Grid3Ds shouldn't manage this themselves...
	Width int
	X0, Y0, Z0 int
}

func NewGrid3D(x0, y0, z0, width int) *Grid3D {
	g := &Grid3D{}
	g.Xs = make([]float64, width * width * width)
	g.X0, g.Y0, g.Z0 = x0, y0, z0
	g.Width = width

	return g
}

func (g *Grid3D) Idx(x, y, z int) int {
	x -= g.X0
	y -= g.Y0
	z -= g.Z0

	return x + y * g.Width + z * g.Width * g.Width
}

func (g *Grid3D) Coords(idx int) (x, y, z int) {
	x = idx % g.Width
	y = idx % (g.Width * g.Width) / g.Width
	z = idx / (g.Width * g.Width)
	return x, y, z
}

func (g *Grid3D) Bounds() *tetra.Bounds {
	b := &tetra.Bounds{}
	b.MinX, b.MaxX = g.X0, g.X0 + g.Width - 1
	b.MinY, b.MaxY = g.Y0, g.Y0 + g.Width - 1
	b.MinZ, b.MaxZ = g.Z0, g.Z0 + g.Width - 1
	return b
}

func main() {
	dir, xStr, yStr, zStr := os.Args[1], os.Args[2], os.Args[3], os.Args[4]
	cellStrs := os.Args[5]

	x, err := strconv.Atoi(xStr)
	if err != nil { panic(err) }
	y, err := strconv.Atoi(yStr)
	if err != nil { panic(err) }
	z, err := strconv.Atoi(zStr)
	if err != nil { panic(err) }
	cells, err := strconv.Atoi(cellStrs)
	if err != nil { panic(err) }

	h, man, ps := helper.ReadCatalogs(dir, x, y, z, 1)

	hGrid := NewGrid3D(0, 0, 0, int(h.GridWidth))
	hX, hY, hZ := hGrid.Coords(int(h.Idx))

	rhoGrid := NewGrid3D(hX * cells, hY * cells, hZ * cells, cells)
	rhoBounds := rhoGrid.Bounds()
	rhos := rhoGrid.Xs

	bufs := NewBuffers()
	cellWidth := h.Width / float64(cells)
	successes, failures := 0, 0

	for i := range ps {
		if i % 100000 == 0 && i > 0 {
			fmt.Printf("# steps = %d\n", i)
		}

		for dir := 0; dir < 6; dir++ {
			h.TetraCorners(ps[i].Id, dir, bufs.Id)

			bufs.T[0] = &ps[i]
			for i := 0; i < 3; i++ { bufs.T[i + 1] = man.Get(bufs.Id[i]) }
			if !bufs.T.Valid() { continue }

			bufs.T.Bounds(cellWidth, &bufs.Bounds)
			ds, df := addVol(rhoGrid, rhoBounds, h, bufs)
			successes += ds
			failures += df
		}
	}

	fmt.Printf("# catalog: (%d %d %d)\n", hX, hY, hZ)
	fmt.Println("#", successes, failures)
	fmt.Println("# cells =", cells)

	for _, rho := range rhos {
		fmt.Println(rho * h.Mass)
	}
}

func addVol(grid *Grid3D, b1 *tetra.Bounds, h *tetra.Header, bufs *Buffers) (int, int) {
	b2 := &bufs.Bounds
	minX, maxX := max(b1.MinX, b2.MinX), min(b1.MaxX, b2.MaxX)
	minY, maxY := max(b1.MinY, b2.MinY), min(b1.MaxY, b2.MaxY)
	minZ, maxZ := max(b1.MinZ, b2.MinZ), min(b1.MaxZ, b2.MaxZ)

	p := &tetra.Particle{ }
	vol := h.Volume(bufs.T, bufs.Vol)
	rho := 1.0 / vol

	successes, failures := 0, 0
	cellWidth := h.Width / float64(grid.Width)

	for x := minX; x <= maxX; x++ {
		for y := minY; y <= maxY; y++ {
			for z := minZ; z <= maxZ; z++ {
				p.Xs[0] = float32((float64(x) + 0.5) * cellWidth)
				p.Xs[1] = float32((float64(y) + 0.5) * cellWidth)
				p.Xs[2] = float32((float64(z) + 0.5) * cellWidth)

				if h.WithinTetra(p, bufs.T, vol, bufs.Vol) {
					grid.Xs[grid.Idx(x, y, z)] += rho
					successes++
				} else {
					failures++
				}
			} 
		}
	}
	return successes, failures
}

func printNormXs(p *tetra.Particle, cellWidth float64) {
	for i := 0; i < 3; i++ {
		fmt.Printf("%.5g ", float64(p.Xs[i]) / cellWidth)
	}
	fmt.Println()
}

func max(x, y int) int {
	if x > y { return x }
	return y
}

func min(x, y int) int {
	if x < y { return x }
	return y
}
