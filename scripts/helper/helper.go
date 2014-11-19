package helper

import (
	"fmt"
	"path"

    "github.com/phil-mansfield/gotetra/catalog"
)

func ReadCatalogs(
	dir string, midX, midY, midZ, layers int,
) (*catalog.Header, *catalog.ParticleManager, []catalog.Particle) {
	man := catalog.NewParticleManager()
	h0 := catalog.ReadHeader(path.Join(dir, "gridcell_0000.dat"))

	var centerPs []catalog.Particle

	for x := midX - layers; x <= midX + layers; x++ {
		for y := midY - layers; y <= midY + layers; y++ {
			for z := midZ - layers; z <= midZ + layers; z++ {
				xIdx := (x + int(h0.GridWidth)) % int(h0.GridWidth)
				yIdx := (y + int(h0.GridWidth)) % int(h0.GridWidth)
				zIdx := (z + int(h0.GridWidth)) % int(h0.GridWidth)
					
				h, ps := readParticles(
					int(h0.GridWidth), xIdx, yIdx, zIdx, dir,
				)

				man.Add(ps)

				if x == midX && y == midY && z == midZ {
					centerPs = ps
					h0 = h
				}
			}
		}
	}

	return h0, man, centerPs
}

func readParticles(
	gridWidth, x, y, z int,
	dir string,
) (*catalog.Header, []catalog.Particle) {
	idx := x + y * gridWidth + z * gridWidth * gridWidth
	name := fmt.Sprintf("gridcell_%04d.dat", idx)
	path := path.Join(dir, name)
	return catalog.Read(path)
}
