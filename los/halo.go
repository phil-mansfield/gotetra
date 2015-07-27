package los

import (
	"fmt"

	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/math/mat"
)

type haloRing struct {
	ProfileRing
	rot, irot mat.Matrix32
}

func (hr *haloRing) Init(norm *geom.Vec, rMin, rMax float64, bins, n int) {
	hr.ProfileRing.Init(rMin, rMax, bins, n)
	
	zAxis := &geom.Vec{0, 0, 1}
	geom.EulerMatrixBetweenAt(norm, zAxis, &hr.rot)
	geom.EulerMatrixBetweenAt(zAxis, norm, &hr.irot)
}

type HaloProfiles struct {
	rs []haloRing
	origin geom.Vec
	rMin, rMax float64
	id int
}

func (hp *HaloProfiles) Init(
	id, rings int, origin *geom.Vec, rMin, rMax float64,
	bins, n int,
) *HaloProfiles {
	
	solid, ok := geom.NewUniquePlatonicSolid(rings)
	if !ok {
		panic(fmt.Sprintf("Cannot uniformly space %d rings.", rings))
	}

	hp.rs = make([]haloRing, rings)
	hp.origin = *origin
	hp.rMin, hp.rMax = rMin, rMax
	hp.id = id

	norms := solid.UniqueNormals()
	for i := 0; i < rings; i++ {
		hp.rs[i].Init(&norms[i], rMin, rMax, bins, n)
	}

	return hp
}
