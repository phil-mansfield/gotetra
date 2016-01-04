package loop_objects

import (
	"math"

	"github.com/phil-mansfield/gotetra/math/rand"
	"github.com/phil-mansfield/gotetra/los/geom"
	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/los/main/gtet_util/loop"
)

type Profile struct {
	loop.Sphere // Handles intersection checks.
	R0 [3]float64
	RMin, RMax float64
	Counts []float64
	
	dlr, lrMin, rMin2, rMax2 float64

	vecBuf []rgeom.Vec
	randBuf []float64
	gen *rand.Generator
}

// type-checking
var (
	_ loop.Object = &Profile{}
)

func NewParticleProfile(
	r0 [3]float64, rMin, rMax float64, rBins int,
) *Profile {
	return newProfile(r0, rMin, rMax, rBins, 0)
}

func NewTetraProfile(
	r0 [3]float64, rMin, rMax float64, rBins, pts int,
) *Profile {
	return newProfile(r0, rMin, rMax, rBins, pts)
}

func newProfile(
	r0 [3]float64, rMin, rMax float64, rBins, pts int,
) *Profile {
	p := &Profile{}

	p.Sphere.Init(r0, rMin, rMax)
	p.R0 = r0
	p.RMin = rMin
	p.RMax = rMax

	p.Counts = make([]float64, rBins)
	p.dlr = (math.Log(p.RMax) - math.Log(p.RMin)) / float64(len(p.Counts))
	p.lrMin = math.Log(p.RMin)
	p.rMin2 = p.RMin * p.RMin
	p.rMax2 = p.RMax * p.RMax

	tetPts := int(math.Ceil(float64(pts*pts*pts) / 6))
	p.vecBuf = make([]rgeom.Vec, tetPts)
	p.randBuf = make([]float64, 3*tetPts)
	p.gen = rand.NewTimeSeed(rand.Xorshift)
	
	return p
}

func (p *Profile) ThreadCopy(id, threads int) loop.Object {
	panic("NYI")
}

func (p *Profile) ThreadMerge(objs []loop.Object) {
	panic("NYI")
}

func (p *Profile) UsePoint(x, y, z float64) {
	x0, y0, z0 := p.R0[0], p.R0[1], p.R0[2]
	dx, dy, dz := x - x0, y - y0, z - z0
	r2 := dx*dx + dy*dy + dz*dz
	if r2 <= p.rMin2 || r2 >= p.rMax2 { return }
	lr := math.Log(r2) / 2
	ir := int(((lr) - p.lrMin) / p.dlr)
	p.Counts[ir]++
}
	
func (p *Profile) UseTetra(t *geom.Tetra) {
	tet := &rgeom.Tetra{}
	t0, t1 := rgeom.Vec(t[0]), rgeom.Vec(t[1])
	t2, t3 := rgeom.Vec(t[2]), rgeom.Vec(t[3])
	tet.Init(&t0, &t1, &t2, &t3)
	tet.RandomSample(p.gen, p.randBuf, p.vecBuf)
	
	for _, pt := range p.vecBuf {
		x, y, z := float64(pt[0]), float64(pt[1]), float64(pt[2])
		p.UsePoint(x, y, z)
	}
}

	
func (p *Profile) Params() loop.Params {
	return loop.Params{
		UsesTetra: len(p.vecBuf) != 0,
	}
}
