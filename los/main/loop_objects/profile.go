package loop_objects

import (
	"math"

	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/los/main/gtet_util/loop"
)

type Profile struct {
	loop.Sphere // Handles intersection checks.
	R0 [3]float64
	RMin, RMax float64
	Counts []float64
	
	dlr, lrMin, rMin2, rMax2 float64
}

// type-checking
var (
	_ loop.Object = &Profile{}
)

func NewProfile(
	r0 [3]float64, rMin, rMax float64, rBins int,
) *Profile {
	p := &Profile{}

	p.Sphere.Init(r0, rMin, rMax)
	p.R0 = r0
	p.RMin = rMin
	p.RMax = rMax

	p.Counts = make([]float64, rBins)
	p.dlr = math.Log(p.RMax) - math.Log(p.RMin)
	p.lrMin = math.Log(p.RMin)
	p.rMin2 = p.RMin * p.RMin
	p.rMax2 = p.RMax * p.RMax

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
	panic("UseTetra called on Profile even though" +
		" loop.Params.UsesTetra == false.")
}

func (p *Profile) Params() loop.Params {
	return loop.Params{
		UsesTetra: false,
	}
}
