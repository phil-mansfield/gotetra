package loop_objects

import (
	"math"
	gorand "math/rand"

	"github.com/phil-mansfield/gotetra/math/rand"
	"github.com/phil-mansfield/gotetra/los/geom"
	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/los/main/gtet_util/loop"
)

type Renderer struct {
	loop.Box
	Origin [3]float64
	Pixels [3]int
	Pw float64
	Counts []int

	vecBuf []rgeom.Vec
	randBuf []float64
	gen *rand.Generator
}

// type-checking
var _ loop.Object = &Renderer{}

func NewParticleRenderer(
	origin [3]float64, pixels [3]int, pw float64,
) *Renderer {
	return newRenderer(origin, pixels, pw, 0)
}

func NewTetraRenderer(
	origin [3]float64, pixels [3]int, pw float64, pts int,
) *Renderer {
	return newRenderer(origin, pixels, pw, pts)
}

func newRenderer(
	origin [3]float64, pixels [3]int, pw float64, pts int,
) *Renderer {
	r := &Renderer{}
	span := [3]float64{
		float64(pixels[0]) * pw,
		float64(pixels[1]) * pw,
		float64(pixels[2]) * pw,
	}
	r.Box.Init(origin, span)
	r.Origin = origin
	r.Pixels = pixels
	r.Pw = pw

	r.Counts = make([]int, pixels[0]*pixels[1]*pixels[2])

    tetPts := int(math.Ceil(float64(pts*pts*pts) / 6))
    r.vecBuf = make([]rgeom.Vec, tetPts)
    r.randBuf = make([]float64, 3*tetPts)
    r.gen = rand.NewTimeSeed(rand.Xorshift)

	return r
}

func (r *Renderer) ThreadCopy(id, threads int) loop.Object {
	panic("NYI")
}

func (r *Renderer) ThreadMerge(objs []loop.Object) {
	panic("NYI")
}

func (r *Renderer) UsePoint(x, y, z float64) {
	if !r.Contains(x, y, z) { return }
	ix := int((x - r.Origin[0]) / r.Pw)
	iy := int((y - r.Origin[1]) / r.Pw)
	iz := int((z - r.Origin[2]) / r.Pw)

	r.Counts[iz*r.Pixels[1]*r.Pixels[0] + iy*r.Pixels[0] + ix]++
}

func (r *Renderer) UseTetra(t *geom.Tetra) {
    tet := &rgeom.Tetra{}
    t0, t1 := rgeom.Vec(t[0]), rgeom.Vec(t[1])
    t2, t3 := rgeom.Vec(t[2]), rgeom.Vec(t[3])
    tet.Init(&t0, &t1, &t2, &t3)
    tet.RandomSample(r.gen, r.randBuf, r.vecBuf)

    for _, pt := range r.vecBuf {
        x, y, z := float64(pt[0]), float64(pt[1]), float64(pt[2])
        r.UsePoint(x, y, z)
    }
}

func (r *Renderer) Params() loop.Params {
    return loop.Params{
        UsesTetra: len(r.vecBuf) != 0,
    }	
}

////////////////////////////
// Kernel Implementations //
////////////////////////////

type BallKernelRenderer struct {
	Renderer
	KernelR float64
	KernelPts int
	origin, span [3]float64
}

type GaussianKernelRenderer struct {
	Renderer
	KernelR, kR2 float64
	KernelPts int
	origin, span [3]float64
}

func NewBallKernelRenderer(
	origin [3]float64, pixels [3]int, pw float64,
	kernelR float64, kernelPts int,
) *BallKernelRenderer {
	r := &BallKernelRenderer{}
	r.Renderer = *NewParticleRenderer(origin, pixels, pw)
	r.KernelR = kernelR
	r.KernelPts = kernelPts
	span := [3]float64{
		pw*float64(pixels[0]), pw*float64(pixels[1]), pw*float64(pixels[2]),
	}
	r.span, r.origin = span, origin
	for i := 0; i < 3; i++ {
		r.span[i] += kernelR*2
		r.origin[i] -= kernelR
	}
	return r
}

func (r *BallKernelRenderer) UsePoint(x, y, z float64) {
	for i := 0; i < r.KernelPts; i++ {
		rx, ry, rz := randomBallPoint(r.KernelR)
		r.Renderer.UsePoint(x+rx, y+ry, z+rz)
	}
}

func randomBallPoint(r float64) (x, y, z float64) {
	x, y, z = gorand.Float64()*2-1, gorand.Float64()*2-1, gorand.Float64()*2-1
	if x*x + y*y + z*z > 1 { return randomBallPoint(r) }
	return x*r, y*r, z*r
}

func (r *BallKernelRenderer) Contains(x, y, z float64) bool {
    lowX, highX := r.origin[0], r.origin[0] + r.span[0]
    lowY, highY := r.origin[1], r.origin[1] + r.span[1]
    lowZ, highZ := r.origin[2], r.origin[2] + r.span[2]
    return lowX < x && x < highX &&
        lowY < y && y < highY &&
        lowZ < z && z < highZ
}

func NewGaussianKernelRenderer(
	origin [3]float64, pixels [3]int, pw float64,
	kernelR float64, kernelPts int,
) *GaussianKernelRenderer {
	r := &GaussianKernelRenderer{}
	r.Renderer = *NewParticleRenderer(origin, pixels, pw)
	r.KernelR = kernelR
	r.KernelPts = kernelPts
	span := [3]float64{
		pw*float64(pixels[0]), pw*float64(pixels[1]), pw*float64(pixels[2]),
	}
	r.span, r.origin = span, origin
	for i := 0; i < 3; i++ {
		r.span[i] += kernelR*2
		r.origin[i] -= kernelR
	}
	return r
}

func (r *GaussianKernelRenderer) Contains(x, y, z float64) bool {
    lowX, highX := r.origin[0], r.origin[0] + r.span[0]
    lowY, highY := r.origin[1], r.origin[1] + r.span[1]
    lowZ, highZ := r.origin[2], r.origin[2] + r.span[2]
    return lowX < x && x < highX &&
        lowY < y && y < highY &&
        lowZ < z && z < highZ
}

func (r *GaussianKernelRenderer) UsePoint(x, y, z float64) {
	panic("NYI")
}
