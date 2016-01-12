package main

import (
	"fmt"
	"log"
	"flag"
	"os"
	"runtime"

	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/render/density"
	"github.com/phil-mansfield/gotetra/render/io"

	obj "github.com/phil-mansfield/gotetra/los/main/loop_objects"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
	"github.com/phil-mansfield/gotetra/los/main/gtet_util/loop"
)

func main() {
	log.Println("gtet_render")

	p := parseCmd()
	ids, snaps, _, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }
	if len(ids) == 0 { log.Fatal("No input IDs.") }

	fnames, err := computeGrids(ids, snaps, p)
	if err != nil { log.Fatal(err.Error()) }


	for _, fname := range fnames { fmt.Println(fname) }
}

type Params struct {
	Pixels int
	SubsampleLength int
	Mult float64

	Sphere float64
	Tetra, Linear, Cubic, SpherePts int
}

// parseCmd parses the command line options and returns them as a struct.
func parseCmd() *Params {
	p := &Params{}

	flag.IntVar(&p.Pixels, "Pixels", 100,
		"Number of pixels on one side of render.")
	flag.IntVar(&p.SubsampleLength, "SubsampleLength", 1,
		"Grid distance between adjacent points used in interpolation.")
	flag.Float64Var(&p.Mult, "Mult", 3, "Size of box as a multiplier of R200m")
	flag.Float64Var(&p.Sphere, "Sphere", 0.0,
		"Setting to a non-zero value, r, indicates that the particles " +
			"supplied to the profile should be treated as constant density " +
			"spheres of radius r instead of particles.")
	flag.IntVar(&p.SpherePts, "SpherePts", 0,
		"The number of points to use when using a spherical kernel.")
	flag.IntVar(&p.Tetra, "Tetra", 0,
		"Setting to a positive value, n, indicates that profiles should be " + 
			"generated from constant density tetrahedra instead of " +
			"particles. n will be the number of particles 'on a side' used " + 
			"during interpolation.")
	flag.IntVar(&p.Linear, "Linear", 0,
		"Setting to a positive value, n, indicates that profiles should be " + 
			"generated from tri-linear interpolation instead of from" +
			"particles. n will be the number of particles 'on a side' used " + 
			"during interpolation.")
	flag.IntVar(&p.Cubic, "Cubic", 0,
		"Setting to a positive value, n, indicates that profiles should be " + 
			"generated from tri-cubic interpolation instead of from" +
			"particles. n will be the number of particles 'on a side' used " + 
			"during interpolation.")
	flag.Parse()

	return p
}

// Determines how the density field is interpolated from the snapshot
// particle locations.
type Method int
const (
	Particle Method = iota
	Sphere
	Tetra
	Linear
	Cubic
)

// This is where the magic happens. Returns a slice of profiles. Each
// profile is represented as a slice of radii with units of comoving Mpc/h
// and densities with units of <rho_m>.
func computeGrids(ids, snaps []int, p *Params) ([]string, error) {
	// Select method.
	method := Particle
	switch {
	case p.Sphere > 0:
		method = Sphere
	case p.Tetra > 0:
		method = Tetra
	case p.Linear > 0:
		method = Linear
	case p.Cubic > 0:
		method = Cubic
	}

	// Set up loop state.
	buf := loop.NewBuffer()

	var (
		err error
		rs []*obj.Renderer
		srs []*obj.BallKernelRenderer
	)
	switch method {
	case Sphere:
		srs, err = newBallKernelRenderers(ids, snaps, p, method)
		rs = make([]*obj.Renderer, len(srs))
		for i := range srs { rs[i] = &srs[i].Renderer }
	default:
		rs, err = newRenderers(ids, snaps, p, method)
	}
	if err != nil { return nil, err }
	workers := runtime.GOMAXPROCS(runtime.NumCPU())

	snapBins, snapIdxs := binBySnap(snaps, ids)
	for snap := range snapBins {
		runtime.GC()

		// Select out all the halos in a given snapshot.
		idxs := snapIdxs[snap]
		snapRs := make([]loop.Object, len(idxs))
		switch method {
		case Sphere:
			for i, idx := range idxs { snapRs[i] = srs[idx] }
		default:
			for i, idx := range idxs { snapRs[i] = rs[idx] }
		}

		// Call loop.Loop().
		switch method {
		case Particle, Tetra, Sphere:
			loop.Loop(
				snap, snapRs, buf, p.SubsampleLength,
				loop.Linear, 1, workers,
			)
		case Linear:
			loop.Loop(
				snap, snapRs, buf, p.SubsampleLength,
				loop.Linear, p.Linear, workers,
			)
		case Cubic:
			loop.Loop(
				snap, snapRs, buf, p.SubsampleLength,
				loop.Cubic, p.Cubic, workers,
			)
		default:
			panic(":3")
		}
	}

	// Convert Profiles to the needed float slices.
	hd, err := util.ReadSnapHeader(snaps[0])
	dx := hd.TotalWidth / float64(hd.CountWidth) * float64(p.SubsampleLength)
	
	var pts int
	switch method {
	case Particle: pts = 1
	case Sphere: pts = p.SpherePts
	case Tetra: pts = p.Tetra
	case Linear: pts = p.Linear
	case Cubic: pts = p.Cubic
	default:
		panic(":3")
	}

	bufs := processCounts(rs, dx, pts)	
	fnames := fileNames(snaps, ids, p, method)
	err = writeBufs(snaps, bufs, rs, p, method, fnames)
	if err != nil { return nil, err }

	return fnames, nil
}

// binBySnap collects the given IDs into groups that share a common snapshot.
// The return format is two maps which accept snap
func binBySnap(snaps, ids []int) (snapBins, idxBins map[int][]int) {
    snapBins = make(map[int][]int)
    idxBins = make(map[int][]int)
    for i, snap := range snaps {
        id := ids[i]
        snapBins[snap] = append(snapBins[snap], id)
		idxBins[snap] = append(idxBins[snap], i)
    }
    return snapBins, idxBins
}

func newRenderers(
	ids, snaps []int, p *Params, method Method,
) ([]*obj.Renderer, error) {
	snapBins, idxBins := binBySnap(snaps, ids)
	rs := make([]*obj.Renderer, len(ids))

	for _, snap := range snaps {
		snapIDs := snapBins[snap]
        idxs := idxBins[snap]

        vals, err := util.ReadRockstar(
            snap, snapIDs, halo.X, halo.Y, halo.Z, halo.Rad200b,
        )
        if err != nil { return nil, err }

        xs, ys, zs, rads := vals[0], vals[1], vals[2], vals[3]

        for i := range xs {
			dr := rads[i] * p.Mult
			origin := [3]float64{ xs[i]-dr, ys[i]-dr, zs[i]-dr }
			pw := 2 * dr / float64(p.Pixels)
			pixels := [3]int{ p.Pixels, p.Pixels, p.Pixels }
			
			switch method {
			case Particle, Linear, Cubic:
				rs[idxs[i]] = obj.NewParticleRenderer(origin, pixels, pw)
			case Tetra:
				rs[idxs[i]] = obj.NewTetraRenderer(
					origin, pixels, pw, p.Tetra,
				)
			case Sphere:
				panic("NYI")
			}
        }
    }
    return rs, nil
}

func newBallKernelRenderers(
	ids, snaps []int, p *Params, method Method,
) ([]*obj.BallKernelRenderer, error) {
	snapBins, idxBins := binBySnap(snaps, ids)
	rs := make([]*obj.BallKernelRenderer, len(ids))

	for _, snap := range snaps {
		snapIDs := snapBins[snap]
        idxs := idxBins[snap]

        vals, err := util.ReadRockstar(
            snap, snapIDs, halo.X, halo.Y, halo.Z, halo.Rad200b,
        )
        if err != nil { return nil, err }

        xs, ys, zs, rads := vals[0], vals[1], vals[2], vals[3]

        for i := range xs {
			dr := rads[i] * p.Mult
			origin := [3]float64{ xs[i]-dr, ys[i]-dr, zs[i]-dr }
			pw := 2 * dr / float64(p.Pixels)
			pixels := [3]int{ p.Pixels, p.Pixels, p.Pixels }
			
			rs[idxs[i]] = obj.NewBallKernelRenderer(
				origin, pixels, pw, p.Sphere * rads[i],
				p.SpherePts*p.SpherePts*p.SpherePts,
			)
        }
    }
    return rs, nil
}

// processCounts takes raw count grids and converts them to density buffers.
func processCounts(rs []*obj.Renderer, dx float64, pts int) []density.Buffer {
	bufs:= make([]density.Buffer, len(rs))
	mp := dx*dx*dx / float64(pts*pts*pts)
	
	for i, r := range rs {
		rhos := make([]float64, len(r.Counts))
		dV := r.Pw*r.Pw*r.Pw
		for j, n := range r.Counts {
			rhos[j] = float64(n) * mp / dV
		}
		bufs[i] = density.WrapperDensityBuffer(rhos)
	}

	return bufs
}

func fileNames(snaps, ids []int, p *Params, method Method) []string {
	names := make([]string, len(ids))

	var paramStr string
	switch method {
	case Particle:
		paramStr = "ngp"
	case Sphere:
		paramStr = "sph"
	case Tetra:
		paramStr = fmt.Sprintf("tet%d", p.Tetra)
	case Linear:
		paramStr = fmt.Sprintf("lin%d", p.Linear)
	case Cubic:
		paramStr = fmt.Sprintf("cbc%d", p.Cubic)
	}

	for i := range snaps {
		names[i] = fmt.Sprintf("s%d_id%d_%s.gtet", snaps[i], ids[i], paramStr)
	}

	return names
}

func writeBufs(
	snaps []int, bufs []density.Buffer, rs []*obj.Renderer,
	p *Params, method Method, fnames []string,
) error {
	var pts int
	switch method {
	case Particle: pts =1
	case Sphere: pts = 1
	case Tetra: pts = p.Tetra
	case Linear: pts = p.Linear
	case Cubic: pts = p.Cubic
	}

	for i := range snaps {
		hd, err := util.ReadSnapHeader(snaps[i])
		if err != nil { return err }

		f, err := os.Create(fnames[i])
		defer f.Close()
		if err != nil { return err }

		// This isn't completely accurate, but that's fine.
		cellOrigin := [3]int{
			int(rs[i].Origin[0] / rs[i].Pw),
			int(rs[i].Origin[1] / rs[i].Pw),
			int(rs[i].Origin[2] / rs[i].Pw),
		}

		loc := io.NewLocationInfo(cellOrigin, rs[i].Pixels, rs[i].Pw)
		cos := io.NewCosmoInfo(
			hd.Cosmo.H100 * 100, hd.Cosmo.OmegaM,
			hd.Cosmo.OmegaL, hd.Cosmo.Z, hd.TotalWidth,
		)
		ren := io.NewRenderInfo(
			pts, int(rs[i].Pw / hd.TotalWidth), p.SubsampleLength, "",
		)

		io.WriteBuffer(bufs[i], cos, ren, loc, f)
	}

	return nil
}
