package main

import (
	"fmt"
	"log"
	"flag"
	"math"
	"runtime"

	"github.com/phil-mansfield/gotetra/render/halo"
	obj "github.com/phil-mansfield/gotetra/los/main/loop_objects"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
	"github.com/phil-mansfield/gotetra/los/main/gtet_util/loop"
)

func main() {
	log.Println("gtet_prof")

	p := parseCmd()
	ids, snaps, _, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }
	if len(ids) == 0 { log.Fatal("No input IDs.") }

	out, err := genProfiles(ids, snaps, p)
	if err != nil { log.Fatal(err.Error()) }
	util.PrintRows(ids, snaps, out)
}

type Params struct {
	RBins, SubsampleLength int
	MaxMult, MinMult float64

	Sphere float64
	Tetra, Linear, Cubic int
}

// parseCmd parses the command line options and returns them as a struct.
func parseCmd() *Params {
	p := &Params{}

	flag.IntVar(&p.RBins, "RBins", 256,
		"Number of radial bins per profile.")
	flag.IntVar(&p.SubsampleLength, "SubsampleLength", 1,
		"Grid distance between adjacent points used in interpolation.")
	flag.Float64Var(&p.MaxMult, "MaxMult", 3,
		"Ending radius of the profiles as a multiplier of R_200m.")
	flag.Float64Var(&p.MinMult, "MinMult", 0.5,
		"Starting radius of the profiles as a multiplier of R_200m.")
	flag.Float64Var(&p.Sphere, "Sphere", 0.0,
		"Setting to a non-zero value, r, indicates that the particles " +
			"supplied to the profile should be treated as constant density " +
			"spheres of radius r instead of particles.")
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
func genProfiles(ids, snaps []int, p *Params) ([][]float64, error) {
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

	buf := loop.NewBuffer()
	profs, err := newProfiles(
		ids, snaps, p.MinMult, p.MaxMult, p.RBins, method,
	)
	if err != nil { return nil, err }

	workers := runtime.GOMAXPROCS(runtime.NumCPU())

	snapBins, snapIdxs := binBySnap(snaps, ids)
	for snap := range snapBins {
		runtime.GC()

		// Select out all the halos in a given snapshot.
		idxs := snapIdxs[snap]
		snapProfiles := make([]loop.Object, len(idxs))
		for i, idx := range idxs { snapProfiles[i] = profs[idx] }

		// Call loop.Loop().
		switch method {
		case Particle:
			log.Println("Calling loop.Loop().")
			loop.Loop(
				snap, snapProfiles, buf, p.SubsampleLength,
				loop.Linear, 1, workers,
			)
		default:
			panic("NYI")
		}
	}


	hd, err := util.ReadSnapHeader(snaps[0])
	dx := hd.TotalWidth / (float64(hd.CountWidth) * float64(p.SubsampleLength))
	
	var pts int
	switch method {
	case Particle, Sphere: pts = 1
	case Tetra: pts = p.Tetra
	case Linear: pts = p.Linear
	case Cubic: pts = p.Cubic
	default:
		panic(":3")
	}
	outs := processCounts(profs, dx, pts)
	
	return outs, nil
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

// newProfiles creates a slice of new empty profiles.
func newProfiles(
	ids, snaps []int, minMult, maxMult float64, rBins int, method Method,
) ([]*obj.Profile, error) {
	snapBins, idxBins := binBySnap(snaps, ids)
	profs := make([]*obj.Profile, len(ids))

	for _, snap := range snaps {
		snapIDs := snapBins[snap]
        idxs := idxBins[snap]

        vals, err := util.ReadRockstar(
            snap, snapIDs, halo.X, halo.Y, halo.Z, halo.Rad200b,
        )
        if err != nil { return nil, err }
        xs, ys, zs, rs := vals[0], vals[1], vals[2], vals[3]
        for i := range xs {
			switch method {
			case Particle:
				profs[idxs[i]] = obj.NewProfile(
					[3]float64{xs[i], ys[i], zs[i]},
					rs[i]*minMult, rs[i]*maxMult, rBins,
				)
			default:
				panic(fmt.Sprintf("Method %d not implemented.", method))
			}
        }
    }
    return profs, nil
}

// processCounts takes raw count profiles and converts them to properly
// formatted output profiles. dx is the distance between adjacent particles and
// pts in the number of points used "one a side" of each Lagrangian unit.
func processCounts(profs []*obj.Profile, dx float64, pts int) [][]float64 {
	outs := make([][]float64, len(profs))

	mp := dx*dx*dx / float64(pts*pts*pts)
	
	for i, prof := range profs {
		n := len(prof.Counts)
		out := make([]float64, 2*n)
		rs, rhos := out[0:n], out[n:2*n]

		dlr := (math.Log(prof.RMax) - math.Log(prof.RMin)) / float64(n)
		lrMin := math.Log(prof.RMin)
		
		for j := range rs {
			rs[j] = math.Exp(lrMin + dlr*(float64(j) + 0.5))

			rLo := math.Exp(dlr*float64(j) + lrMin)
			rHi := math.Exp(dlr*float64(j+1) + lrMin)
			dV := (rHi*rHi*rHi - rLo*rLo*rLo) * 4 * math.Pi / 3
			rhos[j] = prof.Counts[j] * mp / dV
		}
		
		outs[i] = out
	}

	return outs
}
