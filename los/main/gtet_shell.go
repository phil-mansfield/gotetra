
package main

import (
	"flag"
	"log"
	"math"
	"math/rand"
	"runtime"
	"sort"

	"github.com/phil-mansfield/gotetra/los"
	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/los/analyze"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/render/halo"
)

type Params struct {
	// HaloProfiles params
	RBins, Spokes, Rings int
	MaxMult, MinMult float64

	// Splashback params
	Order, Window, Levels int
	Cutoff float64

	// Alternate modes
	MedianProfile, MeanProfile bool
}

func main() {
	// Parse.
	log.Println("gtet_shell")
	p := parseCmd()
	ids, snaps, _, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }

	if len(ids) == 0 { return }

	// Compute coefficients.
	out := make([][]float64, len(ids))
	snapBins, idxBins := binBySnap(snaps, ids)
	buf := make([]analyze.RingBuffer, p.Rings)
	for i := range buf { buf[i].Init(p.Spokes, p.RBins) }

	sortedSnaps := []int{}
	for snap := range snapBins {
		sortedSnaps = append(sortedSnaps, snap)
	}
	sort.Ints(sortedSnaps)

	valids := make([]bool, len(ids))

	var losBuf *los.Buffers

MainLoop:
	for _, snap := range sortedSnaps { 
		log.Println("Snap", snap)
		if snap == -1 { continue }
		snapIDs := snapBins[snap]
		idxs := idxBins[snap]

		// Bin halos
		hds, files, err := util.ReadHeaders(snap)
		if err != nil { err.Error() }
		if losBuf == nil { losBuf = los.NewBuffers(files[0], &hds[0]) }
		halos, err := createHalos(snap, &hds[0], snapIDs, p)
		for i := range halos {
			// Screw it, we're too early in the catalog. Abort!
			if !halos[i].IsValid { continue MainLoop }
		}
		
		if err != nil { log.Fatal(err.Error()) }
		intrBins := binIntersections(hds, halos)

		// Add densities. Done header by header to limit I/O time.
		preMS, postMS := runtime.MemStats{}, runtime.MemStats{}
		hdContainer := make([]io.SheetHeader, 1)
		fileContainer := make([]string, 1)
		for i := range hds {
			log.Printf(
				"gtet_shell: analyzing sheet%d%d%d.dat", i/64, (i/8)%8, i%8,
			)
			runtime.ReadMemStats(&preMS)
			runtime.GC()
			runtime.ReadMemStats(&postMS)
			log.Printf(
				"gtet_shell: Pre GC: %d MB (%d MB), Post GC: %d MB (%d MB)",
				preMS.Alloc / 1000000, preMS.Sys / 1000000,
				postMS.Alloc / 1000000, postMS.Sys / 1000000,
			)
			if len(intrBins[i]) == 0 { continue }
			hdContainer[0] = hds[i]
			fileContainer[0] = files[i]
			los.LoadPtrDensities(
				intrBins[i], hdContainer, fileContainer, losBuf,
			)
		}
		
		if p.MedianProfile {
			// Calculate median profile.
			for i := range halos {
				runtime.GC()
				out[idxs[i]] = calcMedian(&halos[i], p)
				valids[idxs[i]] = true
			}
		} else if p.MeanProfile {
			for i := range halos {
				runtime.GC()
				out[idxs[i]] = calcMean(&halos[i], p)
				valids[idxs[i]] = true
			}
		} else {
			// Calculate Penna coefficients.
			for i := range halos {
				runtime.GC()
				var ok bool
				out[idxs[i]], ok = calcCoeffs(&halos[i], buf, p)
				if ok { valids[idxs[i]] = true }
			}
		}
		
	}
	
	ids = util.Filter(ids, valids)
	snaps = util.Filter(snaps, valids)
	
	util.PrintRows(ids, snaps, out)
}

func parseCmd() *Params {
	// Parse command line.
	p := &Params{}
	flag.IntVar(&p.RBins, "RBins", 256,
		"Number of radial bins used per LoS.")
	flag.IntVar(&p.Spokes, "Spokes", 1024,
		"Number of LoS's used per ring.")
	flag.IntVar(&p.Rings, "Rings", 10,
		"Number of rings used per halo. 3, 4, 6, and 10 rings are\n" + 
			"guaranteed to be uniformly spaced.")
	flag.Float64Var(&p.MaxMult, "MaxMult", 3,
		"Ending radius of LoSs as a multiple of R_200m.")
	flag.Float64Var(&p.MinMult, "MinMult", 0.5,
		"Starting radius of LoSs as a multiple of R_200m.")
	flag.IntVar(&p.Order, "Order", 5,
		"Order of the shell fitting function.")
	flag.IntVar(&p.Window, "Window", 121,
		"Number of bins within smoothign window. Must be odd.")
	flag.IntVar(&p.Levels, "Levels", 4,
		"The number of recurve max-finding levels used by the 2D edge finder.")
	flag.Float64Var(&p.Cutoff, "Cutoff", 0.0,
		"The shallowest slope that can be considered a splashback point.")
	flag.BoolVar(&p.MedianProfile, "MedianProfile", false,
		"Compute the median halo profile instead of the shell. " + 
			"KILL THIS OPTION.")
		flag.BoolVar(&p.MeanProfile, "MeanProfile", false,
		"Compute the mean halo profile instead of the shell. " + 
			"KILL THIS OPTION.")
	flag.Parse()
	return p
}

func createHalos(
	snap int, hd *io.SheetHeader, ids []int, p *Params,
) ([]los.HaloProfiles, error) {
	vals, err := util.ReadRockstar(
		snap, ids, halo.X, halo.Y, halo.Z, halo.Rad200b,
	)
	if err != nil { return nil, err }

	xs, ys, zs, rs := vals[0], vals[1], vals[2], vals[3]
	
	// Initialize halos.
	halos := make([]los.HaloProfiles, len(ids))
	seenIDs := make(map[int]bool)
	for i, id := range ids {
		origin := &geom.Vec{
			float32(xs[i]), float32(ys[i]), float32(zs[i]),
		}

		if rs[i] <= 0 { continue }
		
		// If we've already seen a halo once, randomize its orientation.
		if seenIDs[id] {
			halos[i].Init(
				id, p.Rings, origin, rs[i] * p.MinMult, rs[i] * p.MaxMult,
				p.RBins, p.Spokes, hd.TotalWidth, los.Log(true),
				los.Rotate(float32(2 * math.Pi * rand.Float64()),
                    float32(2 * math.Pi * rand.Float64()),
                    float32(2 * math.Pi * rand.Float64())),
			)
		} else {
			seenIDs[id] = true
			halos[i].Init(
				id, p.Rings, origin, rs[i] * p.MinMult, rs[i] * p.MaxMult,
				p.RBins, p.Spokes, hd.TotalWidth, los.Log(true),
			)
		}
	}

	return halos, nil
}

func binIntersections(
	hds []io.SheetHeader, halos []los.HaloProfiles,
) [][]*los.HaloProfiles {

	bins := make([][]*los.HaloProfiles, len(hds))
	for i := range hds {
		for hi := range halos {
			if (&halos[hi]).SheetIntersect(&hds[i]) {
				bins[i] = append(bins[i], &halos[hi])
			}
		}
	}
	return bins
}

func calcCoeffs(
	halo *los.HaloProfiles, buf []analyze.RingBuffer, p *Params,
) ([]float64, bool) {
	for i := range buf {
		buf[i].Clear()
		buf[i].Splashback(halo, i, p.Window, p.Cutoff)
	}
	pxs, pys, ok := analyze.FilterPoints(buf, p.Levels)
	if !ok { return nil, false }
	cs, _ := analyze.PennaVolumeFit(pxs, pys, halo, p.Order, p.Order)
	return cs, true
}

func calcMedian(halo *los.HaloProfiles, p *Params) []float64 {
	rs := make([]float64, p.RBins)
	halo.GetRs(rs)
	rhos := halo.MedianProfile()
	return append(rs, rhos...)
}

func calcMean(halo *los.HaloProfiles, p *Params) []float64 {
	rs := make([]float64, p.RBins)
	halo.GetRs(rs)
	rhos := halo.MeanProfile()
	return append(rs, rhos...)
}

func intr(xs []float64) []interface{} {
	is := make([]interface{}, len(xs))
	for i, x := range xs { is[i] = x }
	return is
}

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
