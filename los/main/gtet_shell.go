package main

import (
	"flag"
	"log"
	"math"
	"runtime"
	"sort"
	
	"github.com/phil-mansfield/gotetra/los"
	"github.com/phil-mansfield/gotetra/los/geom"
	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/los/analyze"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/math/rand"
)

// TODO: Someone needs to come in and restructure this monstrosity.

type Params struct {
	// HaloProfiles params
	RBins, Spokes, Rings int
	MaxMult, MinMult float64

	// Splashback params
	HFactor float64
	Order, Window, Levels, SubsampleLength int
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

	if len(ids) == 0 { log.Fatal("No input IDs.") }
	
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

	var losBuf *los.Buffers

	var rowLength int
	switch {
	case p.MedianProfile:
		rowLength = p.RBins * 2
	case p.MeanProfile:
		rowLength = p.RBins * 2
	default:
		rowLength = p.Order*p.Order*2
	}
	
	for i := range out {
		out[i] = make([]float64, rowLength)
	}

	
MainLoop:
	for _, snap := range sortedSnaps { 
		if snap == -1 { continue }
		snapIDs := snapBins[snap]
		idxs := idxBins[snap]

		// Bin halos
		hds, files, err := util.ReadHeaders(snap)
		
		if err != nil { log.Fatal(err.Error()) }
		if losBuf == nil {
			losBuf = los.NewBuffers(files[0], &hds[0], p.SubsampleLength)
		}
		halos, err := createHalos(snap, &hds[0], snapIDs, p)
		for i := range halos {
			// Screw it, we're too early in the catalog. Abort!
			if !halos[i].IsValid { continue MainLoop }
		}

		ms := runtime.MemStats{}
		runtime.ReadMemStats(&ms)
		log.Printf(
			"gtet_shell: Alloc: %d MB, Sys: %d MB",
			ms.Alloc / 1000000, ms.Sys / 1000000,
		)
		
		if err != nil { log.Fatal(err.Error()) }
		intrBins := binIntersections(hds, halos)
		
		// Add densities. Done header by header to limit I/O time.
		hdContainer := make([]io.SheetHeader, 1)
		fileContainer := make([]string, 1)
		for i := range hds {
			runtime.GC()

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
			}
		} else if p.MeanProfile {
			for i := range halos {
				runtime.GC()
				out[idxs[i]] = calcMean(&halos[i], p)
			}
		} else {
			// Calculate Penna coefficients.
			for i := range halos {
				runtime.GC()
				var ok bool
				out[idxs[i]], ok = calcCoeffs(&halos[i], buf, p)
				if !ok { log.Fatal("Welp, fix this.") }
			}
		}
		
	}
	
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
	flag.Float64Var(&p.HFactor, "HFactor", 10.0,
		"Factor controling how much an angular wedge can vary from " +
			"its neighbor. (If you don't know what this is, don't change it.)")
	flag.IntVar(&p.Order, "Order", 5,
		"Order of the shell fitting function.")
	flag.IntVar(&p.Window, "Window", 121,
		"Number of bins within smoothign window. Must be odd.")
	flag.IntVar(&p.Levels, "Levels", 4,
		"The number of recurve max-finding levels used by the 2D edge finder.")
	flag.IntVar(&p.SubsampleLength, "SubsampleLength", 1,
		"The number of particle edges per tetrahedron edge. Must be 2^n.")
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
	g := rand.NewTimeSeed(rand.Xorshift)

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
				los.Rotate(float32(g.Uniform(0, 2 * math.Pi)),
                    float32(g.Uniform(0, 2 * math.Pi)),
                    float32(g.Uniform(0, 2 * math.Pi))),
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

type profileRange struct {
	rMin, rMax float64
	v0 rgeom.Vec
}

func newProfileRanges(ids, snaps []int, p *Params) ([]profileRange, error) {
	snapBins, idxBins := binBySnap(snaps, ids)
	ranges := make([]profileRange, len(ids))
	for _, snap := range snaps {
		snapIDs := snapBins[snap]
		idxs := idxBins[snap]

		vals, err := util.ReadRockstar(
			snap, snapIDs, halo.X, halo.Y, halo.Z, halo.Rad200b,
		)
		if err != nil { return nil, err }
		xs, ys, zs, rs := vals[0], vals[1], vals[2], vals[3]
		for i := range xs {
			pr := profileRange{
				p.MinMult * rs[i], p.MaxMult * rs[i],
				rgeom.Vec{ float32(xs[i]), float32(ys[i]), float32(zs[i]) },
			}
			ranges[idxs[i]] = pr
		}
	}
	return ranges, nil
}

// Used for load balancing.
func zCounts(grid []bool, n int) []int {
	counts := make([]int, n)

	i := 0
	for z := 0; z < n; z++ {
		for y := 0; y < n; y++ {
			for x := 0; x < n; x++ {
				if grid[i] { counts[z]++ }
				i++
			}
		}
	}
	
	return counts
}

// Used for load balanacing.
func zSplit(zCounts []int, workers int) [][]int {
	tot := 0
	for _, n := range zCounts { tot += n }

	splits := make([]int, workers + 1)
	si := 1
	splitWidth := tot / workers
	if splitWidth * workers < tot { splitWidth++ }
	target := splitWidth

	sum := 0
	for i, n := range zCounts {
		sum += n
		if sum > target {
			splits[si] = i
			for sum > target { target += splitWidth }
			si++
		}
	}
	for ; si < len(splits); si++ { splits[si] = len(zCounts) }

	splitIdxs := make([][]int, workers)
	for i := range splitIdxs {
		jStart, jEnd := splits[i], splits[i + 1]
		for j := jStart; j < jEnd; j++ {
			if zCounts[j] > 0 { splitIdxs[i] = append(splitIdxs[i], j) }
		}
	}

	return splitIdxs
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

func binRangeIntersections(
	hds []io.SheetHeader, ranges []profileRange, idxs []int,
) [][]int {
	bins := make([][]int, len(hds))
	for i := range hds {
		for _, j := range idxs {
			if sheetIntersect(ranges[idxs[j]], &hds[i]) {
				bins[i] = append(bins[i], j)
			}
		}
	}
	return bins
}

func sheetIntersect(r profileRange, hd *io.SheetHeader) bool {
	tw := float32(hd.TotalWidth)
	return inRange(r.v0[0], float32(r.rMax), hd.Origin[0], hd.Width[0], tw) &&
		inRange(r.v0[1], float32(r.rMax), hd.Origin[1], hd.Width[1], tw) &&
		inRange(r.v0[2], float32(r.rMax), hd.Origin[2], hd.Width[2], tw)
}

func inRange(x, r, low, width, tw float32) bool {
	return wrapDist(x, low, tw) > -r && wrapDist(x, low + width, tw) < r
}

func wrapDist(x1, x2, width float32) float32 {
	dist := x1 - x2
	if dist > width / 2 {
		return dist - width
	} else if dist < width / -2 {
		return dist + width
	} else {
		return dist
	}
}

func prependRadii(rhos [][]float64, ranges []profileRange) [][]float64 {
	out := make([][]float64, len(rhos))
	for i := range rhos {
		rs := make([]float64, len(rhos[i]))
		lrMin, lrMax := math.Log(ranges[i].rMin), math.Log(ranges[i].rMax)
		dlr := (lrMax - lrMin) / float64(len(rhos[i]))
		for j := range rs {
			rs[j] = math.Exp((float64(j) + 0.5) * dlr + lrMin)
		}
		out[i] = append(rs, rhos[i]...)
	}
	return out
}

func calcCoeffs(
	halo *los.HaloProfiles, buf []analyze.RingBuffer, p *Params,
) ([]float64, bool) {
	for i := range buf {
		buf[i].Clear()
		buf[i].Splashback(halo, i, p.Window, p.Cutoff)
	}
	pxs, pys, ok := analyze.FilterPoints(buf, p.Levels, p.HFactor)
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
