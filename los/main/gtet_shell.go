package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"runtime"
	"sort"
	"strings"
	
	"github.com/phil-mansfield/gotetra/los"
	sph "github.com/phil-mansfield/gotetra/los/sphere_halo"
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
	// Halo profile params
	RBins, Spokes, Rings int
	MaxMult, MinMult float64
	// Halo profile params that are only used if set in Sphere mode.
	SphereMult, DefaultRhoMult float64

	// Splashback params
	HFactor float64
	Order, Window, Levels, SubsampleLength int
	Cutoff float64

	// Interpolation
	Interpolation string

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

	err = loop(ids, snaps, p, out)
	if err != nil { log.Fatal(err.Error()) }

	util.PrintRows(ids, snaps, out)
}

func loop(ids, snaps []int, p *Params, out [][]float64) error {
	snapBins, idxBins := binBySnap(snaps, ids)
	ringBuf := make([]analyze.RingBuffer, p.Rings)
	for i := range ringBuf { ringBuf[i].Init(p.Spokes, p.RBins) }

	sortedSnaps := []int{}
	for snap := range snapBins {
		sortedSnaps = append(sortedSnaps, snap)
	}
	sort.Ints(sortedSnaps)

	hds, files, err := util.ReadHeaders(sortedSnaps[0])
	if err != nil { return err }

	// Create buffers for the different methods.
	var (
		losBuf *los.Buffers
		sphBuf *sphBuffers
	)
	switch strings.ToLower(p.Interpolation) {
	case "tetra":
		losBuf = los.NewBuffers(files[0], &hds[0], p.SubsampleLength)
	case "sphere":
		workers := runtime.NumCPU()
		gw := hds[0].GridWidth
		sphBuf = &sphBuffers{
			intr: make([]bool, gw*gw*gw),
			vecs: make([]rgeom.Vec, gw*gw*gw),
			sphWorkers: make([]sph.SphereHalo, workers - 1),
		}
	default:
		panic(fmt.Sprintf("Unknown interpolation mode '%s'", p.Interpolation))
	}

	for _, snap := range sortedSnaps { 
		if snap == -1 { continue }
		snapIDs := snapBins[snap]
		idxs := idxBins[snap]

		// Create Halos
		runtime.GC()
		halos, err := createHalos(snap, &hds[0], snapIDs, p)
		if err != nil { return err }
		printMemStats()

		// Loop over objects.
		switch strings.ToLower(p.Interpolation) {
		case "tetra":
			err := tetraLoop(snap, ids, idxs, halos, p, losBuf, out)
			if err != nil { return err }
		case "sphere":
			err := sphereLoop(snap, ids, idxs, halos, p, sphBuf, out)
			if err != nil { return err }
		default:
			panic(fmt.Sprintf("Unknown interpolation mode '%s'",
				p.Interpolation))
		}

		// Analysis
		haloAnalysis(halos, idxs, p, ringBuf, out)
	}

	return nil
}

func sphereLoop(
	snap int, IDs, ids []int, halos []los.Halo, p *Params,
	sphBuf *sphBuffers, out [][]float64,
) error {
	hds, files, err := util.ReadHeaders(snap)
	if err != nil { return err }
	intrBins := binIntersections(hds, halos)

	for i := range hds {
		runtime.GC()
		if len(intrBins[i]) == 0 { continue }
		
		err := io.ReadSheetPositionsAt(files[i], sphBuf.vecs)
		if err != nil { return err }

		binHs := intrBins[i]
		for j := range binHs {
			h, ok := binHs[j].(*sph.SphereHalo)
			if !ok { panic("Invalid Halo interface given to sphereLoop().") }
			loadSphereVecs(h, sphBuf, &hds[i], p)
		}
		
	}

	return nil
}

type sphBuffers struct {
	sphWorkers []sph.SphereHalo
	vecs []rgeom.Vec
	intr []bool
}

func loadSphereVecs(
	h *sph.SphereHalo, sphBuf *sphBuffers, hd *io.SheetHeader, p *Params,
) {
	workers := runtime.NumCPU()
	runtime.GOMAXPROCS(workers)
	sphWorkers, vecs, intr := sphBuf.sphWorkers, sphBuf.vecs, sphBuf.intr
	if len(sphWorkers) + 1 != workers { panic("impossible")}

	sync := make(chan bool, workers)

	h.Transform(vecs, hd.TotalWidth)
	rad := h.RMax() * p.SphereMult / p.MaxMult
	h.Intersect(vecs, rad, intr)

	counts := zCounts(sphBuf.intr, int(hd.GridWidth))
	zIdxs := zSplit(counts, workers)

	h.Split(sphWorkers)

	for i := range sphWorkers {
		wh := &sphBuf.sphWorkers[i]
		go chanLoadSphereVec(wh, vecs, intr, zIdxs[i], hd, p, sync)
	}
	chanLoadSphereVec(h, vecs, intr, zIdxs[workers-1], hd, p, sync)

	for i := 0; i < workers; i++ { <-sync }

	h.Join(sphWorkers)
}

func chanLoadSphereVec(
	h *sph.SphereHalo, vecs []rgeom.Vec, intr []bool,
	zIdxs []int, hd *io.SheetHeader,  p *Params, sync chan bool,
) {
	gw, sw := int(hd.GridWidth), int(hd.SegmentWidth)

	rad := h.RMax() * p.SphereMult / p.MaxMult
	sphVol := 4*math.Pi/3*rad*rad*rad

	pl := hd.TotalWidth / float64(int(hd.CountWidth) / p.SubsampleLength)
	pVol := pl*pl*pl

	rho := pVol / sphVol
	for _, z := range zIdxs {
		if z >= sw { continue }

		for y := 0; y < sw; y++ {
			for x := 0; x < sw; x++ {

				idx := x + y*gw + z*gw*gw
				if intr[idx] { h.Insert(geom.Vec(vecs[idx]), rad, rho) }
			}
		}
	}

	sync <- true
}

// loopSnap loops over all the halos in a given snapshot.
func tetraLoop(
	snap int, IDs, idxs []int, halos []los.Halo, p *Params,
	losBuf *los.Buffers, out [][]float64,
) error {
	hds, files, err := util.ReadHeaders(snap)
	if err != nil { return err }

	intrBins := binIntersections(hds, halos)
	
	// Add densities. Done header by header to limit I/O time.
	hdContainer := make([]io.SheetHeader, 1)
	fileContainer := make([]string, 1)
	for i := range hds {
		runtime.GC()
		
		if len(intrBins[i]) == 0 { continue }
		hdContainer[0] = hds[i]
		fileContainer[0] = files[i]

		hps := make([]*los.HaloProfiles, len(intrBins[i]))
		for j := range hps {
			var ok bool
			hps[j], ok = intrBins[i][j].(*los.HaloProfiles)
			if !ok { panic("Invalid Halo given to tetraLoop().") }
		}

		los.LoadPtrDensities(hps, hdContainer, fileContainer, losBuf)
	}

	return nil
}

func haloAnalysis(
	halos []los.Halo, idxs []int, p *Params,
	ringBuf []analyze.RingBuffer, out [][]float64,
) {
	switch {
	case p.MedianProfile:
		// Calculate median profile.
		for i := range halos {
			runtime.GC()
			out[idxs[i]] = calcMedian(halos[i], p)
		}
	case p.MeanProfile:
		// Calculate mean profile.
		for i := range halos {
			runtime.GC()
			out[idxs[i]] = calcMean(halos[i], p)
		}
	default:
		// Calculate Penna coefficients.
		for i := range halos {
			runtime.GC()
			var ok bool
			out[idxs[i]], ok = calcCoeffs(halos[i], ringBuf, p)
			if !ok { log.Fatal("Welp, fix this.") }
		}
	}
}

// printMemStats prints out allocation statistics to the log files.
func printMemStats() {
	ms := runtime.MemStats{}
	runtime.ReadMemStats(&ms)
	log.Printf(
		"gtet_shell: Alloc: %d MB, Sys: %d MB",
		ms.Alloc / 1000000, ms.Sys / 1000000,
	)
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
		"Number of bins within smoothing window. Must be odd.")
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
	flag.StringVar(&p.Interpolation, "Interpolation", "Tetra",
		"Interpolation Mode. Must be one of [Tetra | Sphere]. Will " + 
			"eventaully support Cubic and Linear. Case insensitive.")
	flag.Float64Var(&p.SphereMult, "SphereMult", 0.1,
		"The radius of the spherical kernels as a multiplier of the halo's " +
			"R200m. Only valid if Interpolation flag is set to Sphere.")
	flag.Float64Var(&p.DefaultRhoMult, "DefaultRhoMult", 0.5,
		"The density assigned to zero-density cells as a multiplier " + 
			"of the dnesity of a single spherical kernel.")
	flag.Parse()
	return p
}

func createHalos(
	snap int, hd *io.SheetHeader, ids []int, p *Params,
) ([]los.Halo, error) {
	vals, err := util.ReadRockstar(
		snap, ids, halo.X, halo.Y, halo.Z, halo.Rad200b,
	)
	if err != nil { return nil, err }

	xs, ys, zs, rs := vals[0], vals[1], vals[2], vals[3]
	g := rand.New(rand.Xorshift, 0)

	// Initialize halos.
	halos := make([]los.Halo, len(ids))
	seenIDs := make(map[int]bool)
	
	for i, id := range ids {
		if rs[i] <= 0 { continue }
		
		switch strings.ToLower(p.Interpolation) {
		case "tetra":
			origin := &geom.Vec{
				float32(xs[i]), float32(ys[i]), float32(zs[i]),
			}
			
			// If we've already seen a halo once, randomize its orientation.
			halo := &los.HaloProfiles{}
			if seenIDs[id] {
				halo.Init(
					id, p.Rings, origin, rs[i] * p.MinMult, rs[i] * p.MaxMult,
					p.RBins, p.Spokes, hd.TotalWidth, los.Log(true),
					los.Rotate(float32(g.Uniform(0, 2 * math.Pi)),
						float32(g.Uniform(0, 2 * math.Pi)),
						float32(g.Uniform(0, 2 * math.Pi))),
				)
				halos[i] = halo
			} else {
				seenIDs[id] = true
				halo.Init(
					id, p.Rings, origin, rs[i] * p.MinMult, rs[i] * p.MaxMult,
					p.RBins, p.Spokes, hd.TotalWidth, los.Log(true),
				)
				halos[i] = halo
			}
			
		case "sphere":
			norms := normVecs(p.Rings)
			origin := [3]float64{
				float64(xs[i]), float64(ys[i]), float64(zs[i]),
			}
			rMax, rMin := rs[i] * p.MaxMult, rs[i] * p.MinMult
			rad := rs[i]*p.SphereMult
			sphVol := 4*math.Pi/3*rad*rad*rad
			pl := hd.TotalWidth/float64(int(hd.CountWidth)/p.SubsampleLength)
			pVol := pl*pl*pl
			rho := pVol / sphVol
			defaultRho := rho * p.DefaultRhoMult
			
			halo := &sph.SphereHalo{}
			halo.Init(norms, origin, rMin, rMax, p.RBins, p.Spokes, defaultRho)

			halos[i] = halo
		default:
			panic(fmt.Sprintf("Unknown Interpolation mode '%s'.",
				p.Interpolation))
		}
	}

	return halos, nil
}

func normVecs(n int) []geom.Vec {
	var vecs []geom.Vec
	gen := rand.New(rand.Xorshift, 0)
	switch n {
	case 3:
		vecs = []geom.Vec{{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}
	default:
		vecs = make([]geom.Vec, n)
		for i := range vecs {
			for {
				x := gen.Uniform(-1, +1)
				y := gen.Uniform(-1, +1)
				z := gen.Uniform(-1, +1)
				r := math.Sqrt(x*x + y*y + z*z)

				if r < 1 {
					vecs[i] = geom.Vec{
						float32(x/r), float32(y/r), float32(z/r),
					}
					break
				}
			}
		}
	}
	
	return vecs
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
	hds []io.SheetHeader, halos []los.Halo,
) [][]los.Halo {

	bins := make([][]los.Halo, len(hds))
	for i := range hds {
		for hi := range halos {
			if halos[hi].SheetIntersect(&hds[i]) {
				bins[i] = append(bins[i], halos[hi])
			}
		}
	}
	return bins
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
	halo los.Halo, buf []analyze.RingBuffer, p *Params,
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

func calcMedian(halo los.Halo, p *Params) []float64 {
	rs := make([]float64, p.RBins)
	halo.GetRs(rs)
	rhos := halo.MedianProfile()
	return append(rs, rhos...)
}

func calcMean(halo los.Halo, p *Params) []float64 {
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
