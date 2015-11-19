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

type Params struct {
	// HaloProfiles params
	RBins, Spokes, Rings int
	MaxMult, MinMult float64

	// Splashback params
	HFactor float64
	Order, Window, Levels, SubsampleLength int
	Cutoff float64

	// Alternate modes
	MedianProfile, MeanProfile, SphericalProfile bool
	SphericalProfilePoints int
}

func main() {
	// Parse.
	log.Println("gtet_shell")
	p := parseCmd()
	ids, snaps, _, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }

	if len(ids) == 0 { log.Fatal("No input IDs.") }

	// We're just going to do this part separately.
	if p.SphericalProfile {
		out, err := sphericalProfile(ids, snaps, p)
		if err != nil { log.Fatal(err.Error()) }
		util.PrintRows(ids, snaps, out)
		return
	}
	
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
	flag.BoolVar(&p.SphericalProfile, "SphericalProfile", false,
		"Compute the radial profile of a halo using standard particle " +
			"binning. KILL THIS OPTION.")
	flag.IntVar(&p.SphericalProfilePoints, "SphericalProfilePoints", 0,
		"Number of pointers per tetrhedra to use when computing spherical " +
			"If 0, tetrahedra won't be used and the profiles will just be" +
			"computed from the particles.")
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

func sphericalProfile(ids, snaps []int, p *Params) ([][]float64, error) {
	// Normal Set up
	var xs []rgeom.Vec
	counts := make([][]float64, len(ids))
	for i := range counts { counts[i] = make([]float64, p.RBins) }
	ranges, err := newProfileRanges(ids, snaps, p)
	if err != nil { return nil, err }
	
	// Tetra setup.
	var (
		vecBuf []rgeom.Vec
		randBuf []float64
		gen *rand.Generator
	)
	if p.SphericalProfilePoints > 0 {
		gen = rand.NewTimeSeed(rand.Xorshift)
		vecBuf = make([]rgeom.Vec, p.SphericalProfilePoints)
		randBuf = make([]float64, 3 * p.SphericalProfilePoints)
	}


	// Bin particles
	snapBins, idxBins := binBySnap(snaps, ids)
	for snap := range snapBins {
		runtime.GC()
		
		idxs := idxBins[snap]
		hds, files, err := util.ReadHeaders(snap)
		if err != nil { return nil, err }
		if len(xs) == 0 {
			n := hds[0].GridWidth*hds[0].GridWidth*hds[0].GridWidth
			xs = make([]rgeom.Vec, n)
		}
		
		intrBins := binRangeIntersections(hds, ranges, idxs)
		for i := range hds {
			if len(intrBins[i]) == 0 { continue }
			log.Printf("%d%d%d -> (%d)", i / 64, (i / 8) % 8, i % 8,
				len(intrBins[i]))
			err := io.ReadSheetPositionsAt(files[i], xs)
			if err != nil { return nil, err }
			for _, j := range intrBins[i] {
				r := ranges[idxs[j]]
				if p.SphericalProfilePoints > 0 {
					tetraBinParticles(
						&hds[i], xs, counts[idxs[j]], p.SubsampleLength,
						r.rMin, r.rMax, r.v0, vecBuf, randBuf, gen,
						
					)
				} else {
					binParticles(
						&hds[i], xs, counts[idxs[j]], p.SubsampleLength,
						r.rMin, r.rMax, r.v0,
					)
				}
			}
		}
	}
	
	// Convert
	for i := range counts {
		rMin, rMax := ranges[i].rMin, ranges[i].rMax
		hds, _, err := util.ReadHeaders(snaps[i])
		if err != nil { return nil, err }
		countsToRhos(&hds[0], counts[i], rMin, rMax, p.SphericalProfilePoints)
	}

	outs := prependRadii(counts, ranges)
	
	return outs, nil
}

func countsToRhos(
	hd *io.SheetHeader, counts []float64, rMin, rMax float64, tetraPoints int,
) {
	dx := hd.TotalWidth / float64(hd.CountWidth)
	mp := dx*dx*dx
	if tetraPoints > 0 {
		mp /= float64(6 * tetraPoints)
	}

	lrMin, lrMax := math.Log(rMin), math.Log(rMax)
	dlr := (lrMax - lrMin) / float64(len(counts))
	for i := range counts {
		rLo := math.Exp(dlr*float64(i) + lrMin)
		rHi := math.Exp(dlr*float64(i + 1) + lrMin)
		dV := (rHi*rHi*rHi - rLo*rLo*rLo) * 4 * math.Pi / 3
		counts[i] *= mp / dV
	}
}

func tetraBinParticles(
	hd *io.SheetHeader, xs []rgeom.Vec, counts []float64, skip int,
	rMin, rMax float64, v0 rgeom.Vec,
	vecBuf []rgeom.Vec, randBuf []float64, gen *rand.Generator,
) {
	min2, max2 := rMin*rMin, rMax*rMax
	x0, y0, z0 := float64(v0[0]), float64(v0[1]), float64(v0[2])

	lrMin, lrMax := math.Log(rMin), math.Log(rMax)
	dlr := (lrMax - lrMin) / float64(len(counts))
	tw := hd.TotalWidth

	sw, gw := int(hd.SegmentWidth), int(hd.GridWidth)
	for iz := 0; iz < sw; iz += skip {
		for iy := 0; iy < sw; iy += skip {
			for ix := 0; ix < sw; ix += skip {
				idx := ix + gw*iy + gw*gw*iz
				for dir := 0; dir < 6; dir++ {
					tetraPoints(idx, dir, gw, skip, xs, gen, randBuf, vecBuf)
					for _, pt := range vecBuf {
						x := float64(pt[0]) - x0
						y := float64(pt[1]) - y0
						z := float64(pt[2]) - z0

						if x > tw / 2 { x -= tw }
						if y > tw / 2 { y -= tw }
						if z > tw / 2 { z -= tw }
						if x < -tw / 2 { x += tw }
						if y < -tw / 2 { y += tw }
						if z < -tw / 2 { z += tw }
						
						r2 := x*x + y*y + z*z
						if r2 <= min2 || r2 >= max2 { continue }
						
						lr := math.Log(r2) / 2
						ri := int((lr - lrMin) / dlr)
					
						counts[ri]++
					}
				}
			}
		}
	}
}

func tetraPoints(
	idx, dir, gw, skip int, xs []rgeom.Vec,
	gen *rand.Generator, randBuf []float64, vecBuf []rgeom.Vec,
) {
	idxBuf, tet := &rgeom.TetraIdxs{}, &rgeom.Tetra{}
	idxBuf.Init(int64(idx), int64(gw), int64(skip), dir)
	i0, i1, i2, i3 := idxBuf[0], idxBuf[1], idxBuf[2], idxBuf[3]
	tet.Init(&xs[i0], &xs[i1], &xs[i2], &xs[i3])
	tet.RandomSample(gen, randBuf, vecBuf)
}

func binParticles(
	hd *io.SheetHeader, xs []rgeom.Vec, counts []float64, skip int,
	rMin, rMax float64, v0 rgeom.Vec,
) {
	min2, max2 := rMin*rMin, rMax*rMax
	x0, y0, z0 := float64(v0[0]), float64(v0[1]), float64(v0[2])

	lrMin, lrMax := math.Log(rMin), math.Log(rMax)
	dlr := (lrMax - lrMin) / float64(len(counts))
	incr := float64(skip*skip*skip)
	tw := hd.TotalWidth

	sw, gw := int(hd.SegmentWidth), int(hd.GridWidth)
	for iz := 0; iz < sw; iz += skip {
		for iy := 0; iy < sw; iy += skip {
			for ix := 0; ix < sw; ix += skip {
				pt := xs[ix + iy*gw + iz*gw*gw]
				x := float64(pt[0]) - x0
				y := float64(pt[1]) - y0
				z := float64(pt[2]) - z0
				
				if x > tw / 2 { x -= tw }
				if y > tw / 2 { y -= tw }
				if z > tw / 2 { z -= tw }
				if x < -tw / 2 { x += tw }
				if y < -tw / 2 { y += tw }
				if z < -tw / 2 { z += tw }
				
				r2 := x*x + y*y + z*z
				if r2 <= min2 || r2 >= max2 { continue }
				
				lr := math.Log(r2) / 2
				ri := int((lr - lrMin) / dlr)
				
				counts[ri] += incr
			}
		}
	}
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
