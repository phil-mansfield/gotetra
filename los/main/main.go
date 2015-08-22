package main

import (
	"encoding/binary"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"math/rand"
	"os"
	"path"
	"strings"
	"time"
	"runtime/pprof"

	"github.com/phil-mansfield/gotetra/render/io"	
	"github.com/phil-mansfield/gotetra/render/halo"

	"github.com/phil-mansfield/gotetra/los"
	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/los/analyze"

	plt "github.com/phil-mansfield/pyplot"
)

const (
	rType = halo.R200m
	rMaxMult = 3.0
	rMinMult = 0.5

	n = 1024
	bins = 256

	rings = 3
	plotStart = 0
	plotCount = 3

	// SubhaloFinder params.
	finderCells = 150
	overlapMult = 3

	hdSaveFile = "hdSave.dat"
)

var (
	colors = []string{
		"DarkSlateBlue", "DarkSlateGray", "DarkTurquoise",
		"DarkViolet", "DeepPink", "DimGray",
	}
	visProfs = []int{
		rand.Intn(n), rand.Intn(n), rand.Intn(n),
		rand.Intn(n), rand.Intn(n), rand.Intn(n),
	}
)

func loadHeaders(files []string, saveDir string) ([]io.SheetHeader, error) {
	saveFile := path.Join(saveDir, hdSaveFile)
	hds := make([]io.SheetHeader, len(files))

	if _, err := os.Stat(saveFile); err == nil {
		f, err := os.Open(saveFile)
		if err != nil { return nil, err }
		defer f.Close()
		fmt.Println("Loading saved headers.")
		binary.Read(f, binary.LittleEndian, hds)
		
	} else {
		fmt.Print("Loading individual headers: ")

		for i := range files {
			if i % 50 == 0 { fmt.Print(i, " ") }
			err = io.ReadSheetHeaderAt(files[i], &hds[i])
			if err != nil { return nil, err }
		}

		f, err := os.Create(saveFile)
		if err != nil { return nil, err }
		defer f.Close()
		binary.Write(f, binary.LittleEndian, hds)
		fmt.Println()

	}
	return hds, nil
}

func main() {
	// Argument Parsing.
	fmt.Println("Running")
	if len(os.Args) != 5 {
		log.Fatalf("Usage: $ %s input_dir halo_file plot_dir save_dir",
			os.Args[0])
	}

	dirName := os.Args[1]
	haloFileName := os.Args[2]
	plotDir := os.Args[3]
	saveDir := os.Args[4]

	// Do I/O and set up buffers.
	files, err := fileNames(dirName)
	if err != nil { log.Fatal(err.Error()) }
	hds, err := loadHeaders(files, saveDir)
	if err != nil { log.Fatal(err.Error()) }
	buf := los.NewBuffers(files[0], &hds[0])
	h := new(los.HaloProfiles)
	
	// Find halos, subhalos, etc.
	xs, ys, zs, ms, rs, err := halo.ReadRockstar(
		haloFileName, rType, &hds[0].Cosmo,
	)

	if err != nil { log.Fatal(err.Error()) }
	fmt.Printf("%d halos read.\n", len(xs))
	g := halo.NewGrid(finderCells, hds[0].TotalWidth, len(xs))
	g.Insert(xs, ys, zs)
	sf := halo.NewSubhaloFinder(g)
	sf.FindSubhalos(xs, ys, zs, rs, overlapMult)

	// Profiling boilerplate.
	f, err := os.Create("out.pprof")
	if err != nil { log.Fatal(err.Error()) }
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()

	// Analyze each halo.
	plotRs, plotRhos := make([]float64, bins), make([]float64, bins)
	_, _ = plotRs, plotRhos
	r := new(analyze.RingBuffer)
	r.Init(n, bins)
	for i := plotStart; i < plotStart + plotCount; i++ {
		fmt.Println("Hosts:", sf.HostCount(i), "Subhalos:", sf.SubhaloCount(i))
		spheres := subhaloSpheres(sf, i, xs, ys, zs, rs)
		
		origin := &geom.Vec{float32(xs[i]), float32(ys[i]), float32(zs[i])}
		h.Init(i, rings, origin, rs[i] * rMinMult, rs[i] * rMaxMult,
			bins, n, hds[0].TotalWidth, los.Log(false))
		hdIntrs, fileIntrs := intersectingSheets(h, hds, files)
		
		fmt.Printf(
			"Halo mass is: %.3g, intersects are: %d\n", ms[i], len(hdIntrs),
		)
			
		intersectionTest(h, hdIntrs, fileIntrs, buf, spheres)
		//plotExampleProfiles(h, ms[i], plotRs, plotRhos, plotDir, spheres)
		//plotExampleDerivs(h, ms[i], plotRs, plotRhos, plotDir, spheres)
		for ring := 0; ring < rings; ring++ {
			r.Clear()
			r.Splashback(h, ring, 61, -5)
			plotPlane(r, ms[i], h.ID(), ring, plotDir)
			plotExampleProfiles(
				h, ms[i], ring, plotRs, plotRhos, plotDir,
			)
		}
	}
	
	plt.Execute()
}

func subhaloSpheres(
	sf *halo.SubhaloFinder, i int, xs, ys, zs, rs []float64,
) []geom.Sphere {
	shIdxs := sf.Subhalos(i)
	subhalos := make([]geom.Sphere, len(shIdxs))

	for j, idx := range shIdxs {
		subhalos[j].R = float32(rs[idx])
		subhalos[j].C = geom.Vec{
			float32(xs[idx]), float32(ys[idx]), float32(zs[idx]),
		}
	}
	return subhalos
}

func plotExampleProfiles(
	hp *los.HaloProfiles, m float64, ring int,
	rs, rhos []float64, dir string,
) {
	fname := path.Join(dir, fmt.Sprintf("profs_h%d_r%d.png", hp.ID(), ring))

	plt.Figure()
	hp.GetRs(rs)

	r := rs[len(rs) - 1] / rMaxMult
	plt.Plot([]float64{r, r}, []float64{1e-2, 1e3}, "k", plt.LW(2))


	for _, visIdx := range visProfs {
		hp.GetRhos(ring, visIdx, rhos)
		rhoSets, auxSets := analyze.NaNSplit(rhos, analyze.Aux(rs))
		
		for i := range rhoSets {
			rawRs, rawRhos := auxSets[0][i], rhoSets[i]
			smoothRhos, smoothDerivs, ok := analyze.Smooth(rawRs, rawRhos, 61)
			if !ok { continue }
			plt.Plot(rawRs, smoothRhos, plt.LW(3), plt.C(colors[visIdx]))
			r, ok := analyze.SplashbackRadius(rawRs, smoothRhos, smoothDerivs)
			if !ok { continue }
			plt.Plot([]float64{r, r}, []float64{1e3,0.01}, plt.C(colors[visIdx]))
		}
	}

	// Plot specifications.
	plt.Title(fmt.Sprintf(
		`Halo %d: $M_{\rm 200c}$ = %.3g $M_\odot/h$`, hp.ID(), m),
	)
	plt.XLabel(`$R$ $[{\rm Mpc}/h]$`, plt.FontSize(16))
	plt.YLabel(`$\rho$ [$\rho_m$]`, plt.FontSize(16))

	plt.XScale("log")

	plt.YScale("log")
	plt.YLim(1e-2, 1e3)
	setXRange(rs[0], rs[len(rs) - 1])

	plt.Grid(plt.Axis("y"))
	plt.Grid(plt.Axis("x"), plt.Which("both"))
	plt.SaveFig(fname)
}

func plotExampleDerivs(
	hp *los.HaloProfiles, m float64, ring int,
	rs, rhos []float64, dir string,
) {
	fname := path.Join(dir, fmt.Sprintf("derivs_%dh.png", hp.ID()))

	plt.Figure()
	hp.GetRs(rs)

	r := rs[len(rs) - 1] / rMaxMult
	plt.Plot([]float64{r, r}, []float64{-20, +10}, "k", plt.LW(2))

	hp.GetRhos(ring, 13, rhos)
	rhoSets, auxSets := analyze.NaNSplit(rhos, analyze.Aux(rs))

	for _, visIdx := range visProfs {
		for i := range rhoSets {
			rawRs, rawRhos := auxSets[0][i], rhoSets[i]
			smoothRhos, smoothDerivs, ok := analyze.Smooth(rawRs, rawRhos, 61)
			if !ok { continue }
			plt.Plot(rawRs, smoothDerivs, plt.LW(3), plt.C(colors[visIdx]))
			r, ok := analyze.SplashbackRadius(rawRs, smoothRhos, smoothDerivs)
			if !ok { continue }
			plt.Plot([]float64{r, r}, []float64{-20, +10}, plt.C(colors[visIdx]))
		}
	}

	// Plot specifications.
	plt.Title(fmt.Sprintf(
		`Halo %d: $M_{\rm 200c}$ = %.3g $M_\odot/h$`, hp.ID(), m),
	)
	plt.XLabel(`$R$ $[{\rm Mpc}/h]$`, plt.FontSize(16))
	plt.YLabel(`$d \ln{\rho}/ d\ln{r}$ [$\rho_m$]`, plt.FontSize(16))

	plt.XScale("log")
	plt.YLim(-20, +10)
	setXRange(rs[0], rs[len(rs) - 1])

	plt.Grid(plt.Axis("y"))
	plt.Grid(plt.Axis("x"), plt.Which("both"))
	plt.SaveFig(fname)
}

func plotPlane(r *analyze.RingBuffer, m float64, id, ring int, dir string) {
	fname := path.Join(dir, fmt.Sprintf("plane_h%d_r%d.png", id, ring))
	xs := make([]float64, 0, r.N)
	ys := make([]float64, 0, r.N)

	for i := 0; i < r.N; i++ {
		if r.Oks[i] {
			xs = append(xs, r.PlaneXs[i])
			ys = append(ys, r.PlaneYs[i])
		}
	}
	
	validRs, validPhis := []float64{}, []float64{}
	for i := range r.Rs {
		if r.Oks[i] {
			validRs = append(validRs, r.Rs[i])
			validPhis = append(validPhis, r.Phis[i])
		}
	}
	kt := analyze.NewKDETree(validRs, validPhis, 4)
	rf := kt.GetRFunc(4)
	spXs := make([]float64, 200)
	spYs := make([]float64, 200)
	dPhi := 2 * math.Pi / float64(len(spXs))
	for i := range spXs {
		r := rf((float64(i) + 0.5) * dPhi)
		cos, sin := math.Sincos(validPhis[i])
		spXs[i], spYs[i] = r * cos, r * sin
	}
	spXs[len(spXs) - 1], spYs[len(spYs) - 1] = spXs[0], spYs[0]


	plt.Figure(plt.FigSize(8, 8))
	plt.Plot(xs, ys, "ow")
	// Plot the colored profiles.
	for i := 0; i < r.N; i++ {
		if r.Oks[i] {
			for visIdx, j := range visProfs {
				if j == i { 
					plt.Plot(
						r.PlaneXs[i], r.PlaneYs[i],
						"o", plt.Color(colors[visIdx]),
					)
				}
			}
		}
	}

	plt.Plot(spXs, spYs, plt.Color("r"), plt.LW(2))

	plt.Title(fmt.Sprintf(
		`Halo %d: $M_{\rm 200c}$ = %.3g $M_\odot/h$`, id, m),
	)
	plt.XLabel(`$X_1$ $[{\rm Mpc}/h]$`, plt.FontSize(16))
	plt.YLabel(`$X_2$ $[{\rm Mpc}/h]$`, plt.FontSize(16))
	rMax := 0.0
	for _, r := range r.Rs {
		if r > rMax { rMax = r }
	}
	plt.XLim(-rMax, +rMax)
	plt.YLim(-rMax, +rMax)
	plt.SaveFig(fname)
}

func setXRange(xLow, xHigh float64) {
	if (xLow < 1 && xHigh  > 1) ||
		(xLow < 0.1 && xHigh > 0.1) || 
		(xLow < 0.01 && xHigh > 0.01) {
		plt.XLim(xLow, xHigh)
	}
}

func strSlice(xs []float64) string {
	tokens := make([]string, len(xs))
	for i := range tokens { tokens[i] = fmt.Sprintf("%.4g", xs[i]) }
	return fmt.Sprintf("[%s]", strings.Join(tokens, ","))
}

// fileNames returns the names of all the files in a directory.
func fileNames(dirName string) ([]string, error) {
	infos, err := ioutil.ReadDir(dirName)
	if err != nil { return nil, err }

	files := make([]string, len(infos))
	for i := range infos {
		files[i] = path.Join(dirName, infos[i].Name())
	}
	return files, nil
}

// intersectingSheets returns all the SheetHeaders and file names that intersect
// with a given halo.
func intersectingSheets(
	h *los.HaloProfiles, hds []io.SheetHeader, files []string,
) ([]io.SheetHeader, []string) {
	hdOuts, fileOuts := []io.SheetHeader{}, []string{}
	for i := range hds {
		if h.SheetIntersect(&hds[i]) {
			hdOuts = append(hdOuts, hds[i])
			fileOuts = append(fileOuts, files[i])
		}
	}
	return hdOuts, fileOuts
}

// intersectionTest is kind of a BS function.
func intersectionTest(
	h *los.HaloProfiles, hds []io.SheetHeader, files []string,
	buf *los.Buffers, subhalos []geom.Sphere,
) {
	hs := []los.HaloProfiles{*h}
	h = &hs[0]

	for i, file := range files {
		hd := &hds[i]
		fmt.Printf("    Reading %s -> ", path.Base(file))
		los.WrapHalo(hs, hd)

		t1 := float64(time.Now().UnixNano())
		buf.ParallelRead(file, hd)
		t2 := float64(time.Now().UnixNano())
		buf.ParallelDensity(h)
		t3 := float64(time.Now().UnixNano())

		fmt.Printf("Setup: %.3g s  Density: %.3g s\n",
			(t2 - t1) / 1e9, (t3 - t2) / 1e9)
	}
	//fmt.Printf("    Rho: %.3g\n", h.Rho(subhalos...))
}
