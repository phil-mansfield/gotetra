package main
import (
	"flag"
	"math"
	"fmt"
	"path"
	"encoding/binary"
	"strings"
	"strconv"
	"log"
	"runtime"
	"runtime/pprof"
	"os"

	"github.com/phil-mansfield/gotetra/density"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/catalog"
	"github.com/phil-mansfield/gotetra/scripts/helper"

	"github.com/phil-mansfield/gotetra/rand"
)

const (
	vecBufLen = 1<<10
	catalogBufLen = 1<<12
)

var (
	gadgetEndianness = binary.LittleEndian
)

func main() {
	var (
		x, y, z int
		createSheet, createCatalog, compareDensity, densityMode, tetraStats bool
		cells int
	)

	logPath := flag.String("Log", "", "Location to write log statements to. " + 
		"Default is stderr.")
	pprofPath := flag.String("PProf", "", "Location to write profile to.")

	flag.IntVar(&x, "X", -1, "z location of operation.")
	flag.IntVar(&y, "Y", -1, "y location of operation.")
	flag.IntVar(&z, "Z", -1, "z location of operation.")
	flag.IntVar(&cells, "Cells", -1, "Width of grid in cells.")

	flag.BoolVar(&createCatalog, "CreateCatalog", false,
		"Generate gotetra catalogs from gadget catalogs.")
	flag.BoolVar(&createSheet, "CreateSheet", false,
		"Generate gotetra sheets from gadget catalogs.")
	flag.BoolVar(&compareDensity, "CompareDensity", false,
		"Compare different methods of calculating densities.")
	flag.BoolVar(&densityMode, "Density", false,
		"Compute density of a cell. Requires Method flag to be set.")
	flag.BoolVar(&tetraStats, "TetraStats", false,
		"Print basic geometry statistics about the tetrahedra in a given cell.")

	method := flag.String("Method", "", "Estimator method. Can be set to" + 
		"(CloudInCell | NearestGridPoint | MonteCarlo | SubRandom)." + 
		"MonteCarlo and SubRandom also allow the Points flag to be set.")
	points := flag.Int("Points", 100, "Number of points to use for method.")
	pointSelector := flag.String("PointSelector", "Flat", "Method used to " + 
		"select the number of points per tetrahedron in Monte Carlo " + 
		"integration. Can be (Flat | PropToCells)")

	flag.Parse()

	if *logPath != "" {
		if lf, err := os.Create(*logPath); err != nil {
			log.Fatalln(err.Error())
		} else {
			log.SetOutput(lf)
			defer lf.Close()
		}
	}

	if *pprofPath != "" {
		f, err := os.Create(*pprofPath)
		if err != nil { log.Fatal(err) }
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	modeName := checkMode(createCatalog, compareDensity,
		tetraStats, densityMode)

	switch {
	case createCatalog:
		checkCells(cells, modeName)

		args := flag.Args()
		if len(args) <= 1 {
			log.Fatalf(
				"Mode %s requires source files and a target directory.",
				modeName,
			)
		}

		sources := args[1: len(args) - 1]
		targetDir := args[len(args) - 1]

		createCatalogsMain(cells, sources, targetDir)
	case densityMode:
		checkCells(cells, modeName)
		checkCoords(x, y, z, modeName)

		args := flag.Args()
		if len(args) != 2 {
			log.Fatalf("Mode %s requires a source and target directory.",
				modeName)
		}
	
		psFlag := density.PointSelectorFromString(*pointSelector)

		source := args[0]
		target := args[1]
		densityMain(x, y, z, cells, *method, *points, psFlag,
			source, target)
	case compareDensity:
		checkCells(cells, modeName)
		checkCoords(x, y, x, modeName)

		args := flag.Args()
		if len(args) != 1 {
			log.Fatalf("Mode %s requires target directory.", modeName)
		}

		source := args[len(args) - 1]
		compareDensityMain(x, y, z, cells, source)
	case tetraStats:
		checkCells(cells, modeName)
		checkCoords(x, y, z, modeName)
		args := flag.Args()
		if len(args) != 1 {
			log.Fatalf("Mode %s requires a target directory.", modeName)
		}
		tetraStatsMain(x, y, z, cells, args[len(args) - 1])
	}
}

func checkCells(cells int, modeName string) {
	if cells == -1 {
		log.Fatalf(
			"The mode %s requires a positive number of cells.\n", modeName,
		)
	}
}

func checkCoords(x, y, z int, modeName string) {
	if x == -1 {
		log.Fatalf("The mode %s requires an x location.\n", modeName)
	} else if y == -1 {
		log.Fatalf("The mode %s requires a y location.\n", modeName)
	} else if z == -1 {
		log.Fatalf("The mode %s requires a z location.\n", modeName)
	}
}

func checkMode(createCatalog, compareDensity, tetraStats, density bool) string {
	n := 0
	modeStr := ""

	if createCatalog {
		n++
		modeStr = "CreateCatalog"
	}

	if compareDensity {
		n++
		modeStr = "CompareDensity"
	}

	if tetraStats {
		n++
		modeStr = "TetraStats"
	}

	if density {
		n++
		modeStr = "Density"
	}

	if n != 1 {
		log.Fatalf("Given %d mode flags, but exactly 1 is required.\n", n)
	}

	return modeStr
}

func createCatalogsMain(cells int, matches []string, outPath string) {
	snapDir := path.Base(path.Dir(matches[0]))
	paramDir := path.Base(path.Dir(path.Dir(path.Dir(matches[0]))))

	l, n, str, err := parseDir(paramDir)
	if err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}

	outParamDir := fmt.Sprintf("Box_L%04d_N%04d_G%04d_%s", l, n, cells, str)
	outParamPath := path.Join(outPath, outParamDir)
	outSnapPath := path.Join(outParamPath, snapDir)

	if err = os.MkdirAll(path.Join(outParamPath, snapDir), 0777); err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}

	catalogNames := createCatalogs(outSnapPath, matches[0], cells, n)

	rebinParticles(matches, catalogNames, int64(cells))
}

// pasreDir reads one of Benedikt's sim directory names and returns the relevent
// physical information.
func parseDir(dir string) (int, int, string, error) {
	parts := strings.Split(dir, "_")

	if len(parts) != 4 {
		return dirErr(dir)
	} else if len(parts[1]) != 5 {
		return dirErr(dir)
	} else if len(parts[2]) != 5 {
		return dirErr(dir)
	}

	l, err := strconv.Atoi(parts[1][1:5])
	if err != nil {
		return 0, 0, "", err
	}
	n, err := strconv.Atoi(parts[2][1:5])
	if err != nil {
		return 0, 0, "", err
	}

	return l, n, parts[3], nil
}

func dirErr(dir string) (int, int, string, error) {
	return 0, 0, "", fmt.Errorf("Invalid source directory '%s'.", dir)
}

// createCatalogs creates new instances of the needed gridcell files in the
// given directory.
func createCatalogs(outPath, exampleInput string, gridWidth, countWidth int) []string {
	h := catalog.ReadGadgetHeader(exampleInput, gadgetEndianness)
	h.Width = h.TotalWidth / float64(gridWidth)
	h.CountWidth = int64(countWidth)
	h.GridWidth = int64(gridWidth)

	ps := []catalog.Particle{}

	maxIdx := gridWidth * gridWidth * gridWidth
	catalogNames := []string{}
	for i := 0; i < maxIdx; i++ {
		name := path.Join(outPath, fmt.Sprintf("gridcell_%04d.dat", i))
		catalogNames = append(catalogNames, name)
		h.Idx = int64(i)
		catalog.Write(name, h, ps)
	}

	return catalogNames
}

// rebinParticles transfers particles from a slice of Gadget catalogs to a
// slice of gotetra catalogs.
func rebinParticles(inFiles, outFiles []string, gridWidth int64) {
	hs := make([]catalog.Header, len(inFiles))
	log.Println("Creating catalog files")
	for i := range hs {
		hs[i] = *catalog.ReadGadgetHeader(inFiles[i], gadgetEndianness)
	}
	log.Println("Finished creating catalog files")

	width := hs[0].TotalWidth / float64(gridWidth)

	maxLen := int64(0)
	for _, h := range hs {
		if h.Count > maxLen {
			maxLen = h.Count
		}
	}

	floatBufMax := make([]float32, maxLen * 3)
	intBufMax := make([]int64, maxLen)
	psBufMax := make([]catalog.Particle, maxLen)

	bufs := createBuffers(outFiles, catalogBufLen)


	for i, inFile := range inFiles {
		if i % 25 == 0 {
			log.Printf("Converting %d/%d catalogs.\n", i, len(inFiles))
		}

		floatBuf := floatBufMax[0:hs[i].Count*3]
		intBuf := intBufMax[0:hs[i].Count]
		psBuf := psBufMax[0:hs[i].Count]

		catalog.ReadGadgetParticlesAt(
			inFile, gadgetEndianness, floatBuf, intBuf, psBuf,
		)

		rebinSlice(width, gridWidth, bufs, psBuf)
	}

	for _, buf := range bufs {
		buf.Flush()
	}
}

// createBuffers creates buffers for all the catalog files that will be
// witten to.
func createBuffers(paths []string, bufSize int) []*catalog.ParticleBuffer {
	bufs := make([]*catalog.ParticleBuffer, len(paths))

	for i, path := range paths {
		bufs[i] = catalog.NewParticleBuffer(path, bufSize)
	}

	return bufs
}

// rebinSlice rebins the given particles into the grid of files represented
// by bufs.
func rebinSlice(
	width float64,
	gridWidth int64,
	bufs []*catalog.ParticleBuffer,
	ps []catalog.Particle,
) {
	for _, p := range ps {
		xIdx := int64(math.Floor(float64(p.Xs[0]) / width))
		yIdx := int64(math.Floor(float64(p.Xs[1]) / width))
		zIdx := int64(math.Floor(float64(p.Xs[2]) / width))
		idx := xIdx + yIdx*gridWidth + zIdx*gridWidth*gridWidth
		bufs[idx].Append(p)
	}
}

func tetraStatsMain(x, y, z, cells int, sourceDir string) {
	h0, man, centerPs := helper.ReadCatalogs(sourceDir, x, y, z, 0)
	cellWidth := h0.Width / float64(cells)

	idxBuf := &geom.TetraIdxs{}
	tet := &geom.Tetra{}
	misses := 0

	vs := make([]float64, 6 * len(centerPs))
	cvs := make([]float64, 6 * len(centerPs))
	mins := make([]float64, 6 * len(centerPs))
	cb := &geom.CellBounds{}

	for i := range centerPs {
		for dir := 0; dir < 6; dir++ {
			idxBuf.Init(centerPs[i].Id, h0.CountWidth, dir)
			p0 := man.Get(idxBuf[0])
			p1 := man.Get(idxBuf[1])
			p2 := man.Get(idxBuf[2])
			p3 := man.Get(idxBuf[3])
			
			if p0 == nil || p1 == nil || p2 == nil || p3 == nil {
				misses++
				continue
			}

			tet.Init(&p0.Xs, &p1.Xs, &p2.Xs, &p3.Xs, h0.TotalWidth)
			vol := tet.Volume()

			if vol <= 0.0 {
				misses++
				continue
			}

			tet.CellBoundsAt(cellWidth, cb)
			cv := float64((cb.Max[0] - cb.Min[0]) *
				(cb.Max[1] - cb.Min[1]) *
				(cb.Max[2] - cb.Min[2]))

			min, _ := tet.MinMaxLeg()
			if min <= 0.0 {
				misses++
				continue
			}

			vs[6*i + dir - misses] = math.Log10(vol)
			cvs[6*i + dir - misses] = math.Log10(cv)
			mins[6*i + dir - misses] = math.Log10(min)
		}
	}

	vs = vs[0: len(vs) - misses]
	cvs = cvs[0: len(cvs) - misses]
	mins = mins[0: len(mins) - misses]
	vBins, vCounts := binSample(vs)
	cvBins, cvCounts := binSample(cvs)
	minBins, minCounts := binSample(mins)

	for i := range vBins { vBins[i] = math.Pow(10, vBins[i]) }
	for i := range cvBins { cvBins[i] = math.Pow(10, cvBins[i]) }
	for i := range minBins { minBins[i] = math.Pow(10, minBins[i]) }

	fmt.Printf("# Misses: %d Total: %d Frac %g\n",
		misses, 6 * h0.Count, float64(misses) / float64(h0.Count) / 6)

	for i := range vBins {
		fmt.Printf("%10.4g %8d %10.5g %8d %10.4g %8d\n",
			vBins[i], vCounts[i], cvBins[i], cvCounts[i],
			minBins[i], minCounts[i])
	}
}

func binSample(xs []float64) ([]float64, []int) {
	minX, maxX := xs[0], xs[0]
	for _, x := range xs {
		minX = min(minX, x)
		maxX = max(maxX, x)
	}
	
	n := float64(len(xs))
	binNum := int(math.Ceil(math.Pow(n, 1.0/3.0) / 2))
	binWidth := (maxX - minX) / float64(binNum)
	
	bins := make([]float64, binNum)
	counts := make([]int, binNum)

	for i := range bins { bins[i] = (float64(i) + 0.5) * binWidth + minX }

	for _, x := range xs {
		idx := int((x - minX) / binWidth)
		if idx == len(counts) { idx-- }
		counts[idx]++
	}

	return bins, counts
}

func min(x, y float64) float64 {
	if x < y { return x }
	return y
}

func max(x, y float64) float64 {
	if x > y { return x }
	return y
}

func densityMain(x, y, z, cells int, method string,
	points int, psFlag density.PointSelectorFlag, sourceDir, outDir string) {

	if points <= 0 {
		log.Fatalf("Positive value for Points is required.")
	}

	// This is done twice. Don't do this the second time you write this.
	switch method {
	case "NearestGridPoint":
	case "CloudInCell":
	case "MonteCarlo":
	case "SubRandom":
	default:
		log.Fatalf("Unrecognized method string %s.", method)
	}

	runtime.GOMAXPROCS(runtime.NumCPU())
	workers := runtime.GOMAXPROCS(0)
	log.Printf("GOMAXPROCS set to %d.", workers)
	
	h0, man, centerPs := helper.ReadCatalogs(sourceDir, x, y, z, 1)
	c := &density.Cell{cells, x, y, z}
	out := make(chan []density.Grid, workers)

	blockLen := len(centerPs) / workers
	start, end := 0, 0

	for wID := 0; wID < workers; wID++ {
		// This is horrifying:
		gs, intr := setupIntr(method, h0, man, cells, points, psFlag, c)

		start, end = end, end + blockLen
		if wID == 0 { end += len(centerPs) % workers }

		go densityRoutine(h0, man, centerPs[start: end], intr, gs, wID, out)
		start = end
	}

	ds := make([][]float64, 1)
	for i := range ds {
		ds[i] = make([]float64, cells * cells * cells)
	}

	for w := 0; w < workers; w++ {
		addDensityGrid(ds, <- out)
	}

	writeDensity(outDir, x, y, z, ds[0])
}

func addDensityGrid(ds [][]float64, gs []density.Grid) {
	if len(ds) != len(gs) {
		log.Fatalf("Length of ds = %d, but length of gs = %d.",
			len(ds), len(gs))
	}

	for i, d := range ds {
		for j := range d {
			d[j] += gs[i].Rhos[j]
		}
	}
}

func setupIntr(method string, h0 *catalog.Header, man *catalog.ParticleManager,
	cells, points int, psFlag density.PointSelectorFlag,
	c *density.Cell) ([]density.Grid, density.Interpolator) {

	gs := make([]density.Grid, 1)
	for i := range gs {
		gs[i].Init(h0.TotalWidth, int(h0.GridWidth),
			make([]float64, cells * cells * cells), c)
	}
	
	var intr density.Interpolator
	
	switch method {
	case "NearestGridPoint": 
		intr = density.NearestGridPoint()
	case "CloudInCell":
		intr = density.CloudInCell()
	case "MonteCarlo":
		intr = density.MonteCarlo(man, h0.CountWidth,
			rand.NewTimeSeed(rand.Xorshift), points, psFlag)
	case "SubRandom":
		intr = density.SobolSequence(man, h0.CountWidth, points)
	default:
		log.Fatalf("Unrecognized method string %s.", method)
	}
	
	return gs, intr
}

func writeDensity(outDir string, x, y, z int, d []float64) {
	fname := fmt.Sprintf("density_%d%d%d.dat", x, y, z)
	path := path.Join(outDir, fname)

	f, err := os.Create(path)
	if err != nil { log.Fatal(err.Error()) }
	defer f.Close()

	binary.Write(f, binary.LittleEndian, d)
}

func densityRoutine(h0 *catalog.Header, man *catalog.ParticleManager,
	centerPs []catalog.Particle, intr density.Interpolator, gs []density.Grid,
	wID int, out chan<- []density.Grid) {

	xsBuf := make([]geom.Vec, vecBufLen)
	idsBuf := make([]int64, vecBufLen)
		
	log.Printf("Set up interpolation worker %d.", wID)
	checkLen := len(centerPs) / 5

	bufIdx := 0
	for i := range centerPs {
		if i % checkLen == 0 {
			log.Printf("%d Particles interpolated on worker %d (%d/%d).\n",
				i, wID, i / checkLen, len(centerPs) / checkLen)
		}

		xsBuf[bufIdx] = centerPs[i].Xs
		idsBuf[bufIdx] = centerPs[i].Id

		if bufIdx == len(xsBuf) - 1 {
			intr.Interpolate(gs, h0.Mass, idsBuf, xsBuf)
			bufIdx = 0
		} else {
			bufIdx++
		}
	}
	intr.Interpolate(gs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	out <- gs
}

func compareDensityMain(x, y, z, cells int, sourceDir string) {
	log.Printf("Starting to read catalogs at (%d, %d, %d)\n", x, y, z)
	h0, man, centerPs := helper.ReadCatalogs(sourceDir, x, y, z, 1)
	log.Println("Finished reading catalogs.")

	c := &density.Cell{cells, x, y, z}

	ngpGs := make([]density.Grid, 1)
	cicGs := make([]density.Grid, 1)
	cCenterGs := make([]density.Grid, 1)
	mc100Gs := make([]density.Grid, 1)
	seq100Gs := make([]density.Grid, 1)
	mc2500Gs := make([]density.Grid, 1)

	ngpGs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	cicGs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	cCenterGs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	mc100Gs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	seq100Gs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	mc2500Gs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	
	ngpIntr := density.NearestGridPoint()
	cicIntr := density.CloudInCell()
	cCenterIntr := density.CellCenter(man, h0.CountWidth)
	mc100Intr := density.MonteCarlo(man, h0.CountWidth,
		rand.NewTimeSeed(rand.Tausworthe), 100, density.Flat)
	seq100Intr := density.SobolSequence(man, h0.CountWidth, 100)
	mc2500Intr := density.MonteCarlo(man, h0.CountWidth,
		rand.NewTimeSeed(rand.Tausworthe), 2500, density.Flat)

	xsBuf := make([]geom.Vec, vecBufLen)
	idsBuf := make([]int64, vecBufLen)
	
	log.Println("Set up interpolators and buffers.")
	checkLen := len(centerPs) / 20

	bufIdx := 0
	for i := range centerPs {
		if i % checkLen == 0 {
			log.Printf("%d Particles interpolated.\n", i)
		}

		xsBuf[bufIdx] = centerPs[i].Xs
		idsBuf[bufIdx] = centerPs[i].Id

		if bufIdx == len(xsBuf) - 1 {
			ngpIntr.Interpolate(ngpGs, h0.Mass, idsBuf, xsBuf)
			cicIntr.Interpolate(cicGs, h0.Mass, idsBuf, xsBuf)
			cCenterIntr.Interpolate(cCenterGs, h0.Mass, idsBuf, xsBuf)
			mc100Intr.Interpolate(mc100Gs, h0.Mass, idsBuf, xsBuf)
			seq100Intr.Interpolate(seq100Gs, h0.Mass, idsBuf, xsBuf)
			mc2500Intr.Interpolate(mc2500Gs, h0.Mass, idsBuf, xsBuf)
			bufIdx = 0
		} else {
			bufIdx++
		}
	}
	ngpIntr.Interpolate(ngpGs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	cicIntr.Interpolate(cicGs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	cCenterIntr.Interpolate(cCenterGs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	mc100Intr.Interpolate(mc100Gs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	seq100Intr.Interpolate(seq100Gs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	mc2500Intr.Interpolate(mc2500Gs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])


	log.Println("Finished interpolation.")

	fmt.Printf("# %8s %8s %8s %8s %8s %8s\n",
		"NGP", "CIC", "Center", "MC - 100", "SS - 100", "MC - 2500")
	for i := 0; i < cells * cells * cells; i++ {
		fmt.Printf("  %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n",
			ngpGs[0].Rhos[i], cicGs[0].Rhos[i], cCenterGs[0].Rhos[i],
			mc100Gs[0].Rhos[i], seq100Gs[0].Rhos[i], mc2500Gs[0].Rhos[i])
	}

	log.Println("Finished printing density grid.")
}

func createSheetMain(cells int, matches []string, outPath string) {
	hd, xs, vs := createGrids(matches)

	// Parse, parse, parse, parse...
	snapDir := path.Base(path.Dir(matches[0]))
	paramDir := path.Base(path.Dir(path.Dir(path.Dir(matches[0]))))

	l, n, str, err := parseDir(paramDir)
	if err != nil {
		log.Fatalf(err.Error())
	}

	outParamDir := fmt.Sprintf("Box_L%04d_N%04d_G%04d_%s", l, n, cells, str)
	outParamPath := path.Join(outPath, outParamDir)
	outSnapPath := path.Join(outParamPath, snapDir)

	outDir := path.Join(outParamPath, snapDir)
	if err = os.MkdirAll(outDir, 0777); err != nil {
		log.Fatalf(err.Error())
	}

	writeGrids(outDir, hd, cells, xs, vs)
}

func createGrids(catalogs []string) (hd *catalog.Header, xs, vs []float64) {
	hd := catalog.ReadGadgetheader(catalogs[0], gadgetEndianness)
	xs := make([]geom.Vec, hd.TotalCount)
	vs := make([]geom.Vec, hd.TotalCount)
	
	fmt.Println("Okay, we're good to go.")
	os.Exit(1)
}

func writeGrids(outDir string, hd *catalog.Header, cells int, xs, vs []float64) {
}
