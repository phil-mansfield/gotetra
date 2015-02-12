package main
import (
	"flag"
	"fmt"
	"path"
	"encoding/binary"
	"strings"
	"strconv"
	"log"
	"runtime"
	"runtime/pprof"
	"os"
	"io/ioutil"

	"github.com/phil-mansfield/gotetra"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/io"
)

const (
	catalogBufLen = 1<<12
)

var (
	gadgetEndianness = binary.LittleEndian
)

func main() {
	var (
		logPath, pprofPath, boundsFile string
		minSheet, maxSheet int
		createSheet, fullDensity, boundedDensity bool
		points, cells, skip int
	)

	flag.StringVar(&logPath, "Log", "",
		"Location to write log statements to. Default is stderr.")
	flag.StringVar(&pprofPath, "PProf", "",
		"Location to write profile to. Default is no profiling.")
	flag.StringVar(&boundsFile, "BoundsFile", "",
		"Location to read in grid bounds from.")

	flag.IntVar(&cells, "Cells", -1, "Width of grid in cells.")
	flag.IntVar(&skip, "Skip", 1, "Width of a tetrahedron in particles.")
	flag.IntVar(&points, "Points", 100, "Number of points to use for method.")

	flag.IntVar(&minSheet, "MinSheet", 0,
		"Lowest index of a sheet file that will be read.")
	flag.IntVar(&maxSheet, "MaxSheet", 511,
		"Highest index of a sheet file that will be read.")


	flag.BoolVar(&fullDensity, "FullDensity", false,
		"Compute full-box density.")
	flag.BoolVar(&boundedDensity, "BoundedDensity", false,
		"Computer density for  sequence of sub-boxes.")
	flag.BoolVar(&createSheet, "CreateSheet", false,
		"Generate gotetra sheets from gadget catalogs.")

	flag.Parse()

	if logPath != "" {
		if lf, err := os.Create(logPath); err != nil {
			log.Fatalln(err.Error())
		} else {
			log.SetOutput(lf)
			defer lf.Close()
		}
	}

	if pprofPath != "" {
		f, err := os.Create(pprofPath)
		if err != nil { log.Fatal(err) }
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	// This is terribly designed.
	modeName := checkMode(createSheet, fullDensity, boundedDensity)

	switch {
	case createSheet:
		checkCells(cells, modeName)

		args := flag.Args()
		if len(args) <= 1 {
			log.Fatalf(
				"Mode %s requires source files and a target directory.",
				modeName,
			)
		}

		sources := args[0: len(args) - 1]
		targetDir := args[len(args) - 1]

		createSheetMain(cells, sources, targetDir)
	case fullDensity:
		checkCells(cells, modeName)
		args := flag.Args()
		if len(args) != 2 {
			log.Fatalf("Mode %s requires a source and target directory.",
				modeName)
		}

		source := args[0]
		target := args[1]
		fullDensityMain(cells, points, skip, source, target,
			minSheet, maxSheet)
	case boundedDensity:
		checkCells(cells, modeName)
		args := flag.Args()
		if len(args) != 2 {
			log.Fatalf("Mode %s requires a source and target directory.",
				modeName)
		}

		source := args[0]
		target := args[1]
		boundedDensityMain(cells, points, skip, source, target,
			boundsFile, minSheet, maxSheet)
	}
}

func checkCells(cells int, modeName string) {
	if cells == -1 {
		log.Fatalf(
			"The mode %s requires a positive number of cells.\n", modeName,
		)
	}
}

func checkMode(createSheet, fullDensity, boundedDensity bool) string {
	// This function is not a good idea.
	n := 0
	modeStr := ""

	if createSheet {
		n++
		modeStr = "CreateSheet"
	}

	if fullDensity {
		n++
		modeStr = "FullDensity"
	}

	if boundedDensity {
		n++
		modeStr = "BoundedDensity"
	}

	if n != 1 {
		log.Fatalf("Given %d mode flags, but exactly 1 is required.\n", n)
	}

	return modeStr
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

func createSheetMain(cells int, matches []string, outPath string) {
	hd, xs, vs := createGrids(matches)

	// Parse, parse, parse, parse...
	snapDir := path.Base(path.Dir(matches[0]))
	paramDir := path.Base(path.Dir(path.Dir(path.Dir(matches[0]))))

	l, n, str, err := parseDir(paramDir)
	if err != nil {
		log.Fatalf(err.Error())
	}
	hd.CountWidth = int64(n)

	outParamDir := fmt.Sprintf("Box_L%04d_N%04d_G%04d_%s", l, n, cells, str)
	outParamPath := path.Join(outPath, outParamDir)
	outDir := path.Join(outParamPath, snapDir)
	if err = os.MkdirAll(outDir, 0777); err != nil {
		log.Fatalf(err.Error())
	}

	writeGrids(outDir, hd, cells, xs, vs)
}

func createGrids(catalogs []string) (hd *io.CatalogHeader, xs, vs []geom.Vec) {
	hs := make([]io.CatalogHeader, len(catalogs))
	for i := range hs {
		hs[i] = *io.ReadGadgetHeader(catalogs[i], gadgetEndianness)
	}

	maxLen := int64(0)
	for _, h := range hs {
		if h.Count > maxLen {
			maxLen = h.Count
		}
	}

	xs = make([]geom.Vec, hs[0].TotalCount)
	vs = make([]geom.Vec, hs[0].TotalCount)
	
	pBuf := make([]io.Particle, maxLen)
	intBuf := make([]int64, maxLen)
	floatBuf := make([]float32, maxLen * 3)

	buf := io.NewParticleBuffer(xs, vs, catalogBufLen)

	for i, cat := range catalogs {

		if i % 25 == 0 {
			log.Printf("Read %d/%d catalogs", i, len(catalogs))
		}

		N := hs[i].Count
		floatBuf = floatBuf[0: 3*N]
		intBuf = intBuf[0: N]
		pBuf = pBuf[0: N]

		io.ReadGadgetParticlesAt(cat, gadgetEndianness,
			floatBuf, intBuf, pBuf)
		runtime.GC()
		buf.Append(pBuf)
	}

	if len(catalogs) % 25 != 0 {
		log.Printf("Read %d/%d catalogs", len(catalogs), len(catalogs))
	}

	return &hs[0], xs, vs
}

func writeGrids(outDir string, hd *io.CatalogHeader,
	cells int, xs, vs []geom.Vec) {

	log.Println("Writing to directory", outDir)

	segmentWidth := int(hd.CountWidth) / cells
	gridWidth := segmentWidth + 1

	xsSeg := make([]geom.Vec, gridWidth * gridWidth * gridWidth)
	vsSeg := make([]geom.Vec, gridWidth * gridWidth * gridWidth)

	shd := &io.SheetHeader{}
	shd.Cosmo = hd.Cosmo
	shd.CountWidth = hd.CountWidth
	shd.Mass = hd.Mass
	shd.TotalWidth = hd.TotalWidth

	shd.SegmentWidth = int64(segmentWidth)
	shd.GridWidth = int64(gridWidth)
	shd.GridCount = int64(shd.GridWidth * shd.GridWidth * shd.GridWidth)
	shd.Cells = int64(cells)

	for z := int64(0); z < shd.Cells; z++ {
		for y := int64(0); y < shd.Cells; y++ {
			for x := int64(0); x < shd.Cells; x++ {
				copyToSegment(shd, xs, vs, xsSeg, vsSeg)
				file := path.Join(outDir, fmt.Sprintf("sheet%d%d%d.dat", x, y,z))
				io.WriteSheet(file, shd, xsSeg, vsSeg)
				runtime.GC()

				if shd.Idx % 25 == 0 {
					log.Printf("Wrote %d/%d sheet segments.",
						shd.Idx, shd.Cells * shd.Cells * shd.Cells,
					)
				}
				shd.Idx++
			}
		}
	}
	if shd.Idx % 25 != 0 {
		log.Printf("Wrote %d/%d sheet segments.",
			shd.Idx, shd.Cells * shd.Cells * shd.Cells,
		)
	}
}

// Note, this only works for the collections of points where each point is
// relatively close to the existing bounding box.
type boundingBox struct {
	Width float64
	Center geom.Vec
	ToMax, ToMin, ToPt geom.Vec
}

func (box *boundingBox) Init(pt *geom.Vec, width float64) {
	box.Width = width
	box.Center = *pt
}

func (box *boundingBox) Add(pt *geom.Vec) {
	pt.SubAt(&box.Center, box.Width, &box.ToPt)

	for i := 0; i < 3; i++ {
		box.ToMin[i], box.ToMax[i] = minMax(
			box.ToMin[i], box.ToMax[i], box.ToPt[i],
		)
	}

	//ToPt is now a buffer for the neccesary shift in the center
	box.ToMax.AddAt(&box.ToMin, &box.ToPt)
	box.ToPt.ScaleSelf(0.5)
	box.Center.AddSelf(&box.ToPt)
    box.ToMax.SubSelf(&box.ToPt, box.Width)
    box.ToMin.SubSelf(&box.ToPt, box.Width)
}


func minMax(min, max, x float32) (outMin, outMax float32) {
	if x > max {
		return min, x
	} else if x < min {
		return x, max
	} else {
		return min, max
	}
}

func copyToSegment(shd *io.SheetHeader, xs, vs, xsSeg, vsSeg []geom.Vec) {
	xStart := shd.SegmentWidth * (shd.Idx % shd.Cells)
	yStart := shd.SegmentWidth * ((shd.Idx / shd.Cells) % shd.Cells)
	zStart := shd.SegmentWidth * (shd.Idx / (shd.Cells * shd.Cells))

	N, N2 := shd.CountWidth, shd.CountWidth * shd.CountWidth

	box := &boundingBox{}
	box.Init(&xs[xStart + N * yStart + N2 * zStart], shd.TotalWidth)

	smallIdx := 0
	for z := zStart; z < zStart + shd.GridWidth; z++ {
		zIdx := z
		if zIdx == shd.CountWidth { zIdx = 0 }
		for y := yStart; y < yStart + shd.GridWidth; y++ {
			yIdx := y
			if yIdx == shd.CountWidth { yIdx = 0 }
			for x := xStart; x < xStart + shd.GridWidth; x++ {
				xIdx := x
				if xIdx == shd.CountWidth { xIdx = 0 }

				largeIdx := xIdx + yIdx * N + zIdx * N2
				
				xsSeg[smallIdx] = xs[largeIdx]
				vsSeg[smallIdx] = vs[largeIdx]

				box.Add(&xsSeg[smallIdx])

				smallIdx++
			}
		}	
	}
	
	box.Center.AddAt(&box.ToMin, &shd.Origin)
	shd.Origin.ModSelf(shd.TotalWidth)
	box.ToMax.ScaleAt(2.0, &shd.Width)
}

func validCellNum(cells int) bool {
	for cells > 1 { cells /= 2 }
	return cells == 1
}

func fullDensityMain(
	cells, points, skip int,
	sourceDir, outDir string,
	minSheet, maxSheet int,
) {
	log.Println("Running FullDensity main.")

	infos, err := ioutil.ReadDir(sourceDir)
	if err != nil { log.Fatalf(err.Error()) }
	files := make([]string, len(infos))
	for i := range infos { files[i] = path.Join(sourceDir, infos[i].Name()) }

	// gotetra code starts here

	boxes := make([]gotetra.Box, 1)
	boxes[0].InitFullFromFile(files[0], cells)	
	man := gotetra.NewFullManager(files, cells, points)
	man.Subsample(skip)

	for i := minSheet; i <= maxSheet; i++ {
		log.Println("Analyzing file", infos[i].Name())
		man.LoadPositions(files[i])
		man.Density(boxes)
	}

	// gotetra code ends here

    out := path.Join(outDir,
        fmt.Sprintf("G%d_NP%d_S%d.dat", cells, points, skip))
    log.Printf("Writing %s", out)
    writeDensity(out, boxes[0].Vals) 
}

func boundedDensityMain(
	cells, points, skip int,
	sourceDir, outDir, boundsFile string,
	minSheet, maxSheet int,
) {
	log.Println("Running BoundedDensity main.")

	infos, err := ioutil.ReadDir(sourceDir)
	if err != nil { log.Fatalf(err.Error()) }
	files := make([]string, len(infos))
	for i := range infos { files[i] = path.Join(sourceDir, infos[i].Name()) }

	bs, err := ioutil.ReadFile(boundsFile)
	if err != nil { log.Fatalf(err.Error()) }
	origins, widths := readCols(string(bs))

	// Biggest:
	// radius := 2.170
	// center := [3]float64{4.602, 100.704, 80.701}
	// 100th Biggest:
	// radius := 0.791
	// center := [3]float64{29.609, 64.128, 42.843}

	// gotetra code starts here
	boxes := make([]gotetra.Box, len(origins))
	for i := range boxes {
		boxes[i].Init(125, cells, origins[i], widths[i])
	}

	man := gotetra.NewBoundedManager(files, boxes, points)
	man.Subsample(skip)

	for i := minSheet; i <= maxSheet; i++ {
		man.LoadPositions(files[i])
		man.Density(boxes)
	}

	// gotetra code ends here

	for i := range boxes {
		boxWidth := boxes[i].CellWidth() * float64(cells)
		out := path.Join(outDir,
			fmt.Sprintf("L%g_G%d_NP%d_S%d_XW%d_YW%d_ZW%d.dat",
				boxWidth, cells, points, skip,
				boxes[i].Width(0), boxes[i].Width(1), boxes[i].Width(2)))
		writeDensity(out, boxes[i].Vals) 
	}
}

func readCols(s string) (origins, widths [][3]float64) {	
	lines := strings.Split(s, "\n")

	if len(lines[len(lines) - 1]) == 0 { lines = lines[:len(lines) - 1] }

	origins = make([][3]float64, len(lines))
	widths = make([][3]float64, len(lines))

	for i, line := range lines {
		tokens := strings.Split(line, " ")
		words := []string{}

		for _, tok := range tokens {
			if len(tok) > 0 { words = append(words, tok)}
		}

		if len(words) != 6 {
			log.Fatalf("Found %d tokens on line %d of BoundsFile. " + 
				"Require 6 tokens.", len(words), i)
		}

		ox, err := strconv.ParseFloat(words[0], 64)
		if err != nil { log.Fatalf(err.Error()) }		
		oy, err := strconv.ParseFloat(words[1], 64)
		if err != nil { log.Fatalf(err.Error()) }
		oz, err := strconv.ParseFloat(words[2], 64)
		if err != nil { log.Fatalf(err.Error()) }

		wx, err := strconv.ParseFloat(words[3], 64)
		if err != nil { log.Fatalf(err.Error()) }
		wy, err := strconv.ParseFloat(words[4], 64)
		if err != nil { log.Fatalf(err.Error()) }
		wz, err := strconv.ParseFloat(words[5], 64)
		if err != nil { log.Fatalf(err.Error()) }

		origins[i] = [3]float64{float64(ox), float64(oy), float64(oz)}
		widths[i] = [3]float64{float64(wx), float64(wy), float64(wz)}
	}

	return origins, widths
}

func writeDensity(fname string, d []float64) {
	f, err := os.Create(fname)
	if err != nil { log.Fatal(err.Error()) }
	defer f.Close()

	binary.Write(f, binary.LittleEndian, d)
}
