package main
import (
	"flag"
	"math"
	"fmt"
	"encoding/binary"
	"strings"
	"strconv"
	"path"
	"log"
	"os"
	"io/ioutil"

	"github.com/phil-mansfield/gotetra/density"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/catalog"
	"github.com/phil-mansfield/gotetra/scripts/helper"

	"github.com/phil-mansfield/num/rand"
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
		createCatalog, compareDensity bool
		cells int
	)

	outPath := flag.String("Log", "", "Location to write log statements to. " + 
		"Default is /dev/null.")

	flag.IntVar(&x, "X", -1, "z location of operation.")
	flag.IntVar(&y, "Y", -1, "y location of operation.")
	flag.IntVar(&z, "Z", -1, "z location of operation.")
	flag.IntVar(&cells, "Cells", -1, "Width of grid in cells.")

	flag.BoolVar(&createCatalog, "CreateCatalog", false,
		"Generate gotetra catalogs from gadget catalogs.")
	flag.BoolVar(&compareDensity, "CompareDensity", false,
		"Compare different methods of calculating densities.")

	flag.Parse()

	if *outPath == "" {
		log.SetOutput(ioutil.Discard)
	} else {
		if lf, err := os.Create(*outPath); err != nil {
			log.Fatalln(err.Error())
		} else {
			log.SetOutput(lf)
			defer lf.Close()
		}
	}

	modeName := checkMode(createCatalog, compareDensity)

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
		
	case compareDensity:
		checkCells(cells, modeName)
		checkCoords(x, y, x, modeName)

		args := flag.Args()
		if len(args) != 1 {
			log.Fatalf("Mode %s requires target directory.", modeName)
		}

		source := args[len(args) - 1]
		compareDensityMain(x, y, z, cells, source)
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

func checkMode(createCatalog, compareDensity bool) string {
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

func compareDensityMain(x, y, z, cells int, sourceDir string) {
	log.Printf("Starting to read catalogs at (%d, %d, %d)\n", x, y, z)
	h0, man, centerPs := helper.ReadCatalogs(sourceDir, x, y, z, 1)
	log.Println("Finished reading catalogs.")

	c := &density.Cell{cells, x, y, z}

	ngpGs := make([]density.Grid, 1)
	cicGs := make([]density.Grid, 1)
	cCenterGs := make([]density.Grid, 1)
	mc10Gs := make([]density.Grid, 1)
	mc100Gs := make([]density.Grid, 1)
	mc1000Gs := make([]density.Grid, 1)

	ngpGs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	cicGs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	cCenterGs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	mc10Gs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	mc100Gs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	mc1000Gs[0].Init(h0.TotalWidth, int(h0.GridWidth),
		make([]float64, cells * cells * cells), c)
	
	ngpIntr := density.NearestGridPoint()
	cicIntr := density.CloudInCell()
	cCenterIntr := density.CellCenter(man, h0.CountWidth)
	mc10Intr := density.MonteCarlo(man, h0.CountWidth,
		rand.NewTimeSeed(rand.Tausworthe), 10)
	mc100Intr := density.MonteCarlo(man, h0.CountWidth,
		rand.NewTimeSeed(rand.Tausworthe), 100)
	mc1000Intr := density.MonteCarlo(man, h0.CountWidth,
		rand.NewTimeSeed(rand.Tausworthe), 1000)

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
			//ngpIntr.Interpolate(ngpGs, h0.Mass, idsBuf, xsBuf)
			//cicIntr.Interpolate(cicGs, h0.Mass, idsBuf, xsBuf)
			cCenterIntr.Interpolate(cCenterGs, h0.Mass, idsBuf, xsBuf)
			//mc10Intr.Interpolate(mc10Gs, h0.Mass, idsBuf, xsBuf)
			//mc100Intr.Interpolate(mc100Gs, h0.Mass, idsBuf, xsBuf)
			//mc1000Intr.Interpolate(mc1000Gs, h0.Mass, idsBuf, xsBuf)
			bufIdx = 0
		} else {
			bufIdx++
		}
	}
	ngpIntr.Interpolate(ngpGs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	cicIntr.Interpolate(cicGs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	cCenterIntr.Interpolate(cCenterGs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	mc10Intr.Interpolate(mc10Gs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	mc100Intr.Interpolate(mc100Gs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])
	mc1000Intr.Interpolate(mc1000Gs, h0.Mass, idsBuf[0: bufIdx], xsBuf[0: bufIdx])


	log.Println("Finished interpolation.")

	fmt.Printf("# %12s %12s %12s %12s %12s %12s\n",
		"NGP", "CIC", "Center", "MC - 10", "MC - 100", "MC - 1000")
	for i := 0; i < cells * cells * cells; i++ {
		fmt.Printf("  %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",
			ngpGs[0].Rhos[i], cicGs[0].Rhos[i], cCenterGs[0].Rhos[i],
			mc10Gs[0].Rhos[i], mc100Gs[0].Rhos[i], mc1000Gs[0].Rhos[i])
	}

	log.Println("Finished printing density grid.")
}
