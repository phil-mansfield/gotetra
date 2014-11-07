package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"math"
	"os"
	"path"
	"strconv"
	"strings"

	tetra "github.com/phil-mansfield/gotetra"
	"github.com/phil-mansfield/gotetra/catalog"
)

var (
	defaultBufSize   = 1 << 12
	gadgetEndianness = binary.LittleEndian
)

func main() {
	// Parse input information
	flag.Parse()
	args := flag.Args()

	if len(args) < 4 {
		fmt.Printf("%s requires three arguments, a regex indicating the "+
			"target\ncatalogs, the output directory, the width of the "+
			"generated catalog\ngrid.", args[0])
		os.Exit(1)
	}

	matches := args[1 : len(args)-2]

	outPath := args[len(args)-2]
	gridWidth, err := strconv.Atoi(args[len(args)-1])

	if err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	} else if gridWidth <= 0 {
		fmt.Println("Grid width must be positive.")
		os.Exit(1)
	}

	snapDir := path.Base(path.Dir(matches[0]))
	paramDir := path.Base(path.Dir(path.Dir(path.Dir(matches[0]))))

	l, n, str, err := parseDir(paramDir)
	if err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}

	// Create files and rebin.

	outParamDir := fmt.Sprintf("Box_L%04d_N%04d_G%04d_%s", l, n, gridWidth, str)
	outParamPath := path.Join(outPath, outParamDir)
	outSnapPath := path.Join(outParamPath, snapDir)

	if err = os.MkdirAll(path.Join(outParamPath, snapDir), 0777); err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}

	catalogNames := createCatalogs(outSnapPath, matches[0], gridWidth, n)

	rebinParticles(matches, catalogNames, int64(gridWidth))
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

	ps := []tetra.Particle{}

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
// slice of tetra catalogs.
func rebinParticles(inFiles, outFiles []string, gridWidth int64) {
	hs := make([]tetra.Header, len(inFiles))
	for i := range hs {
		hs[i] = *catalog.ReadGadgetHeader(inFiles[i], gadgetEndianness)
	}

	width := hs[0].TotalWidth / float64(gridWidth)

	maxLen := int64(0)
	for _, h := range hs {
		if h.Count > maxLen {
			maxLen = h.Count
		}
	}

	floatBufMax := make([]float32, maxLen * 3)
	intBufMax := make([]int64, maxLen)
	psBufMax := make([]tetra.Particle, maxLen)

	bufs := createBuffers(outFiles, defaultBufSize)

	for i, inFile := range inFiles {
		fmt.Println(inFile)

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
	ps []tetra.Particle,
) {
	for _, p := range ps {
		xIdx := int64(math.Floor(p.Xs[0] / width))
		yIdx := int64(math.Floor(p.Xs[1] / width))
		zIdx := int64(math.Floor(p.Xs[2] / width))
		idx := xIdx + yIdx*gridWidth + zIdx*gridWidth*gridWidth
		bufs[idx].Append(p)
	}
}
