package main

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"path"
	"strconv"
	"strings"

	tetra "github.com/phil-mansfield/gotetra"
	"github.com/phil-mansfield/gotetra/catalog"
)

func main() {
	// Parse input information

	if len(os.Args) < 4 {
		fmt.Printf("%s requires three arguments, a regex indicating the " + 
			"target\ncatalogs, the output directory, the width of the" + 
			"generated catalog\ngrid.", os.Args[0])
		os.Exit(1)
	}

	matches := os.Args[1:len(os.Args) - 2]

	outPath := os.Args[len(os.Args) - 2]
	gridWidth, err := strconv.Atoi(os.Args[len(os.Args) - 1])

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

	// Create requisite directories.

	outParamDir := fmt.Sprintf("Box_L%04d_N%04d_G%04d_%s", l, n, gridWidth, str)
	outParamPath := path.Join(outPath, outParamDir)
	outSnapPath := path.Join(outParamPath, snapDir)
	if err = os.MkdirAll(path.Join(outParamPath, snapDir), 0777); err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}	

	// Rebin catalogs.
	
	if err = createHeaders(outSnapPath, matches[0], gridWidth, n); err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}

	width := float64(l) / float64(gridWidth)

	for _, match := range matches {
		_, ps := catalog.ReadGadget(match, binary.LittleEndian)
		bins := rebinParticles(width, int64(gridWidth), ps)
		for i, ps := range bins {
			if len(ps) > 0 { 
				name := path.Join(
					outSnapPath, fmt.Sprintf("gridcell_%04d.dat", i),
				)
				catalog.Append(name, ps)
			}
		}
	}
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
	if err != nil { return 0, 0, "", err }
	n, err := strconv.Atoi(parts[2][1:5])
	if err != nil { return 0, 0, "", err }

	return l, n, parts[3], nil
}

func dirErr(dir string) (int, int, string, error) {
	return 0, 0, "", fmt.Errorf("Invalid source directory '%s'.", dir)
}

func createHeaders(outPath string, exampleInput string, gridWidth int, countWidth int) error {
	h := catalog.ReadGadgetHeader(exampleInput, binary.LittleEndian)
	h.Width = h.TotalWidth / float64(gridWidth)
	h.CountWidth = int64(countWidth)
	h.GridWidth = int64(gridWidth)

	ps := make([]tetra.Particle, 0)

	maxIdx := gridWidth * gridWidth * gridWidth
	for i := 0; i < maxIdx; i++ {
		name := path.Join(outPath, fmt.Sprintf("gridcell_%04d.dat", i))
		h.Idx = int64(i)
		catalog.Write(name, h, ps)
	}

	return nil
}

func rebinParticles(width float64, gridWidth int64, ps []tetra.Particle) ([][]tetra.Particle) {
	totalCells := gridWidth * gridWidth * gridWidth
	bins := make([][]tetra.Particle, totalCells)

	for i := 0; i < len(bins); i++ { 
		bins[i] = make([]tetra.Particle, 0)
	}

	for _, p := range ps {
		xIdx := int64(math.Floor(p.Xs[0] / width))
		yIdx := int64(math.Floor(p.Xs[1] / width))
		zIdx := int64(math.Floor(p.Xs[2] / width))
		idx :=  xIdx + yIdx * gridWidth + zIdx * gridWidth * gridWidth
		bins[idx] = append(bins[idx], p)
	}

	return bins
}
