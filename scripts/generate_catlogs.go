package main

import (
	"encoding/binary"
	"fmt"
	"io/ioutil"
	"os"
	"path"
	"strconv"
	"strings"

	tetra "github.com/phil-mansfield/gotetra"
	"github.com/phil-mansfield/gotetra/catalog"
)

func main() {
	// Parse input information

	if len(os.Args) != 4 {
		fmt.Printf("%s requires three arguments, a regex indicating the target"+
			" catalogs,\n the output directory, the width of the" + 
			" generated catalog grid.\n", os.Args[0])
		os.Exit(1)
	}

	pattern := os.Args[1]
	outPath := os.Args[2]
	gridWidth, err := strconv.Atoi(os.Args[3])
	if err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	} else if gridWidth <= 0 {
		fmt.Println("Grid width must be positive.")
		os.Exit(1)
	}

	snapDir := path.Base(path.Dir(pattern))
	paramDir := path.Base(path.Dir(path.Dir(pattern)))
	
	l, n, str, err := parseDir(paramDir)
	if err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}

	// Create requisite directories.

	outParamDir := fmt.Sprintf("Box_L%04d_N%04d_G%04d_%s", l, n, gridWidth,str)
	outParamPath := path.Join(outPath, outParamDir)
	if err = os.Mkdir(outParamPath, 0666); err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}

	outSnapPath := fmt.Sprintf(outParamPath, snapDir)
	if err = os.Mkdir(path.Join(outParamPath, snapDir), 0666); err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}	

	// Rebin catalogs.

	matches, err := allMatches(pattern)
	if err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}
	
	if err = createHeaders(outSnapPath, matches[0], gridWidth, n); err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}
}

// allMatches returns a slice of all filenames which path the given filepath.
// Standard unix regex syntax applies to the final element of the path.
func allMatches(pattern string) ([]string, error) {
	dir := path.Dir(pattern)
	matches := make([]string, 0)

	fs, err := ioutil.ReadDir(dir)
	if err != nil { return nil, err }

	for _, f := range fs {
		name := path.Join(dir, f.Name())
		if matched, err := path.Match(pattern, name); err != nil {
			return nil, err
		} else if matched {
			matches = append(matches, name)
		}
	}

	return matches, nil
}

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
		name := path.Join(outPath, fmt.Sprintf("gridcell_%04d.dat"))
		h.Idx = int64(i)
		fmt.Println(path.Base(name))
		fmt.Println(h)
		catalog.Write(name, h, ps)
	}

	return nil
}
