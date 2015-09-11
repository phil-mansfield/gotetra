package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"log"
	"io/ioutil"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/los/geom"
)

type Params struct {
	MaxMult float64
}

func main() {
	p := parseCmd()
	ids, snaps, coeffs, err := parseStdin()
	if err != nil { log.Fatal(err.Error()) }
	snapBins, coeffBins, idxBins := binBySnap(snaps, ids, coeffs)

	masses := make([]float64, len(ids))

	for snap, snapIDs := range snapBins {
		idxs := idxBins[snap]
		snapCoeffs := coeffBins[snap]
		if snap == -1 { continue }

		hds, files, err := readHeaders(snap)
		if err != nil { err.Error() }
		hBounds, err := boundingSpheres(snap, snapIDs, p)

		intrBins := binIntersections(hds, hBounds)

		for i := range hds {
			if len(intrBins[i]) == 0 { continue }
			for j := range idxs {
				masses[idxs[j]] += massContained(
					&hds[i], files[i], snapCoeffs[j], hBounds[j],
				)
			}
		}
	}

	printMasses(ids, snaps, masses)
}

func parseStdin() (ids, snaps []int, coeffs [][]float64, err error) {
	ids, snaps, coeffs = []int{}, []int{}, [][]float64{}
	lines, err := stdinLines()
	if err != nil { return nil, nil, nil, err }
	for i, line := range lines {
		rawTokens := strings.Split(line, " ")
		tokens := make([]string, 0, len(rawTokens))
		for _, tok := range rawTokens {
			if len(tok) != 0 { tokens = append(tokens, tok) }
		}

		var (
			id, snap int
			hCoeffs []float64
			err error
		)
		switch {
		case len(tokens) == 0:
			continue
		case len(tokens) <= 2:
			if tokens[0] == "" { continue }
			return nil, nil, nil, fmt.Errorf(
				"Line %d of stdin has 1 token, but >2 are required.", i + 1,
			)
		case len(tokens) > 2:
			id, err = strconv.Atoi(tokens[0])
			if err != nil {
				return nil, nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[0],
				)
			} 
			snap, err = strconv.Atoi(tokens[1]) 
			if err != nil {
				return nil, nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[1],
				)
			}
			
			hCoeffs = make([]float64, len(tokens) - 2) 
			for i := range hCoeffs {
				hCoeffs[i], err = strconv.ParseFloat(tokens[i + 2], 64)
			}
		}

		ids = append(ids, id)
		snaps = append(snaps, snap)
		coeffs = append(coeffs, hCoeffs)
	}

	return ids, snaps, coeffs, nil
}

func stdinLines() ([]string, error) {
	bs, err := ioutil.ReadAll(os.Stdin)
	if err != nil {
		return nil, fmt.Errorf(
			"Error reading stdin: %s.", err.Error(),
		)
	}

	text := string(bs)
	return strings.Split(text, "\n"), nil
}

func parseCmd() *Params {
	p := &Params{}
	flag.Float64Var(&p.MaxMult, "MaxMult", 3, 
		"Ending radius of LoSs as a multiple of R_200m. " + 
			"Should be the same value as used in gtet_shell.")
	flag.Parse()
	return p
}

func binBySnap(
	snaps, ids []int, coeffs [][]float64,
) (snapBins map[int][]int,coeffBins map[int][][]float64,idxBins map[int][]int) {
	snapBins = make(map[int][]int)
	coeffBins = make(map[int][][]float64)
	idxBins = make(map[int][]int)
	for i, snap := range snaps {
		snapBins[snap] = append(snapBins[snap], ids[i])
		coeffBins[snap] = append(coeffBins[snap], coeffs[i])
		idxBins[snap] = append(idxBins[snap], i)
	}
	return snapBins, coeffBins, idxBins
}

func readHeaders(snap int) ([]io.SheetHeader, []string, error) {
	memoDir := os.Getenv("GTET_MEMO_DIR")
	if memoDir == "" {
		// You don't want to memoize? Fine. Deal with the consequences.
		return readHeadersFromSheet(snap)
	}
	if _, err := os.Stat(memoDir); err != nil {
		return nil, nil, err
	}

	memoFile := path.Join(memoDir, fmt.Sprintf("hd_snap%d.dat", snap))

	if _, err := os.Stat(memoFile); err != nil {
		// File not written yet.
		hds, files, err := readHeadersFromSheet(snap)
		if err != nil { return nil, nil, err }
		
        f, err := os.Create(memoFile)
        if err != nil { return nil, nil, err }
        defer f.Close()
        binary.Write(f, binary.LittleEndian, hds)

		return hds, files, nil
	} else {
		// File exists: read from it instead.

		f, err := os.Open(memoFile)
        if err != nil { return nil, nil, err }
        defer f.Close()
		
		n, err := sheetNum(snap)
		if err != nil { return nil, nil, err }
		hds := make([]io.SheetHeader, n)
        binary.Read(f, binary.LittleEndian, hds) 

		gtetFmt := os.Getenv("GTET_FMT")
		dir := fmt.Sprintf(gtetFmt, snap)
		files, err := dirContents(dir)
		if err != nil { return nil, nil, err }

		return hds, files, nil
	}
}

func sheetNum(snap int) (int, error) {
	gtetFmt := os.Getenv("GTET_FMT")
	dir := fmt.Sprintf(gtetFmt, snap)
	files, err := dirContents(dir)
	if err != nil { return 0, err }
	return len(files), nil
}

func readHeadersFromSheet(snap int) ([]io.SheetHeader, []string, error) {
	gtetFmt := os.Getenv("GTET_FMT")
	dir := fmt.Sprintf(gtetFmt, snap)
	files, err := dirContents(dir)
	if err != nil { return nil, nil, err }

	hds := make([]io.SheetHeader, len(files))
	for i := range files {
		err = io.ReadSheetHeaderAt(files[i], &hds[i])
		if err != nil { return nil, nil, err }
	}
	return hds, files, nil
}

func dirContents(dir string) ([]string, error) {
	infos, err := ioutil.ReadDir(dir)
	if err != nil { return nil, err }
	
	files := make([]string, len(infos))
	for i := range infos {
		files[i] = path.Join(dir, infos[i].Name())
	}

	return files, nil
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

func inRange(x, r, low, width, tw float32) bool {
	return wrapDist(x, low, tw) > -r && wrapDist(x, low + width, tw) < r
}

// SheetIntersect returns true if the given halo and sheet intersect one another
// and false otherwise.
func sheetIntersect(s geom.Sphere, hd *io.SheetHeader) bool {
	tw := float32(hd.TotalWidth)
	return inRange(s.C[0], s.R, hd.Origin[0], hd.Width[0], tw) &&
		inRange(s.C[1], s.R, hd.Origin[1], hd.Width[1], tw) &&
		inRange(s.C[2], s.R, hd.Origin[2], hd.Width[2], tw)
}

func binIntersections(
	hds []io.SheetHeader, spheres []geom.Sphere,
) [][]geom.Sphere {
	bins := make([][]geom.Sphere, len(hds))
	for i := range hds {
		for si := range spheres {
			if sheetIntersect(spheres[si], &hds[i]) {
				bins[i] = append(bins[i], spheres[si])
			}
		}
	}
	return bins
}

func boundingSpheres(snap int, ids []int, p *Params) ([]geom.Sphere, error) {
	panic("NYI")
}

func massContained(
	hd *io.SheetHeader, file string, coeffs []float64, sphere geom.Sphere,
) float64 {
	panic("NYI")
}

func printMasses(ids, snaps []int, masses []float64) {
	idWidth, snapWidth, massWidth := 0, 0, 0
	for i := range ids {
		iWidth := len(fmt.Sprintf("%d", ids[i]))
		sWidth := len(fmt.Sprintf("%d", snaps[i]))
		mWidth := len(fmt.Sprintf("%.5g", masses[i]))
		if iWidth > idWidth { idWidth = iWidth }
		if sWidth > snapWidth { snapWidth = sWidth }
		if mWidth > massWidth { massWidth = mWidth }
	}

	rowFmt := fmt.Sprintf("%%%dd %%%dd %%%d.5g\n",
		idWidth, snapWidth, massWidth)
	for i := range ids { fmt.Printf(rowFmt, ids[i], snaps[i], masses[i]) }
}
