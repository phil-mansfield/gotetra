package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/phil-mansfield/gotetra/los"
	"github.com/phil-mansfield/gotetra/los/analyze"
	"github.com/phil-mansfield/gotetra/render/io"
)

type Params struct {
	// HaloProfiles params
	RBins, Spokes, Rings int
	MaxMult, MinMult float64

	// Splashback params
	Order, Window int
}

func main() {
	// Parse.
	p := parseCmd()
	ids, snaps, err := parseStdin()
	if err != nil { log.Fatal(err.Error()) }

	if len(ids) == 0 { return }

	// Compute coefficients.
	coeffs := make([][]float64, len(ids))
	snapBins, idxBins := binBySnap(snaps, ids)
	buf := make([]analyze.RingBuffer, p.Rings)
	for i := range buf { buf[i].Init(p.Spokes, p.RBins) }

	for snap, snapIDs := range snapBins {
		idxs := idxBins[snap]

		if snap == -1 { continue }

		// Bin halos
		hds, err := readHeaders(snap)
		if err != nil { err.Error() }
		halos := createHalos(&hds[0], snapIDs, p)
		intrBins := binIntersections(hds, halos)

		// Add densities. Done header by header to limit I/O time.
		for i := range hds {
			addSheet(i, &hds[i], intrBins[i], p)
		}

		for i := range halos {
			coeffs[idxs[i]] = calcCoeffs(&halos[i], buf, p)
		}
	}

	printCoeffs(ids, snaps, coeffs)
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
	flag.IntVar(&p.Order, "Order", 5,
		"Order of the shell fitting function.")
	flag.IntVar(&p.Window, "Window", 121,
		"Number of bins within smoothign window. Must be odd.")
	flag.Parse()
	return p
}

func parseStdin() (ids, snaps []int, err error) {
	ids, snaps = []int{}, []int{}
	lines, err := stdinLines()
	if err != nil { return nil, nil, err }
	for i, line := range lines {
		tokens := strings.Split(line, " ")

		var (
			id, snap int
			err error
		)
		switch len(tokens) {
		case 0:
			continue
		case 2:
			id, err = strconv.Atoi(tokens[0])
			if err != nil {
				return nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[0],
				)
			} 
			snap, err = strconv.Atoi(tokens[1]) 
			if err != nil {
				return nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[1],
				)
			} 
		case 1:
			if tokens[0] == "" { continue }
			fallthrough
		default:
			return nil, nil, fmt.Errorf(
				"Line %d of stdin has %d tokens, but 2 are required.",
				i + 1, len(tokens),
			)
		}

		ids = append(ids, id)
		snaps = append(snaps, snap)
	}

	return ids, snaps, err
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

func readHeaders(snap int) ([]io.SheetHeader, error) {
	memoDir := os.Getenv("GTET_MEMO_DIR")
	if memoDir == "" {
		// You don't want to memoize? Fine. Deal with the consequences.
		return readHeadersFromSheet(snap)
	}
	if _, err := os.Stat(memoDir); err != nil {
		return nil, err
	}

	memoFile := path.Join(memoDir, fmt.Sprintf("hd_snap%d.dat", snap))

	if _, err := os.Stat(memoFile); err != nil {
		// File not written yet.
		hds, err := readHeadersFromSheet(snap)
		if err != nil { return nil, err }
		
        f, err := os.Create(memoFile)
        if err != nil { return nil, err }
        defer f.Close()
        binary.Write(f, binary.LittleEndian, hds)

		return hds, nil
	} else {
		// File exists: read from it instead.

		f, err := os.Open(memoFile)
        if err != nil { return nil, err }
        defer f.Close()
		
		n, err := sheetNum(snap)
		if err != nil { return nil, err }
		hds := make([]io.SheetHeader, n)
        binary.Read(f, binary.LittleEndian, hds) 
		return hds, nil
	}
	
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

func sheetNum(snap int) (int, error) {
	gtetFmt := os.Getenv("GTET_FMT")
	dir := fmt.Sprintf(gtetFmt, snap)
	files, err := dirContents(dir)
	if err != nil { return 0, err }
	return len(files), nil
}

func readHeadersFromSheet(snap int) ([]io.SheetHeader, error) {
	gtetFmt := os.Getenv("GTET_FMT")
	dir := fmt.Sprintf(gtetFmt, snap)
	files, err := dirContents(dir)
	hds := make([]io.SheetHeader, len(files))
	if err != nil { return nil, err }
	for i := range files {
		err = io.ReadSheetHeaderAt(files[i], &hds[i])
		if err != nil { return nil, err }
	}
	return hds, nil
}

func createHalos(hd *io.SheetHeader, ids []int, p *Params) []los.HaloProfiles {
	return nil
}

func binIntersections(
	hds []io.SheetHeader, halos []los.HaloProfiles,
) [][]*los.HaloProfiles {
	return nil
}

func addSheet(
	i int, hd *io.SheetHeader, halos []*los.HaloProfiles, p *Params,
) {
}

func calcCoeffs(
	halo *los.HaloProfiles, buf []analyze.RingBuffer, p *Params,
) []float64 {
	return []float64{1, 2, 3, 4, 5}
}

func printCoeffs(ids, snaps []int, coeffs [][]float64) {
	idWidth, snapWidth := 0, 0
	coeffWidths := make([]int, len(coeffs[0]))
	for i := range ids {
		iWidth := len(fmt.Sprintf("%d", ids[i]))
		sWidth := len(fmt.Sprintf("%d", snaps[i]))
		if iWidth > idWidth { idWidth = iWidth }
		if sWidth > snapWidth { snapWidth = sWidth }
		for j := range coeffs[i] {
			width := len(fmt.Sprintf("%.5g", coeffs[i][j]))
			if width > coeffWidths[j] { coeffWidths[j] = width }
		}
	}

	rowFmt := fmt.Sprintf("%%%dd %%%dd", idWidth, snapWidth)
	coefFmts := make([]string, len(coeffs[0]))
	for i := range coefFmts {
		rowFmt += fmt.Sprintf(" %%%d.5g", coeffWidths[i])
	}
	rowFmt += "\n"
	for i := range ids {
		args := append([]interface{}{ids[i], snaps[i]}, intr(coeffs[i])...)
		fmt.Printf(rowFmt, args...)
	}
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
