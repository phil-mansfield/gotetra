package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"math/rand"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/phil-mansfield/gotetra/los"
	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/los/analyze"
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/render/halo"
)

type Params struct {
	// HaloProfiles params
	RBins, Spokes, Rings int
	MaxMult, MinMult float64

	// Splashback params
	Order, Window, Levels int
	Cutoff float64
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

	var losBuf *los.Buffers
	for snap, snapIDs := range snapBins {
		log.Println(snap)
		idxs := idxBins[snap]

		if snap == -1 { continue }

		// Bin halos
		hds, files, err := readHeaders(snap)
		if err != nil { err.Error() }
		if losBuf == nil { losBuf = los.NewBuffers(files[0], &hds[0]) }
		halos, err := createHalos(snap, &hds[0], snapIDs, p)
		if err != nil { log.Fatal(err.Error()) }
		intrBins := binIntersections(hds, halos)

		// Add densities. Done header by header to limit I/O time.
		for i := range hds {
			if len(intrBins[i]) == 0 { continue }
			hdContainer := []io.SheetHeader{hds[i]}
			fileContainer := []string{files[i]}
			los.LoadPtrDensities(
				intrBins[i], hdContainer, fileContainer, losBuf,
			)
		}
		
		// Calculate Penna coefficients.
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
	flag.IntVar(&p.Levels, "Levels", 4,
		"The number of recurve max-finding levels used by the 2D edge finder.")
	flag.Float64Var(&p.Cutoff, "Cutoff", 0.0,
		"The shallowest slope that can be considered a splashback point.")
	flag.Parse()
	return p
}

func parseStdin() (ids, snaps []int, err error) {
	ids, snaps = []int{}, []int{}
	lines, err := stdinLines()
	if err != nil { return nil, nil, err }
	for i, line := range lines {
		rawTokens := strings.Split(line, " ")
		tokens := make([]string, 0, len(rawTokens))
		for _, tok := range rawTokens {
			if len(tok) != 0 { tokens = append(tokens, tok) }
		}

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

	return ids, snaps, nil
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

func createHalos(
	snap int, hd *io.SheetHeader, ids []int, p *Params,
) ([]los.HaloProfiles, error) {
	// Read coordinates, radii and IDs of halos.
	rockstarDir := os.Getenv("GTET_ROCKSTAR_DIR")
	if rockstarDir == "" {
		return nil, fmt.Errorf("$GTET_ROCKSTAR_DIR not set.")
	} 

	files, err := dirContents(rockstarDir)
	if err != nil { return nil, err }
	// There are no halos in the 0th snapshot
	file := files[snap - 1]

	rids, vals, err := halo.ReadRockstarVals(
		file, &hd.Cosmo, halo.X, halo.Y, halo.Z, halo.M200b,
	)
	if err != nil { return nil, err }

	xs, ys, zs, ms := vals[0], vals[1], vals[2], vals[3]
	halo.R200m.Radius(&hd.Cosmo, ms, ms)
	rs := ms

	// Initialize halos.
	halos := make([]los.HaloProfiles, len(ids))
	seenIDs := make(map[int]bool)
	for i, id := range ids {
		idx := -1
		// TODO: Sort first if the ID count is large enough.
		for j, rid := range rids {
			if rid == id { 
				idx = j
				break
			}
		}

		if idx == -1 {
			return nil, fmt.Errorf("Halo ID %d not in halo catalog.", id)
		}
		
		origin := &geom.Vec{
			float32(xs[idx]), float32(ys[idx]), float32(zs[idx]),
		}

		// If we've already seen a halo once, randomize its orientation.
		if seenIDs[id] {
			halos[i].Init(
				id, p.Rings, origin, rs[idx] * p.MinMult, rs[idx] * p.MaxMult,
				p.RBins, p.Spokes, hd.TotalWidth, los.Log(true),
				los.Rotate(float32(2 * math.Pi * rand.Float64()),
                    float32(2 * math.Pi * rand.Float64()),
                    float32(2 * math.Pi * rand.Float64())),
			)
		} else {
			seenIDs[id] = true
			halos[i].Init(
				id, p.Rings, origin, rs[idx] * p.MinMult, rs[idx] * p.MaxMult,
				p.RBins, p.Spokes, hd.TotalWidth, los.Log(true),
			)
		}
	}

	return halos, nil
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

func calcCoeffs(
	halo *los.HaloProfiles, buf []analyze.RingBuffer, p *Params,
) []float64 {
	for i := range buf {
		buf[i].Clear()
		buf[i].Splashback(halo, i, p.Window, p.Cutoff)
	}
	pxs, pys := analyze.FilterPoints(buf, p.Levels)
	cs, _ := analyze.PennaVolumeFit(pxs, pys, halo, p.Order, p.Order)
	return cs
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

func dirContents(dir string) ([]string, error) {
	infos, err := ioutil.ReadDir(dir)
	if err != nil { return nil, err }
	
	files := make([]string, len(infos))
	for i := range infos {
		files[i] = path.Join(dir, infos[i].Name())
	}

	return files, nil
}
