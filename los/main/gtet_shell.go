
package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/phil-mansfield/gotetra/los"
	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/los/analyze"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
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

	// Alternate modes
	MedianProfile bool
}

func main() {
	// Parse.
	log.Println("gtet_shell")
	p := parseCmd()
	ids, snaps, err := parseStdin()
	if err != nil { log.Fatal(err.Error()) }

	if len(ids) == 0 { return }

	// Compute coefficients.
	out := make([][]float64, len(ids))
	snapBins, idxBins := binBySnap(snaps, ids)
	buf := make([]analyze.RingBuffer, p.Rings)
	for i := range buf { buf[i].Init(p.Spokes, p.RBins) }

	sortedSnaps := []int{}
	for snap := range snapBins {
		sortedSnaps = append(sortedSnaps, snap)
	}
	sort.Ints(sortedSnaps)

	valids := make([]bool, len(ids))

	var losBuf *los.Buffers
	for _, snap := range sortedSnaps { 
		log.Println("Snap", snap)
		if snap == -1 { continue }
		snapIDs := snapBins[snap]
		idxs := idxBins[snap]

		// Bin halos
		hds, files, err := util.ReadHeaders(snap)
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
		
		if p.MedianProfile {
			// Calculate median profile.
			for i := range halos {
				out[idxs[i]] = calcMedian(&halos[i], p)
			}
		} else {
			// Calculate Penna coefficients.
			for i := range halos {
				var ok bool
				out[idxs[i]], ok = calcCoeffs(&halos[i], buf, p)
				if ok { valids[idxs[i]] = true }
			}
		}
		
	}
	
	ids = util.Filter(ids, valids)
	snaps = util.Filter(snaps, valids)
	printCoeffs(ids, snaps, out)
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
	flag.BoolVar(&p.MedianProfile, "MedianProfile", false,
		"Compute the median halo profile instead of the shell. " + 
			"KILL THIS OPTION.")
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

func createHalos(
	snap int, hd *io.SheetHeader, ids []int, p *Params,
) ([]los.HaloProfiles, error) {
	vals, err := util.ReadRockstar(
		snap, ids, halo.X, halo.Y, halo.Z, halo.Rad200b,
	)
	if err != nil { return nil, err }

	xs, ys, zs, rs := vals[0], vals[1], vals[2], vals[3]

	// Initialize halos.
	halos := make([]los.HaloProfiles, len(ids))
	seenIDs := make(map[int]bool)
	for i, id := range ids {
		origin := &geom.Vec{
			float32(xs[i]), float32(ys[i]), float32(zs[i]),
		}

		// If we've already seen a halo once, randomize its orientation.
		if seenIDs[id] {
			halos[i].Init(
				id, p.Rings, origin, rs[i] * p.MinMult, rs[i] * p.MaxMult,
				p.RBins, p.Spokes, hd.TotalWidth, los.Log(true),
				los.Rotate(float32(2 * math.Pi * rand.Float64()),
                    float32(2 * math.Pi * rand.Float64()),
                    float32(2 * math.Pi * rand.Float64())),
			)
		} else {
			seenIDs[id] = true
			halos[i].Init(
				id, p.Rings, origin, rs[i] * p.MinMult, rs[i] * p.MaxMult,
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
) ([]float64, bool) {
	for i := range buf {
		buf[i].Clear()
		buf[i].Splashback(halo, i, p.Window, p.Cutoff)
	}
	pxs, pys, ok := analyze.FilterPoints(buf, p.Levels)
	if !ok { return nil, false }
	cs, _ := analyze.PennaVolumeFit(pxs, pys, halo, p.Order, p.Order)
	return cs, true
}

func calcMedian(halo *los.HaloProfiles, p *Params) []float64 {
	rs := make([]float64, p.RBins)
	halo.GetRs(rs)
	rhos := halo.MedianProfile()
	return append(rs, rhos...)
}


func printCoeffs(ids, snaps []int, coeffs [][]float64) {
	idWidth, snapWidth := 0, 0
	coeffWidths := make([]int, len(coeffs[len(coeffs) - 1]))
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
