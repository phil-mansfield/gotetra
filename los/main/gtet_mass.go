package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"log"
	"io/ioutil"
	"math"
	"os"
	"path"
	"sort"
	"strconv"
	"strings"

	"github.com/phil-mansfield/gotetra/cosmo"
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/render/halo"
	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/los/geom"
	"github.com/phil-mansfield/gotetra/los/analyze"
)

const (
	minSnap = 30
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
	rads := make([]float64, len(ids))

	sortedSnaps := []int{}
	for snap := range snapBins {
		sortedSnaps = append(sortedSnaps, snap)
	}
	sort.Ints(sortedSnaps)

	log.Println("gtet_mass")
	for _, snap := range sortedSnaps {
		log.Println("Snap", snap)
		snapIDs := snapBins[snap]
		snapCoeffs := coeffBins[snap]
		idxs := idxBins[snap]

		if snap < minSnap { continue }

		hds, files, err := readHeaders(snap)
		if err != nil { log.Fatal(err.Error()) }
		hBounds, err := boundingSpheres(snap, &hds[0], snapIDs, p)
		if err != nil { log.Fatal(err.Error()) }
		intrBins := binIntersections(hds, hBounds)

		xs := []rgeom.Vec{}
		for i := range hds {
			if len(intrBins[i]) == 0 { continue }
			hd := &hds[i]

			n := hd.GridWidth*hd.GridWidth*hd.GridWidth
			if len(xs) == 0 { xs = make([]rgeom.Vec, n) }
			err := io.ReadSheetPositionsAt(files[i], xs)
			if err != nil { log.Fatal(err.Error()) }

			for j := range idxs {
				masses[idxs[j]] += massContained(
					&hds[i], xs, snapCoeffs[j], hBounds[j],
				)
			}
		}

		for j := range idxs {
			rads[idxs[j]] = rSp(&hds[0], snapCoeffs[j])
		}
	}

	printMasses(ids, snaps, masses, rads)
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

func boundingSpheres(
	snap int, hd *io.SheetHeader, ids []int, p *Params,
) ([]geom.Sphere, error) {
	rockstarDir := os.Getenv("GTET_ROCKSTAR_DIR")
	if rockstarDir == "" { 
		return nil, fmt.Errorf("$GTET_ROCKSTAR_DIR not set.")
	}
	
	hlists, err := dirContents(rockstarDir)
	if err != nil { return nil, err }
	rids, vals, err := halo.ReadRockstarVals(
		hlists[snap - 1], &hd.Cosmo, halo.X, halo.Y, halo.Z, halo.Rad200b,
	)
	xs, ys, zs, rs := vals[0], vals[1], vals[2], vals[3]

	spheres := make([]geom.Sphere, len(ids))
	for i := range spheres {
		j := -1
		for idx := range xs {
			if rids[idx] == ids[i] {
				j = idx
				break
			}
		}

		if j == -1 {
			return nil, fmt.Errorf("Halo %d not found in snap %d.",
				ids[i], snap)
		}
		spheres[i].C = geom.Vec{float32(xs[j]), float32(ys[j]), float32(zs[j])}
		spheres[i].R = float32(rs[j])
	}

	return spheres, nil
}

func findOrder(coeffs []float64) int {
	i := 1
	for {
		if 2*i*i == len(coeffs) {
			return i
		} else if 2*i*i > len(coeffs) {
			panic("Impossible")
		}
		i++
	}
}

func wrap(x, tw2 float32) float32 {
	if x > tw2 {
		return x - tw2
	} else if x < -tw2 {
		return x + tw2
	}
	return x
}

func coords(idx, cells int64) (x, y, z int64) {
    x = idx % cells
    y = (idx % (cells * cells)) / cells
    z = idx / (cells * cells)
    return x, y, z
}

func rSp(hd *io.SheetHeader, coeffs []float64) float64 {
	order := findOrder(coeffs)
	shell := analyze.PennaFunc(coeffs, order, order, 2)
	vol := shell.Volume(10 * 1000)

	r := math.Pow(vol / (math.Pi * 4 / 3), 0.33333)
	return r
}

func massContained(
	hd *io.SheetHeader, xs []rgeom.Vec, coeffs []float64, sphere geom.Sphere,
) float64 {
	c := &hd.Cosmo
	rhoM := cosmo.RhoAverage(c.H100 * 100, c.OmegaM, c.OmegaL, c.Z )
	dx := hd.TotalWidth / float64(hd.CountWidth) / (1 + c.Z)
	ptMass := rhoM * (dx*dx*dx)
	tw2 := float32(hd.TotalWidth) / 2

	order := findOrder(coeffs)
	shell := analyze.PennaFunc(coeffs, order, order, 2)

	// This prevents excess calls to the shell function:
	low, high := shell.RadialRange(10 * 1000)
	low2, high2 := float32(low*low), float32(high*high)

	sum := 0.0
	sw := hd.SegmentWidth
	for si := int64(0); si < sw*sw*sw; si++ {
		xi, yi, zi := coords(si, hd.SegmentWidth)
		i := xi + yi*sw + zi*sw*sw
		x, y, z := xs[i][0], xs[i][1], xs[i][2]
		x, y, z = x - sphere.C[0], y - sphere.C[1], z - sphere.C[2]
		x = wrap(x, tw2)
		y = wrap(y, tw2)
		z = wrap(z, tw2)

		r2 := x*x + y*y +z*z

		if r2 < low2 || ( r2 < high2 &&
			shell.Contains(float64(x), float64(y), float64(z))) {
			sum += ptMass
		}
	}
	return sum
}

func printMasses(ids, snaps []int, masses, rads []float64) {
	idWidth, snapWidth, massWidth, radWidth := 0, 0, 0, 0
	for i := range ids {
		iWidth := len(fmt.Sprintf("%d", ids[i]))
		sWidth := len(fmt.Sprintf("%d", snaps[i]))
		mWidth := len(fmt.Sprintf("%.5g", masses[i]))
		rWidth := len(fmt.Sprintf("%.5g", rads[i]))
		if iWidth > idWidth { idWidth = iWidth }
		if sWidth > snapWidth { snapWidth = sWidth }
		if mWidth > massWidth { massWidth = mWidth }
		if rWidth > radWidth { radWidth = rWidth }
	}

	rowFmt := fmt.Sprintf("%%%dd %%%dd %%%d.5g %%%d.5g\n",
		idWidth, snapWidth, massWidth, radWidth)
	for i := range ids { fmt.Printf(rowFmt, ids[i], snaps[i], masses[i], rads[i]) }
}
