package main

import (
	"fmt"
	goio "io"
	"log"
	"math"
	"math/rand"
	"os"
	"time"

	"github.com/phil-mansfield/table"
	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/render/io"
)

type Halo struct {
	ih int
	pos [3]float64
	r float64
}

func main() {
	if len(os.Args) != 4 {
		log.Fatal(
			"Required file use: $ %s gtet_file sh_file out_prefix", os.Args[0],
		)
	}

	gtetFile, shFile, outPrefix := os.Args[1], os.Args[2], os.Args[3]

	_, shs := readSubhalos(shFile)

	hd, err := io.ReadGridHeader(gtetFile)
	if err != nil { log.Fatal(err.Error()) }
	grid, err := io.ReadGrid(gtetFile)
	if err != nil { log.Fatal(err.Error()) }

	bh, bsh := &halo.Bounds{}, &halo.Bounds{}
	for i := 0; i < 3; i++ {
		bh.Origin[i] = int(hd.Loc.PixelOrigin[i])
		bh.Span[i] = int(hd.Loc.PixelSpan[i])
	}
	cw := hd.Loc.PixelWidth
	width := hd.Cosmo.BoxWidth
	totalPixels := round(hd.Cosmo.BoxWidth / cw)

	for _, sh := range shs {
		bsh.SphereBounds(sh.pos, sh.r, cw, width)
		removeSubhalo(bh, bsh, grid, totalPixels)
	}

	rs := radii(bh, cw)
	valid := make([]bool, len(rs))
	for i := range valid { valid[i] = true }

	rand.Seed(time.Now().UnixNano())
	vecs, rProfs, rhoProfs, _ := lineProfiles(grid, bh, 1000, cw, -1)
	xVecs, xRProfs, xRhoProfs, _ := lineProfiles(grid, bh, 1000, cw, 0)
	yVecs, yRProfs, yRhoProfs, _ := lineProfiles(grid, bh, 1000, cw, 1)
	zVecs, zRProfs, zRhoProfs, _ := lineProfiles(grid, bh, 1000, cw, 2)
	mainOut := fmt.Sprintf("%s_lines.txt", outPrefix)
	xOut := fmt.Sprintf("%s_x_lines.txt", outPrefix)
	yOut := fmt.Sprintf("%s_y_lines.txt", outPrefix)
	zOut := fmt.Sprintf("%s_z_lines.txt", outPrefix)
	printData(mainOut, vecs, rProfs, rhoProfs)
	printData(xOut, xVecs, xRProfs, xRhoProfs)
	printData(yOut, yVecs, yRProfs, yRhoProfs)
	printData(zOut, zVecs, zRProfs, zRhoProfs)
}

func printData(fname string, vecs [][3]float64, rProfs, rhoProfs [][]float64) {
	f, err := os.Create(fname)
	if err != nil { log.Fatal(err.Error()) }
	fmt.Fprintln(f, "# The first three rows are unit vector components.")
	printVecs(f, vecs)
	fmt.Fprintln(f, "# The rest are individual profiles.")
	printProfiles(f, rProfs, rhoProfs)
}

func lineProfiles(
	grid []float64, bh *halo.Bounds, profCount int, cw float64, axis int,
) (vecs [][3]float64, rProfs, rhoProfs [][]float64, attempts int) {
	vecs = make([][3]float64, 0, profCount)
	rProfs = make([][]float64, 0, profCount)
	rhoProfs = make([][]float64, 0, profCount)

	attempts = 0
	for len(rProfs) != profCount {
		unit := randomUnit(axis)
		rs, rhos, ok := extractLine(grid, bh, unit, cw)
		attempts++
		if !ok { continue }
		rProfs = append(rProfs, rs)
		rhoProfs = append(rhoProfs, rhos)
		vecs = append(vecs, unit)
	}

	return vecs, rProfs, rhoProfs, attempts
}


func printVecs(f goio.Writer, vecs [][3]float64) {
	for dim := 0; dim < 3; dim++ {
		for _, vec := range vecs {
			fmt.Fprintf(f, "%9.4g %9.4g ", vec[dim], math.NaN())
		}
		fmt.Fprintln(f)
	}
}

func printProfiles(f goio.Writer, rProfs, rhoProfs [][]float64) {
	minLen := len(rProfs[0])
	for _, prof := range rProfs {
		if minLen > len(prof) { minLen = len(prof) }
	}

	for col := 0; col < minLen; col++ {
		for profNum := range rProfs {
			r, rho := rProfs[profNum][col], rhoProfs[profNum][col]
			fmt.Fprintf(f, "%9.4g %9.4g ", r, rho)
		}
		fmt.Fprintln(f)
	}
}

func randomUnit(axis int) [3]float64 {
	x := rand.Float64() * 2 - 1
	y := rand.Float64() * 2 - 1
	z := rand.Float64() * 2 - 1

	switch axis {
	case 0:
		x = 0
	case 1:
		y = 0
	case 2:
		z = 0
	}

	norm := math.Sqrt(x*x + y*y + z*z)
	vec := [3]float64{}
	vec[0], vec[1], vec[2] = x/norm, y/norm, z/norm
	return vec
}

func gridMult(b *halo.Bounds, axis int) int {
	switch axis {
	case 0:
		return 1
	case 1:
		return b.Span[0]
	case 2:
		return b.Span[0] * b.Span[1]
	}
	panic(":3")
}

func extractLine(
	grid []float64, b *halo.Bounds, unit [3]float64, cw float64,
) (rs, rhos []float64, ok bool) {
	// Assumes a cubical box. We'll do some checks to make sure
	// this doens't break anything.
	axis0 := 0 // Major axis.
	if math.Abs(unit[1]) > math.Abs(unit[axis0]) { axis0 = 1 }
	if math.Abs(unit[2]) > math.Abs(unit[axis0]) { axis0 = 2 }
	axis1, axis2 := (axis0 + 1) % 3, (axis0 + 2) % 3
	mult0 := gridMult(b, axis0)
	mult1 := gridMult(b, axis1)
	mult2 := gridMult(b, axis2)

	rs = make([]float64, b.Span[axis0] / 2)
	rhos = make([]float64, b.Span[axis0] / 2)

	// Centers and conversion factors.
	c0 := float64(b.Span[0]) / 2
	c1 := float64(b.Span[1]) / 2
	c2 := float64(b.Span[2]) / 2
	dx01 := unit[axis1] / unit[axis0]
	dx02 := unit[axis2] / unit[axis0]

	n := 0
	smoothStarted := false
	for j := 0; j < len(rs); j++ {
		// Points are spaced by cell on the major axis and are interpolated
		// using CIC on the the other two axes.
		dxp0 := float64(j)
		if unit[axis0] < 0 { dxp0 = -dxp0 }
		dxp1, dxp2 := dxp0 * dx01, dxp0 * dx02

		xp0, xp1, xp2 := c0 + dxp0, c1 + dxp1, c2 + dxp2
		i0, i1, i2 := int(xp0), int(xp1), int(xp2)
		if i0 < 0 || i1 < 0 || i2 < 0 || i0 + 1 >= b.Span[axis0] || 
			i1 + 1 >= b.Span[axis1] || i2 + 1 >= b.Span[axis2] {
			break
		}
		xc1, xc2 := float64(i1), float64(i2)
		d1, d2 := xp1 - xc1, xp2 - xc2
		t1, t2 := 1 - d1, 1 - d2

		rho := grid[i0*mult0 + i1*mult1 + i2*mult2]*t1*t2 +
			grid[i0*mult0 + i1*mult1 + (i2 + 1)*mult2]*t1*d2 +
			grid[i0*mult0 + (i1 + 1)*mult1 + i2*mult2]*d1*t2 +
			grid[i0*mult0 + (i1 + 1)*mult1 + (i2 + 1)*mult2]*d1*d2
		if math.IsNaN(rho) {
			if smoothStarted {
				return nil, nil, false
			} else {
				//smoothStarted = true
			}
		}
		rhos[j], rs[j] = rho, math.Sqrt(dxp0*dxp0 + dxp1*dxp1 + dxp2*dxp2) * cw
		n++
	}

	return rs[:n], rhos[:n], true
}

func round(x float64) int {
	if math.Abs(math.Floor(x) - x) < 0.5 {
		return int(math.Floor(x))
	} else {
		return int(math.Ceil(x))
	}
}

func radii(b *halo.Bounds, cw float64) []float64 {
	n := b.Span[0] * b.Span[1] * b.Span[2]
	rs := make([]float64, n)

	// This is a slightly suboptimal way to write this.
	midX := (float64(b.Span[0]) - 0.5) * cw / 2 
	midY := (float64(b.Span[1]) - 0.5) * cw / 2 
	midZ := (float64(b.Span[2]) - 0.5) * cw / 2 

	idx := 0
	for z := 0; z < b.Span[2]; z++ {
		dz := float64(z) * cw - midZ
		dz2 := dz*dz
		for y := 0; y < b.Span[1]; y++ {
			dy := float64(y) * cw - midY
			dy2 := dy*dy
			for x := 0; x < b.Span[0]; x++ {
				dx := float64(x) * cw - midX
				dx2 := dx*dx

				rs[idx] = math.Sqrt(dx2 + dy2 + dz2)
				idx++
			}
		}
	}

	return rs
}

// This is slow, but I don't particularly care.
func removeSubhalo(bh, bsh *halo.Bounds, grid []float64, width int) {
	loX, hiX := bsh.Origin[0], bsh.Origin[0] + bsh.Span[0]
	loY, hiY := bsh.Origin[1], bsh.Origin[1] + bsh.Span[1]
	loZ, hiZ := bsh.Origin[2], bsh.Origin[2] + bsh.Span[2]
	nan := math.NaN()

	n := 0

	for z := loZ; z < hiZ; z++ {
		if !bh.Inside(z, width, 2) { continue }
		for y := loY; y < hiY; y++ {
			if !bh.Inside(y, width, 1) { continue }
			for x := loX; x < hiX; x++ {
				if !bh.Inside(x, width, 0) { continue }
				bx, by, bz := bh.ConvertIndices(x, y, z, width)
				grid[bz*bh.Span[0]*bh.Span[1] + by*bh.Span[0] + bx] = nan
				n++
			}
		}
	}
}

type binnedMeanOptions struct {
	valid []bool
	lowLim, highLim float64
	useLog bool
}

var defaultBinnedMeanOptions = binnedMeanOptions{
	valid: nil,
	lowLim: math.NaN(),
	highLim: math.NaN(),
	useLog: false,
}

func binnedMean(
	xs, ys []float64, binCount int, opt *binnedMeanOptions,
) (vals, means []float64) {
	if opt == nil { opt = &defaultBinnedMeanOptions }

	var maxX, minX float64
	if math.IsNaN(opt.lowLim) {
		minX = math.MaxFloat64
		for i, x := range xs {
			if opt.valid != nil && !opt.valid[i] { continue }
			if math.IsNaN(ys[i]) { continue }
			if x < minX { minX = x }
		}
	} else {
		minX = opt.lowLim
	}

	if math.IsNaN(opt.lowLim) {
		maxX = math.SmallestNonzeroFloat64
		for i, x := range xs {
			if opt.valid != nil && !opt.valid[i] { continue }
			if math.IsNaN(ys[i]) { continue }
			if x > maxX { maxX = x }
		}
	} else {
		maxX = opt.highLim
	}

	if opt.useLog { minX, maxX = math.Log(minX), math.Log(maxX) }
	binWidth := (maxX - minX) / float64(binCount)

	counts, means := make([]int, binCount), make([]float64, binCount)
	vals = make([]float64, binCount)

	nanCounts := make([]int, binCount)
	for i := 0; i < len(xs); i++ {
		if opt.valid != nil && !opt.valid[i] { continue }

		var idx int
		if opt.useLog {
			lx := math.Log(xs[i])
			if lx < minX || lx > maxX { continue }
			idx = int((lx - minX) / binWidth)
		} else {
			if xs[i] < minX || xs[i] > maxX { continue }
			idx = int((xs[i] - minX) / binWidth)
		}

		if idx == binCount { idx-- }
		if math.IsNaN(ys[i]) {
			nanCounts[idx]++
			continue 
		}
		counts[idx]++
		means[idx] += ys[i]
	}

	sum := 0
	for _, n := range counts { sum += n }

	shiftIdx := 0
	for idx := range counts {
		if counts[idx] == 0 { continue }

		means[shiftIdx] = means[idx] / float64(counts[idx])
		vals[shiftIdx] = (float64(idx) + 0.5) * binWidth + minX
		if opt.useLog { vals[shiftIdx] = math.Exp(vals[shiftIdx]) }
		shiftIdx++
	}

	return vals[:shiftIdx], means[:shiftIdx]
}

func readSubhalos(fname string) (h *Halo, shs []Halo) {
	cols, err := table.ReadTable(fname, []int{0, 1, 2, 3, 4}, nil)
	if err != nil { log.Fatal(err.Error()) }

	ihs, xs, ys, zs, rs := cols[0], cols[1], cols[2], cols[3], cols[4]
	h = &Halo{ int(ihs[0]), [3]float64{xs[0], ys[0], zs[0]}, rs[0] }
	shs = make([]Halo, len(ihs) - 1)
	for i := range shs {
		shs[i] = Halo{
			int(ihs[i+1]), [3]float64{xs[i+1], ys[i+1], zs[i+1]}, rs[i+1],
		}
	}

	return h, shs
}
