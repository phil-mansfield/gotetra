package main

import (
	"math"
	"path"
	"strings"
	"fmt"
	"time"
	"encoding/binary"
	"os"
	"log"
	"sort"
	"math/rand"

	"github.com/phil-mansfield/table"
	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/render/io"
)

const (
	finderCells = 150
	xCol = 17
	yCol = 18
	zCol = 19

	rType = halo.R200m
	mLow, mHigh = 8e11, 1e12 

	overlapMult = 3
	renderMult = 3
	pixels = 400

	// Remember: 7 arguments
	shellText = `#!/bin/sh
#SBATCH --job-name=tetra.Render_%s
#SBATCH --output=%s/Render_%s.out
#SBATCH --error=%s/Render_%s.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --time=2:00:00
#SBATCH --mem=29GB
#SBATCH --account=pi-kravtsov
module load go/1.3
go build -gcflags=-B github.com/phil-mansfield/gotetra && time ~/code/go/src/github.com/phil-mansfield/gotetra/main/main -Render %s %s`
	// Remember: 3 arguments
	renderText=`[Render]
Quantity = Density
Input  = %s
Output = %s
ImagePixels = %d
Particles = %d`
	// Remember: 5 arguments
	ballText=`[Ball "%s"]
X = %g
Y = %g
Z = %g
Radius = %g`
)

type Halos struct {
	xs, ys, zs, m200ms, r200ms, rPrints []float64
}

type Halo struct {
	ih int
	x, y, z, r float64
}

func (h *Halo) String() string {
	return fmt.Sprintf("%10d %13g %13g %13g %13g",
		h.ih, h.x, h.y, h.z, h.r)
}

func (hs *Halos) Len() int { return len(hs.r200ms) }
func (hs *Halos) Less(i, j int) bool { return hs.r200ms[i] < hs.r200ms[j] }
func (hs *Halos) Swap(i, j int) {
	hs.r200ms[i], hs.r200ms[j] = hs.r200ms[j], hs.r200ms[i]
	hs.m200ms[i], hs.m200ms[j] = hs.m200ms[j], hs.m200ms[i]
	hs.rPrints[i], hs.rPrints[j] = hs.rPrints[j], hs.rPrints[i]
	hs.xs[i], hs.xs[j] = hs.xs[j], hs.xs[i]
	hs.ys[i], hs.ys[j] = hs.ys[j], hs.ys[i]
	hs.zs[i], hs.zs[j] = hs.zs[j], hs.zs[i]
}

func main() {
	if len(os.Args) != 7 {
		log.Fatalf(
			"Usage: $ %s input_file catalog_file " + 
				"in_dir out_dir script_dir r_def", os.Args[0],
		)
	}

	fileName, catName := os.Args[1], os.Args[2]
	inDir, outDir, scriptDir := os.Args[3], os.Args[4], os.Args[5]
	rDef := os.Args[6]

	// Read catalog.

	hd := io.ReadGadgetHeader(catName, binary.LittleEndian)
	cosmo := &hd.Cosmo

	// Read Table.

	r200mCol := halo.R200m.RockstarColumn()
	rPrintType, ok := halo.RadiusFromString(rDef)
	if !ok { log.Fatalf("%s is not a recognized halo radius type.", rDef) }
	rPrintCol := rPrintType.RockstarColumn()

	colIdxs := []int{ xCol, yCol, zCol, r200mCol, rPrintCol }
	cols, err := table.ReadTable(fileName, colIdxs, nil)
	if err != nil { log.Fatalf(err.Error()) }

	xs, ys, zs, m200ms := cols[0], cols[1], cols[2], cols[3]
	r200ms := make([]float64, len(m200ms))
	rType.Radius(cosmo, m200ms, r200ms)

	// Find rPrint and mPrint using the correct conversion.

	var mPrints, rPrints []float64
	if rPrintType.RockstarMass() {
		mPrints := cols[4]
		rPrints = make([]float64, len(mPrints))
		rType.Radius(cosmo, mPrints, rPrints)
	} else {
		rPrints = cols[4]
		mPrints = make([]float64, len(rPrints))
		for i := range rPrints { rPrints[i] /= 1000 } // kpc -> Mpc
		rType.Mass(cosmo, rPrints, mPrints)
	}

	sort.Sort(sort.Reverse(&Halos{ xs, ys, zs, m200ms, r200ms, rPrints }))

	// Find Subhalos.

	g := halo.NewGrid(finderCells, hd.TotalWidth, len(xs))
	g.Insert(xs, ys, zs) 
	sf := halo.NewSubhaloFinder(g)
	sf.FindSubhalos(xs, ys, zs, r200ms, overlapMult)

	// Find some halos within a mass range that I'm interested in:

	iLow, iHigh := findIdx(m200ms, mLow), findIdx(m200ms, mHigh)
	rand.Seed(time.Now().UnixNano())
	
	idxs := []int{ 1051, 1092, 1115 }

	for i := 0; i < 10; i++ {
		ih := rand.Intn(iLow - iHigh) + iHigh
		if i >= 3 { break }
		ih = idxs[i]

		if sf.HostCount(ih) > 0 { continue }

		fmt.Printf("Halo %d  Mass: %.3g\n", ih, m200ms[ih])
		fmt.Printf(
			"Subhalos: %5d Hosts: %5d\n",
			sf.SubhaloCount(ih), sf.HostCount(ih),
		)

		name := fmt.Sprintf("%dh", ih)
		genScripts(
			name, inDir, outDir, scriptDir, hd.TotalWidth, rPrintType,
			newHalo(ih, xs, ys, zs, rPrints),
			newHalos(sf.Subhalos(ih), xs, ys, zs, rPrints),
		)
	}
}

func newHalo(ih int, xs, ys, zs, rs []float64) *Halo {
	return &Halo{ ih, xs[ih], ys[ih], zs[ih], rs[ih] }
}

func newHalos(ihs []int, xs, ys, zs, rs []float64) []Halo {
	hs := make([]Halo, len(ihs))
	for i, ih := range ihs {
		hs[i] = Halo{ ih, xs[ih], ys[ih], zs[ih], rs[ih] }
	}
	return hs
}

func particles(L, l float64, pixels int) int {
	n := int(500 * math.Pow(L / (3 * l), 3) *
		math.Pow(float64(pixels) / 500, 3))
	if n < 20 { n = 20 }
	return n
}

func genScripts(
	name, inDir, outDir, scriptDir string, L float64, 
	rType halo.Radius, h *Halo, intrs []Halo,
) {
	R := h.r * renderMult

	// Create names

	rs := rType.String()
	shellName := path.Join(scriptDir, fmt.Sprintf("%s_shell.sh", name))
	renderName := path.Join(scriptDir, fmt.Sprintf("%s_render.txt", name))
	ballName := path.Join(scriptDir, fmt.Sprintf("%s_ball.txt", name))
	shName := path.Join(scriptDir, fmt.Sprintf("%s_%s_sh.txt", name, rs))

	// Create bodies

	shellBody := fmt.Sprintf(
		shellText,
		name, scriptDir, name, scriptDir, name, renderName, ballName,
	)
	renderBody := fmt.Sprintf(
		renderText, inDir, outDir, pixels, particles(L, 2*R, pixels),
	)
	ballBody := fmt.Sprintf(ballText, name, h.x, h.y, h.z, R)

	shLines := make([]string, len(intrs) + 1)
	shLines[0] = h.String()
	for i, sh := range intrs { shLines[i + 1] = sh.String() }
	shBody := strings.Join(shLines, "\n")
	
	// Write

	write(shellName, shellBody)
	write(renderName, renderBody)
	write(ballName, ballBody)
	write(shName, shBody)
}

func write(name, body string) {
	f, err := os.Create(name)
	if err != nil { log.Fatal(err.Error()) }
	defer f.Close()
	_, err = f.Write(([]byte)(body))
	if err != nil { log.Fatal(err.Error()) }
}

func findIdx(xs []float64, x float64) int {
	for i, xx := range xs { if xx < x { return i } }
	return len(xs)
}
