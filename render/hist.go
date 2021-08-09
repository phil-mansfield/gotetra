package render

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"runtime"
	"strings"

	"github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/io"
)

type HistInfo struct {
	Min, Max float64
	Bins int
	Scale string
}

type HistBox struct {
	Origin, Span [3]float64

	Centers []float64
	Counts []int
}

func NewHistBox(con *io.BoxConfig, bins int) HistBox {
	box := HistBox{
		Origin: [3]float64{con.X, con.Y, con.Z},
		Span: [3]float64{con.XWidth, con.YWidth, con.ZWidth},
		Centers: make([]float64, bins),
		Counts: make([]int, bins),
	}
	return box
}

func(box *HistBox) Contains(v geom.Vec, L float64) bool {
	for k := 0; k < 3; k++ {
		if !contains1D(box.Origin[k], box.Span[k], float64(v[k]), L) {
			return false
		}
	}
	return true
}

// HistManager is a struct which manages constructing 
type HistManager struct {
	xs []geom.Vec
	hd io.SheetHeader
	files []string

	boxes []HistBox
	unitBufs [][]geom.Vec

	// workspaces
	hists [][]int
	qs [][]float64
	inBox [][]bool
	vecBufs [][]geom.Vec
	
	skip int

	quantity string

	workers int

	gridHd *io.GridHeader
	grid []float64
}

// NewHistmanager creates a new HistManager.
func NewHistManager(
	files []string, boxes []HistBox, points int,
	quantity string, gridFile string,
) (*HistManager, error) {
	man := &HistManager{
		files: files,
		quantity: quantity,
		skip: 1,
	}
	err := io.ReadSheetHeaderAt(files[0], &man.hd)
	if err != nil { return nil, err }
	man.xs = make([]geom.Vec, man.hd.GridCount)

	// Create the unit cubes used to populate tetrahedra with points.
	man.unitBufs = unitBufs(UnitBufCount, points)

	// Set number of workers to number of available cores.
	man.workers = NumCores
	runtime.GOMAXPROCS(man.workers)

	man.gridHd, err = io.ReadGridHeader(gridFile)
	if err != nil { return nil, err }
	man.grid, err = io.ReadGrid(gridFile)
	if err != nil { return nil, err }

	man.boxes = boxes

	return man, nil
}

func (man *HistManager) Subsample(skip int) {
	if !isPowTwo(man.skip) {
		log.Fatalf("Skip ingrement is %d, must be power of two.", man.skip)
	}

	man.skip = skip
}

// Hist uses HistManager to compute a histogram with the given properties.
func (man *HistManager) Hist(info *HistInfo) error {
	// Set up workspaces.
	man.hists = make([][]int, man.workers)
	man.qs = make([][]float64, man.workers)
	man.inBox = make([][]bool, man.workers)
	man.vecBufs = make([][]geom.Vec, man.workers)

	pts := len(man.unitBufs[0])
	for i := 0; i < man.workers; i++ {
		man.hists[i] = make([]int, info.Bins)
		man.qs[i] = make([]float64, pts)
		man.inBox[i] = make([]bool, pts)
		man.vecBufs[i] = make([]geom.Vec, pts)
	}

	// Initialize Box output.
	for i := range man.boxes {
		man.boxes[i].Counts = make([]int, info.Bins)
		man.boxes[i].Centers = histCenters(info)
	}

	// Loop over files and do work.
	for i, file := range man.files {
		log.Printf("Analyzed files %d/%d", i, len(man.files))
		err := man.HistFromFile(file, info)
		if err != nil { return err }
	}
	return nil
}

// histCenters returns the centers of a histogram.
func histCenters(info *HistInfo) []float64 {
	min, max := info.Min, info.Max

	isLog := strings.ToLower(info.Scale) == "log"
	if isLog { min, max = math.Log10(min), math.Log10(max) }
	
	dx := (max - min) / float64(info.Bins)

	centers := make([]float64, info.Bins)
	for i := range centers {
		centers[i] = min + dx * (float64(i) + 0.5)
		if isLog { centers[i] = math.Pow(10, centers[i]) }
	}

	return centers
}

// HistFromFile updates the histograms of each box using only the particles in
// the given file.
func (man *HistManager) HistFromFile(file string, info *HistInfo) error {
	err := io.ReadSheetHeaderAt(file, &man.hd)
	if err != nil { return  err }
	out := make(chan int, man.workers)

	for bi := range man.boxes {
		if !histIntersect(&man.hd, &man.boxes[bi]) { continue }

		io.ReadSheetPositionsAt(file, man.xs)

		for id := 0; id < man.workers; id++ {
			go man.chanHistogram(id, &man.boxes[bi], info, out)
		}

		// merge worker histograms into the box histogram.
		for i := 0; i < man.workers; i++ {
			id := <-out
			for j := 0; j < info.Bins; j++ {
				man.boxes[bi].Counts[j] += man.hists[id][j]
			}			
		}

		sum := 0
		for i := range man.boxes[bi].Counts {
			sum += man.boxes[bi].Counts[i]
		}
	}
	return nil
}

// Returns true if sheet intersects with box and false otherwise.
func histIntersect(sheet *io.SheetHeader, box *HistBox) bool {
	for k := 0; k < 3; k++ {
		if !intersect1D(
			float64(sheet.Origin[k]), float64(sheet.Width[k]),
			box.Origin[k], box.Span[k], sheet.TotalWidth,
		) {
			return false
		}
	}
	return true
}

// intersect1D returns true if a segment starting at a with widht aSpan and 
// a second segment starting at b with width bSpan intersect within a periodic
// box of width L.
func intersect1D(a, aSpan, b, bSpan, L float64) bool {
	return contains1D(a, aSpan, b, L) ||
		contains1D(a, aSpan, b + bSpan, L) ||
		contains1D(b, bSpan, a, L) ||
		contains1D(b, bSpan, a + aSpan, L)
}

// contains1D returns true if b is within a segment that starts at a and has
// width aSpan within a periotic box of width L. a must be in [0, L), but b
// can be within [0, 2*L).
func contains1D(a, aSpan, b, L float64) bool {
	if b >= L { b -= L }
	return (b > a && a + aSpan > b) ||
		(b < a && a + aSpan - L > b)
}

// chanHistogram is a worker function run on a single thread which constructs
// a histogram using the buffers associated with the worker ID. This ID is sent
// to the out channel. 
func (man *HistManager) chanHistogram(
	worker int, box *HistBox, info *HistInfo, out chan<- int,
) {
	// Clear the histogram buffer.
	for i := range man.hists[worker] { man.hists[worker][i] = 0 }	

	// Convenience variables.
	gridWidth := int(man.hd.GridWidth)
	segWidth := int(man.hd.SegmentWidth)
	hist, qs, inBox := man.hists[worker], man.qs[worker], man.inBox[worker]
	vecBuf := man.vecBufs[worker]

	// Evaluate whatever quantity is being measured.
	switch strings.ToLower(man.quantity) {
	case strings.ToLower("Density"):
		for z := 0; z < segWidth; z += man.skip {
			for y := 0; y < segWidth; y += man.skip {
				for x := 0; x < segWidth; x += man.skip {
					idx := x + y*gridWidth + z*gridWidth*gridWidth
					if idx%(man.skip*man.skip*man.skip*man.workers) != worker {
						continue
					}

					if !man.cubeIntersects(idx) { continue }
					for dir := 0; dir < geom.TetraDirCount; dir++ {
						man.getDensities(idx, dir, box, vecBuf, qs, inBox)
						histogram(qs, inBox, info, hist)
					}
				}
			}
		}
	default:
		panic("Non-implemented quantity: " + man.quantity)
	}

	out <- worker
}

// tetIntersects returns true if the cube at idx intersects with the currently
// loaded sheet.
func (man *HistManager) cubeIntersects(idx int) bool {
	w := int(man.hd.GridWidth)
	ix := idx % w
	iy := (idx / w) % w
	iz := idx / (w*w)

	// Sanity check.
	if ix < 0 || man.skip + ix >= w ||
		iy < 0 || man.skip + iy >= w ||
		iz < 0 || man.skip + iz >= w {
		panic(fmt.Sprint("Internal inconsistency: attempting to analyze " + 
			"particle (%d %d %d), but SegmentWidth = %d", ix, iy, iz, w))
	}

	origin, span := man.cubeBoundingBox(ix, iy, iz)

	for k := 0; k < 3; k++ {
		if !intersect1D(
			float64(man.hd.Origin[k]), float64(man.hd.Width[k]),
			float64(origin[k]), float64(span[k]), man.hd.TotalWidth,
		) {
			return false
		}
	}
	return true
}

func (man *HistManager) cubeBoundingBox(
	ix, iy, iz int,
) (origin, span geom.Vec) {
	jump := man.skip
	w, L := int(man.hd.GridWidth), float32(man.hd.TotalWidth)
	origin, span = man.xs[ix + iy*w + iz*w*w], geom.Vec{ }

	for dz := 0; dz < jump; dz += jump {
		for dy := 0; dy < jump; dy += jump {
			for dx := 0; dx < jump; dx += jump {
				// Get the vector in each corner of the Lagrangian cube.
				x, y, z := ix + dx, iy + dy, iz + dz
				vec := man.xs[x + y*w + z*w*w]

				// Update origin and span in each dimension.
				for k := 0; k < 3; k++ {
					delta := vec[k] - origin[k]
					if delta > L/2 {
						delta -= L 
					} else if delta < -L/2 {
						delta += L 
					}

					if delta > 0 && delta > span[k] {
						span[k] = delta
					} else if delta < 0 {
						origin[k] -= delta
						if origin[k] < 0 { origin[k] += L }
						span[k] += delta
					}
				}
			}
		}
	}
	return origin, span
}

// getDensities computes the mass-weighted densities of the dir tetrahedron
// associated with the ith particle in the file. densities are written to the
// densities array and inBox specifies whether the Monte Carlo samples are
// inside box. If inBox[i] is false, densities[i] is set to -1.
func (man *HistManager) getDensities(
	i, dir int, box *HistBox,
	vecBuf []geom.Vec,
	densities []float64, inBox []bool,
) {
	// Precompute float32 conversions.
	L, origin, span := float32(man.hd.TotalWidth), geom.Vec{ }, geom.Vec{ }
	for k := 0; k < 3; k++ {
		origin[k], span[k] = float32(box.Origin[k]), float32(box.Span[k])
	}

	// Buffers
	tet, idxBuf := geom.Tetra{ }, geom.TetraIdxs{ }

	// Generate points from within the tetrahedra.
	idxBuf.Init(int64(i), man.hd.GridWidth, int64(man.skip), dir)
	tet.Init(
		&man.xs[idxBuf[0]], &man.xs[idxBuf[1]],
		&man.xs[idxBuf[2]], &man.xs[idxBuf[3]],
	)
	
	bufIdx := rand.Int63n(int64(len(man.unitBufs)))
	tet.DistributeTetra(man.unitBufs[bufIdx], vecBuf)

	for i := range vecBuf {
		// Put the vectors back into the periodic box.
		for k := 0; k < 3; k++ {
			if vecBuf[i][k] >= L { vecBuf[i][k] -= L }
			if vecBuf[i][k] < 0 { vecBuf[i][k] += L }
		}

		// Check if the point is in range, otherwise throw it out.
		if !box.Contains(vecBuf[i], man.hd.TotalWidth) {
			densities[i] = -1
			inBox[i] = false
			continue
		}

		// Finally, interpolate off the density grid.
		idx := gridIndex(man.gridHd, vecBuf[i])
		if idx >= 0 {
			densities[i] = man.grid[idx]
			//densities[i] = cloudInCell(man.grid, man.gridHd, idx, vecBuf[i])
			inBox[i] = true
		} else {
			densities[i] = -1
			inBox[i] = false
		}
	}
}

func gridIndex(hd *io.GridHeader, vec geom.Vec) int {
	idx := [3]int{ }
	L := hd.Cosmo.BoxWidth
	for k := 0; k < 3; k++ {
		delta := float64(vec[k]) - hd.Loc.Origin[k]
		if delta < 0 { delta += L }
		if delta >= hd.Loc.Span[k] { delta -= L }

		if delta < 0 || delta >= hd.Loc.Span[k] { return -1 }
		idx[k] = int(delta / hd.Loc.PixelWidth)
	}
	
	
	return idx[0] + idx[1]*int(hd.Loc.PixelSpan[0]) +
		idx[2]*int(hd.Loc.PixelSpan[0]*hd.Loc.PixelSpan[1])
}

func (man *HistManager) periodizeTetra(tet *geom.Tetra, L float32) {
	for i := 1; i < 4; i++ {
		for k := 0; k < 3; k++ {
			delta := tet.Corners[i][k] - tet.Corners[i][0]
			if delta > L/2 { tet.Corners[i][k] -= L }
			if delta < -L/2 { tet.Corners[i][k] += L }
		}
	}
} 

func histogram(x []float64, ok []bool, info *HistInfo, counts []int) {
	// Avoid unneeded dereferences and conversions
	min, max := info.Min, info.Max
	fBins := float64(info.Bins)

	if strings.ToLower(info.Scale) == "log" {
		min, max := math.Log10(min), math.Log10(max)
		dx := (max - min) / fBins

		for i := range x {
			if !ok[i] { continue }

			idx := (math.Log10(x[i]) - min) / dx
			if idx < 0 || idx >= fBins { continue }

			counts[int(idx)]++
			
		}
	} else {
		dx := (max - min) / fBins
		
		for i := range x {
			if !ok[i] { continue }

			idx := (x[i] - min) / dx
			if idx < 0 || idx >= fBins { continue }
			counts[int(idx)]++	
		}
	}
}
