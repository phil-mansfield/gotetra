package gotetra

import (
	"log"
	"math"
	"runtime"

	"github.com/phil-mansfield/gotetra/density"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/io"
	"github.com/phil-mansfield/gotetra/rand"
)

type Direction int

// Box represents 
type Box struct {
	Vals []float64
	cb geom.CellBounds
	cells int
	cellWidth float64
	full bool
}

type baseManager struct {
	hd io.SheetHeader
	xs, vs []geom.Vec
	cb *geom.CellBounds
	rhoBufs [][]float64
	unitBufs [][]geom.Vec
	intrs []density.Interpolator
	pts, skip, workers int
	ptRho float64
}

type FullManager struct {
	baseManager
}

type BoundedManager struct {
	baseManager
	file string
	loaded bool
}

type Tetra struct { geom.Tetra }

type Manager interface {
	LoadPositions(file string)
	LoadVelocities(file string)
	Subsample(skip int)
	Density(boxes []Box)
	CurlDiv(curls [3][]Box, divs []Box)
}

// Type-checking
var (
	_ Manager = new(FullManager)
	_ Manager = new(BoundedManager)
)

const (
	MaxBufferLen = 1<<30
	UnitBufferLen = 1<<12
)

func (box *Box) InitFullFromFile(file string, cells int) {
	hd := io.SheetHeader{}
	io.ReadSheetHeaderAt(file, &hd)
	box.InitFull(hd.TotalWidth, cells)
}

func (box *Box) InitFull(boxWidth float64, cells int) {
	box.Vals = make([]float64, cells * cells * cells)
	box.cb.Origin = [3]int{0, 0, 0}
	box.cb.Width = [3]int{cells, cells, cells}
	box.cells = cells
	box.cellWidth = boxWidth / float64(cells)
	box.full = true
}

func (box *Box) Init(
	boxWidth float64, cells int,
	origin, width [3]float64,
) {
	cellWidth := boxWidth / float64(cells)

	for j := 0; j < 3; j++ {
		box.cb.Origin[j] = int(math.Floor(float64(origin[j]) / cellWidth))
		box.cb.Width[j] = 1 + int(math.Floor(
			float64(width[j] + origin[j]) / cellWidth),
		)
		box.cb.Width[j] -= box.cb.Origin[j]
	}

	box.Vals = make([]float64, box.cb.Width[0] * box.cb.Width[1] * box.cb.Width[2])

	box.cells = cells
	box.cellWidth = cellWidth
	box.full = false
}

func NewFullManager(files []string, maxCells, pts int) *FullManager {
	man := new(FullManager)

    max := 0
    for i, file := range files {
		io.ReadSheetHeaderAt(file, &man.hd)
		cb := man.hd.CellBounds(maxCells)

        vol := cb.Width[0] * cb.Width[1] * cb.Width[2]
        if vol > MaxBufferLen {
            log.Fatalf(
				"Header %d would have more than %d cells in it.",i,MaxBufferLen,
			)
		}

        if vol > max {
            max = vol
        }
    }

	man.init(files[0], max, pts)

	return man
}

func NewBoundedManager(files []string, boxes []Box, pts int) *BoundedManager {
	man := new(BoundedManager)

    max := 0
    for _, box := range boxes {
		vol := box.cb.Width[0] * box.cb.Width[1] * box.cb.Width[2]
		if vol > max { max = vol }
    }

	man.init(files[0], max, pts)
	man.loaded = false
	man.file = ""

	return man
}

func (man *baseManager) init(file string, max, pts int) {
	io.ReadSheetHeaderAt(file, &man.hd)

	gridLen := man.hd.GridWidth * man.hd.GridWidth * man.hd.GridWidth
	man.xs = make([]geom.Vec, gridLen)
	man.vs = make([]geom.Vec, gridLen)

	man.workers = runtime.NumCPU()
	man.rhoBufs = make([][]float64, man.workers)
	for i := range man.rhoBufs {
		man.rhoBufs[i] = make([]float64, max)
	}

	xs := make([]float64, pts)
	ys := make([]float64, pts)
	zs := make([]float64, pts)
	man.unitBufs = make([][]geom.Vec, UnitBufferLen)
	gen := rand.NewTimeSeed(rand.Golang)

	for i := range man.unitBufs {
		man.unitBufs[i] = make([]geom.Vec, pts)
		gen.UniformAt(0.0, 1.0, xs)
		gen.UniformAt(0.0, 1.0, ys)
		gen.UniformAt(0.0, 1.0, zs)
		geom.DistributeUnit(xs, ys, zs, man.unitBufs[i])
	}

	man.skip = 1
	man.pts = pts

	man.intrs = make([]density.Interpolator, man.workers)
	for i := range man.intrs {
		man.intrs[i] = density.MonteCarlo(
			man.hd.SegmentWidth, pts, int64(man.skip), man.unitBufs,
		)
	}

	runtime.GOMAXPROCS(man.workers)
}

func (man *FullManager) LoadPositions(file string) {
	io.ReadSheetHeaderAt(file, &man.hd)
	io.ReadSheetPositionsAt(file, man.xs)
	runtime.GC()
}

func (man *BoundedManager) LoadPositions(file string) {
	io.ReadSheetHeaderAt(file, &man.hd)

	// We might not actually have to read this.
	man.loaded = false
	man.file = file
}

func (man *baseManager) LoadVelocities(file string) {
	log.Fatalf("LoadVelocities is broken and breaks everything. Stahp.")
	io.ReadSheetHeaderAt(file, &man.hd)
	io.ReadSheetVelocitiesAt(file, man.vs)
	runtime.GC()
}

func (man *baseManager) Subsample(skip int) {
	if !isPowTwo(skip) {
		log.Fatalf("Skip ingrement is %d, must be power of two.", skip)
	}
	man.skip = skip

	for i := range man.intrs {
		man.intrs[i] = density.MonteCarlo(
			man.hd.SegmentWidth, man.pts, int64(man.skip), man.unitBufs,
		)
	}
}

func isPowTwo(x int) bool {
	for x & 1 == 0 && x > 0 { x >>= 1 }
	return x == 1
}

func (man *FullManager) Density(rhos []Box) {
	out := make(chan int, man.workers)
	segLen := man.hd.SegmentWidth * man.hd.SegmentWidth * man.hd.SegmentWidth
	chunkLen :=  int(segLen) / man.workers / (man.skip * man.skip * man.skip)

	for _, grid := range rhos {
		if !grid.full {
			log.Fatalf("Supplied a non-full Box to a FullManager.")
		}

		man.cb = man.hd.CellBounds(grid.cells)

		frac := float64(grid.cells) / float64(man.hd.CountWidth)
		man.ptRho = frac * frac * frac
		
		for id := 0; id < man.workers - 1; id++ {
			low, high := chunkLen * id, chunkLen * (id + 1)
			go man.chanInterpolate(id, low, high, grid.cells, out)
		}
		
		id := man.workers - 1
		low, high := chunkLen * id, int(segLen)
		man.chanInterpolate(id, low, high, grid.cells, out)
		
		for i := 0; i < man.workers; i++ {
			id := <-out
			buf := man.rhoBufs[id]
			density.AddBuffer(grid.Vals, buf, man.cb, grid.cells)
		}
	}
}

func (man *FullManager) chanInterpolate(id, low, high, cells int, out chan<- int) {
	skipVol := man.skip * man.skip * man.skip
	man.cb.ScaleVecs(man.xs[low * skipVol: high * skipVol], cells, man.hd.TotalWidth)

	buf := man.rhoBufs[id]
	for i := range buf { buf[i] = 0.0 }
	intr := man.intrs[id]
	intr.Interpolate(buf, man.cb, man.ptRho, man.xs, low, high)
	out <- id
}

func (man *BoundedManager) Density(rhos []Box) {
	out := make(chan int, man.workers)
	segLen := man.hd.SegmentWidth * man.hd.SegmentWidth * man.hd.SegmentWidth
	chunkLen :=  int(segLen) / man.workers / (man.skip * man.skip * man.skip)

	for _, grid := range rhos {
		if grid.full {
			log.Fatalf("Supplied a full Box to a BoundedManager.")
		}

		man.cb = man.hd.CellBounds(grid.cells)
		if !man.loaded && man.cb.Intersect(&grid.cb, grid.cells) {
			if man.file == "" {
				log.Fatalf("No file has been loaded prior to ")
			}
			io.ReadSheetPositionsAt(man.file, man.xs)
			man.loaded = true
		}

		frac := float64(grid.cells) / float64(man.hd.CountWidth)
		man.ptRho = frac * frac * frac
		
		for id := 0; id < man.workers - 1; id++ {
			low, high := chunkLen * id, chunkLen * (id + 1)
			go man.chanBoundedInterpolate(id, low, high, grid.cells, &grid.cb, out)
		}
		
		id := man.workers - 1
		low, high := chunkLen * id, int(segLen)
		man.chanBoundedInterpolate(id, low, high, grid.cells, &grid.cb, out)
		
		bufLen := grid.cb.Width[0] * grid.cb.Width[1] * grid.cb.Width[2]
		for i := 0; i < man.workers; i++ {
			id := <-out
			buf := man.rhoBufs[id]
			for j := 0; j < bufLen; j++ { grid.Vals[j] += buf[j] }
		}
	}
}

func (man *BoundedManager) chanBoundedInterpolate(
	id, low, high, cells int, boxCb *geom.CellBounds, out chan<- int,
) {

	skipVol := man.skip * man.skip * man.skip
	man.cb.ScaleVecs(man.xs[low * skipVol: high * skipVol], cells, man.hd.TotalWidth)

	buf := man.rhoBufs[id]
	for i := range buf { buf[i] = 0.0 }
	intr := man.intrs[id]
	intr.BoundedInterpolate(buf, boxCb, cells, man.ptRho, man.xs, man.cb, low, high)
	out <- id
}

func (man *baseManager) CurlDiv(curls [3][]Box, divs []Box) {
	log.Fatalf("CurlDiv not yet implemented")
}
