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

type Manager struct {
	hd io.SheetHeader
	xs, vs []geom.Vec
	cb *geom.CellBounds
	rhoBufs [][]float64
	unitBufs [][]geom.Vec
	intrs []density.Interpolator
	pts, skip, workers int
}
type Tetra struct { geom.Tetra }

const (
	MaxBufferLen = 1<<10
	UnitBufferLen = 1<<12
)

func (box *Box) InitFullFromFile(file string, cells int) {
	hd := io.SheetHeader{}
	io.ReadSheetHeaderAt(file, &hd)
	box.InitFull(hd.TotalWidth, cells)
}

func (box *Box) InitFull(boxWidth float64, cells int) {
	for i := range box.Vals { box.Vals[i] = 0.0 }
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

	box.cells = cells
	box.cellWidth = cellWidth
	box.full = false
}

func NewManager(files []string, cells, pts int) *Manager {
	man := new(Manager)

    max := 0
    for i, file := range files {
		io.ReadSheetHeaderAt(file, &man.hd)
		cb := man.hd.CellBounds(cells)

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

	return man
}

func (man *Manager) Load(file string) {
	man.LoadPositions(file)
	man.LoadVelocities(file)
}
func (man *Manager) LoadPositions(file string) {
	io.ReadSheetHeaderAt(file, &man.hd)
	io.ReadSheetPositionsAt(file, man.xs)
	runtime.GC()
}
func (man *Manager) LoadVelocities(file string) {
	io.ReadSheetHeaderAt(file, &man.hd)
	io.ReadSheetVelocitiesAt(file, man.vs)
	runtime.GC()
}

func (man *Manager) Skip(skip int) {
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

func (man *Manager) Tetra(x, y, z int, dir Direction) *Tetra {
	log.Fatalf("Tetra not yet implemented")
	return nil
}

func (man *Manager) TetraAt(x, y, z int, dir Direction, tet *Tetra) {
	log.Fatalf("TetraAt not yet implemented")
}

func (man *Manager) Density(rhos []Box) {
	if len(rhos) != 1 {
		log.Fatalf("For now you can only interpolate onto one box.")
	} else if !rhos[0].full {
		log.Fatalf("For now, you need to only interpolate onto full boxes.")
	}

	grid := rhos[0]

	out := make(chan int, man.workers)
	segLen := man.hd.SegmentWidth * man.hd.SegmentWidth * man.hd.SegmentWidth
	chunkLen :=  int(segLen) / man.workers

	man.cb = man.hd.CellBounds(grid.cells)
	man.cb.ScaleVecs(man.xs, grid.cells, man.hd.TotalWidth)

	for id := 0; id < man.workers - 1; id++ {
		low, high := chunkLen * id, chunkLen * (id + 1)
		go man.chanInterpolate(id, low, high, out)
	}
	id := man.workers - 1
	low, high := chunkLen * id, int(segLen)
	man.chanInterpolate(id, low, high, out)

	for i := 0; i < man.workers; i++ {
		id := <-out
		buf := man.rhoBufs[id]
		density.AddBuffer(grid.Vals, buf, man.cb, grid.cells)
	}
}

func (man *Manager) chanInterpolate(id, low, high int, out chan<- int) {
	buf := man.rhoBufs[id]
	intr := man.intrs[id]
	intr.Interpolate(buf, man.cb, 1.0, man.xs, low, high)
	out <- id
}

func (man *Manager) CurlDiv(curls [3][]Box, divs []Box) {
	log.Fatalf("CurlDiv not yet implemented")
}
