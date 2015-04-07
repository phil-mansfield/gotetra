package gotetra

import (
	"fmt"
	"log"
	"path"
	"runtime"

	"github.com/phil-mansfield/gotetra/density"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/io"
	"github.com/phil-mansfield/gotetra/rand"
)

const (
	UnitBufCount = 1 << 12
	tetraIntr = true
)

type Manager struct {
	// The currently loaded seet segment.
	xs, scaledXs []geom.Vec
	xCb geom.CellBounds
	hd io.SheetHeader

	renderers []renderer
	skip int
	unitBufs [][]geom.Vec

	// io related things
	log bool
	files []string
	ms runtime.MemStats

	// workspaces
	workers int
	workspaces []workspace
}

type renderer struct {
	box Box
	cb geom.CellBounds
	over Overlap
	validSegs map[string]bool
}

type workspace struct {
	buf []float64
	intr density.Interpolator
	lowX, highX int
}

func NewManager(files []string, boxes []Box, logFlag bool) (*Manager, error) {
	man := new(Manager)
	man.log = logFlag

	gen := rand.NewTimeSeed(rand.Tausworthe)
	man.unitBufs = make([][]geom.Vec, UnitBufCount)

	maxPoints := 0
	for _, b := range boxes {
		if b.Points() > maxPoints { maxPoints = b.Points() }
	}

	for bi := range man.unitBufs {
		man.unitBufs[bi] = make([]geom.Vec, maxPoints)
		buf := man.unitBufs[bi]
		for j := range buf {
			for k := 0; k < 3; k++ {
				buf[j][k] = float32(gen.Uniform(0, 1))
			}
		}
		geom.DistributeUnit(buf)
	}
	
	man.skip = 1

	man.workers = runtime.NumCPU()
	runtime.GOMAXPROCS(man.workers)
	man.workspaces = make([]workspace, man.workers)

	man.renderers = make([]renderer, len(boxes))
	for i := range man.renderers {
		man.renderers[i].box = boxes[i]
		man.renderers[i].cb = geom.CellBounds{
			boxes[i].CellOrigin(), boxes[i].CellSpan(),
		}
		man.renderers[i].validSegs = make(map[string]bool)
	}

	man.files = make([]string, 0)
	maxBufSize := 0

	for _, file := range files {
		err := io.ReadSheetHeaderAt(file, &man.hd)
		if err != nil { return nil, err }

		intersect := false
		for i := range boxes {
			cells := boxes[i].Cells()
			hCb := man.hd.CellBounds(cells)

			if hCb.Intersect(&man.renderers[i].cb, cells) {
				bufSize := boxes[i].Overlap(&man.hd).BufferSize()
				if bufSize > maxBufSize { maxBufSize = bufSize }
				
				man.renderers[i].validSegs[file] = true
				intersect = true
			}
		}

		if intersect {
			man.files = append(man.files, file)
		}
	}

	err := io.ReadSheetHeaderAt(files[0], &man.hd)
	if err != nil { return nil, err }
	man.xs = make([]geom.Vec, man.hd.GridCount)
	man.scaledXs = make([]geom.Vec, man.hd.GridCount)

	if man.log {
		log.Printf(
			"Workspace buffer size: %d. Number of workers: %d",
			maxBufSize, man.workers,
		)
	}

	for i := range man.workspaces {
		man.workspaces[i].buf = make([]float64, maxBufSize)
	}

	if man.log {
		runtime.ReadMemStats(&man.ms)
		log.Printf(
			"Alloc: %5d MB, Sys: %5d MB",
			man.ms.Alloc >> 20, man.ms.Sys >> 20,
		)
	}

	return man, nil
}

func (r *renderer) requiresFile(file string) bool {
	return r.validSegs[file]
}

func (r *renderer) scaleXs(man *Manager) {
	copy(man.scaledXs, man.xs)
	r.over.ScaleVecs(man.scaledXs)
}

func (r *renderer) ptVal(man *Manager) float64 {
	frac := float64(r.box.Cells()) / float64(man.hd.CountWidth)
	return frac * frac * frac
}

func (r *renderer) initWorkspaces(man *Manager) {
	segFrac := int(man.hd.SegmentWidth) / man.skip
	segLen := segFrac * segFrac * segFrac
	// chunkLen := segLen / man.workers

	for i := range man.unitBufs {
		man.unitBufs[i] = man.unitBufs[i][0: r.box.Points()]
	}

	// TODO: fix this int64 silliness
	for id := range man.workspaces {
		man.workspaces[id].intr = density.MonteCarlo(
			man.hd.SegmentWidth, r.box.Points(), r.box.Cells(),
			int64(man.skip), man.unitBufs, r.over,
		)

		man.workspaces[id].buf =
			man.workspaces[id].buf[0: r.over.BufferSize()]
		man.workspaces[id].lowX = id * man.skip
		man.workspaces[id].highX = segLen
		// man.workspaces[id].lowX = chunkLen * id
		// man.workspaces[id].highX = chunkLen * (id + 1)
		// if id == man.workers - 1 {
		//	man.workspaces[id].highX = segLen
		//}
	}
}

func (man *Manager) Log(flag bool) { man.log = flag }

func (man *Manager) Subsample(subsampleLength int) {
	if !isPowTwo(man.skip) {
		log.Fatalf("Skip ingrement is %d, must be power of two.", man.skip)
	}

	man.skip = subsampleLength
}

func isPowTwo(x int) bool {
	for x & 1 == 0 && x > 0 { x >>= 1 }
	return x == 1
}

func (man *Manager) RenderDensity() error {
	for _, file := range man.files {
		err := man.RenderDensityFromFile(file)
		if err != nil { return err }
	}
	return nil
}

func (man *Manager) RenderDensityFromFile(file string) error {
	if man.log {
		log.Printf("Rendering file %s", path.Base(file))
	}

	err := man.loadFile(file)
	if err != nil { return err }

	out := make(chan int, man.workers)

	for ri := range man.renderers {
		r := &man.renderers[ri]

		if !r.requiresFile(file) { continue }
		man.xCb = *man.hd.CellBounds(r.box.Cells())

		r.over = r.box.Overlap(&man.hd)
		r.cb = geom.CellBounds{ r.box.CellOrigin(), r.box.CellSpan() }
		r.scaleXs(man)
		r.initWorkspaces(man)

		for id := 0; id < man.workers - 1; id++ {
			go man.chanInterpolate(id, r, out)
		}
		id := man.workers - 1
		man.chanInterpolate(id, r, out)

		for i := 0; i < man.workers; i++ {
			id := <-out
			r.over.Add(man.workspaces[id].buf, r.box.Vals())
		}
	}
	
	if man.log {
		runtime.ReadMemStats(&man.ms)
		log.Printf(
			"Alloc: %5d MB, Sys: %5d MB",
			man.ms.Alloc >> 20, man.ms.Sys >> 20,
		)
	}
	return nil
}

func printBox(xs []float64, cb *geom.CellBounds) {
	idx := 0
	for z := 0; z < cb.Width[2]; z++ {
		for y := 0; y < cb.Width[1]; y++ {
			for x := 0; x < cb.Width[0]; x++ {
				fmt.Print(xs[idx])
				idx++
			}
			fmt.Println()
		}
		fmt.Println()
	}
}

func sum(xs []float64) float64 {
	s := 0.0
	for _, x := range xs { s += x }
	return s
}

func (man *Manager) loadFile(file string) error {
	err := io.ReadSheetHeaderAt(file, &man.hd)
	if err != nil { return err }
	err = io.ReadSheetPositionsAt(file, man.xs)
	if err != nil { return err }
	runtime.GC()

	return nil
}

func (man *Manager) chanInterpolate(id int, r *renderer, out chan<- int) {
	w := &man.workspaces[id]

	for i := range w.buf { w.buf[i] = 0.0 }
	if tetraIntr {
		w.intr.Interpolate(
			w.buf, man.scaledXs,
			r.ptVal(man), w.lowX, w.highX,
			man.workers * man.skip,
		)
	} else {
		r.over.Interpolate(
			w.buf, man.scaledXs,
			r.ptVal(man), w.lowX, w.highX,
			man.workers * man.skip,
		)
	}
	
	out <- id
}
