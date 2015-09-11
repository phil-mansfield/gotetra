package los

import (
	"runtime"

	"github.com/phil-mansfield/gotetra/render/io"
	rGeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/los/geom"
)

type Buffers struct {
	xs []rGeom.Vec
	ts []geom.Tetra
	ss []geom.Sphere
	rhos []float64
	intr []bool
	bufHs []HaloProfiles
}

func NewBuffers(file string, hd *io.SheetHeader) *Buffers {
	buf := new(Buffers)

    sw := hd.SegmentWidth
    buf.xs = make([]rGeom.Vec, hd.GridCount)
    buf.ts = make([]geom.Tetra, 6*sw*sw*sw)
    buf.ss = make([]geom.Sphere, 6*sw*sw*sw)
    buf.rhos = make([]float64, 6*sw*sw*sw)
	buf.intr = make([]bool, 6*sw*sw*sw)

	buf.Read(file, hd)
	return buf
}

func (buf *Buffers) ParallelRead(file string, hd *io.SheetHeader) {
	workers := runtime.NumCPU()
	runtime.GOMAXPROCS(workers)
	buf.read(file, hd, workers)
}

func (buf *Buffers) Read(file string, hd *io.SheetHeader) {
	buf.read(file, hd, 1)
}

func (buf *Buffers) read(file string, hd *io.SheetHeader, workers int) {
	io.ReadSheetPositionsAt(file, buf.xs)
	tw := float32(hd.TotalWidth)
	// This can only be parallelized if we sychronize afterwards. This
	// is insignificant compared to the serial I/O time.
	for i := range buf.xs {
		for j := 0; j < 3; j++ {
			if buf.xs[i][j] < hd.Origin[j] {
				buf.xs[i][j] += tw
			}
		}
	}

	out := make(chan int, workers)
	for id := 0; id < workers - 1; id++ {
		go buf.chanRead(hd, id, workers, out)
	}
	buf.chanRead(hd, workers - 1, workers, out)

	for i := 0; i < workers; i++ { <- out }
}

func (buf *Buffers) chanRead(
	hd *io.SheetHeader, id, workers int, out chan<- int,
) {
	// Remember: Grid -> All particles; Segment -> Particles that can be turned
	// into tetrahedra.
	n := hd.SegmentWidth*hd.SegmentWidth*hd.SegmentWidth
	tw := hd.TotalWidth
	tFactor := tw*tw*tw / float64(hd.Count * 6)
	idxBuf := new(rGeom.TetraIdxs)

	jump := int64(workers)
	for segIdx := int64(id); segIdx < n; segIdx += jump {
		x, y, z := coords(segIdx, hd.SegmentWidth)
		for dir := int64(0); dir < 6; dir++ {
			ti := 6 * segIdx + dir
			idxBuf.InitCartesian(x, y, z, hd.GridWidth, int(dir))
			unpackTetra(idxBuf, buf.xs, &buf.ts[ti])
			buf.ts[ti].Orient(+1)

			buf.rhos[ti] = tFactor / buf.ts[ti].Volume()

			buf.ts[ti].BoundingSphere(&buf.ss[ti])
		}
	}

	out <- id
}

func (buf *Buffers) ParallelDensity(h *HaloProfiles) {
	workers := runtime.NumCPU()
	out := make(chan int, workers)

	for id := 0; id < workers - 1; id++ {
		go buf.chanIntersect(h, id, workers, out)
	}
	buf.chanIntersect(h, workers - 1, workers, out)
	for i := 0; i < workers; i++ { <-out }

	idxs, ok := splits(buf.intr, workers)
	if ! ok { return }

	if workers > len(h.rs) { workers = len(h.rs) }
	for id := 0; id < workers - 1; id++ {
		go buf.chanDensity(h, idxs[id], idxs[id+1], out)
	}
	buf.chanDensity(h, idxs[workers-1], idxs[workers], out)
	for i := 0; i < workers; i++ { <-out }
}

func (buf *Buffers) chanDensity(
	h *HaloProfiles, start, end int, out chan <- int,
) {
	for ri := start; ri < end; ri++ {
		r := &h.rs[ri]
		for ti := 0; ti < len(buf.ts); ti++ {
			if buf.intr[ti] { r.Density(&buf.ts[ti], buf.rhos[ti]) }
		}
	}

	out <- 0
}

func (buf *Buffers) chanIntersect(
	h *HaloProfiles, id, workers int, out chan <- int,
) {
	bufLen := len(buf.ts) / workers
	bufStart, bufEnd := id * bufLen, (id + 1) * bufLen
	if id == workers - 1 { bufEnd = len(buf.ts) }
	for i := bufStart; i < bufEnd; i++ {
		buf.intr[i] = h.Sphere.SphereIntersect(&buf.ss[i]) &&
			!h.minSphere.TetraContain(&buf.ts[i])
	}
	out <- id
}

func splits(intr []bool, workers int) (idxs []int, ok bool) {
	n := 0
	for _, ok := range intr {
		if ok { n++ }
	}

	if n == 0 { return nil, false }
	idxs = make([]int, workers + 1)

	spacing := n / workers
	// When m == spacing, insert a split at index j and reset m.
	m, j := 0, 1
	for i, ok := range intr {
		if ok {
			m++
			if m == spacing {
				idxs[j] = i + 1
				j++
				m = 0
				if j == workers { break }
			}
		}
	}

	for j = j; j <= workers ; j++ {
		idxs[j] = len(intr)
	}

	return idxs, true
}

func coords(idx, cells int64) (x, y, z int64) {
    x = idx % cells
    y = (idx % (cells * cells)) / cells
    z = idx / (cells * cells)
    return x, y, z
}

func index(x, y, z, cells int64) int64 {
    return x + y * cells + z * cells * cells
}

func unpackTetra(idxs *rGeom.TetraIdxs, xs []rGeom.Vec, t *geom.Tetra) {
    for i := 0; i < 4; i++ {
		t[i] = geom.Vec(xs[idxs[i]])
    }
}

// WrapHalo updates the coordinates of a slice of HaloProfiles so that they
// as close to the given sheet as periodic boundary conditions will allow.
func WrapHalo(hps []*HaloProfiles, hd *io.SheetHeader) {
	tw := float32(hd.TotalWidth)
	newC := &geom.Vec{}
	for i := range hps {
		h := hps[i]
		for j := 0; j < 3; j++ {
			if h.cCopy[j] + h.R < hd.Origin[j] {
				newC[j] = h.cCopy[j] + tw
			} else {
				newC[j] = h.cCopy[j]
			}
		}
		h.ChangeCenter(newC)
	}
}
