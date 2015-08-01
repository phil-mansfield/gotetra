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
}

func NewBuffers(file string, hd *io.SheetHeader) *Buffers {
	buf := new(Buffers)

    sw := hd.SegmentWidth
    buf.xs = make([]rGeom.Vec, hd.GridCount)
    buf.ts = make([]geom.Tetra, 6*sw*sw*sw)
    buf.ss = make([]geom.Sphere, 6*sw*sw*sw)
    buf.rhos = make([]float64, 6*sw*sw*sw)

	buf.Read(file, hd)
	return buf
}

func (buf *Buffers) Read(file string, hd *io.SheetHeader) {
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

	workers := runtime.NumCPU()
	runtime.GOMAXPROCS(workers)
	out := make(chan int, workers)
	for offset := 0; offset < workers - 1; offset++ {
		go buf.process(hd, offset, workers, out)
	}
	buf.process(hd, workers - 1, workers, out)

	for i := 0; i < workers; i++ { <- out }
}

func (buf *Buffers) process(
	hd *io.SheetHeader, offset, jump int, out chan<- int,
) {
	// Remember: Grid -> All particles; Segment -> Particles that can be turned
	// into tetrahedra.
	n := hd.SegmentWidth*hd.SegmentWidth*hd.SegmentWidth
	tw := hd.TotalWidth
	tFactor := tw*tw*tw / float64(hd.Count * 6)
	idxBuf := new(rGeom.TetraIdxs)
	jump64 := int64(jump)
	for segIdx := int64(offset); segIdx < n; segIdx+=jump64 {
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

	out <- offset
}

// CountAll computes profiles for all the given halos which count the number
// of tetrahedra overlapping points at a given radius.
func (buf *Buffers) CountAll(hs []HaloProfiles) {
	for hi := range hs {
		h := &hs[hi]
		for ti := range buf.ts {
			t := &buf.ts[ti]
			s := &buf.ss[ti]

			if h.Sphere.SphereIntersect(s) &&  
				!h.minSphere.TetraContain(t) {
				h.Count(t) 
			}
		}
	}
}

// DensityAll computes profiles for all the given halos which give the density
// of points at a given radius.
func (buf *Buffers) DensityAll(hs []HaloProfiles) {
	for hi := range hs {
		h := &hs[hi]
		for ti := range buf.ts {
			t := &buf.ts[ti]
			s := &buf.ss[ti]
			if h.Sphere.SphereIntersect(s) &&  
				!h.minSphere.TetraContain(t) {
				h.Density(t, buf.rhos[ti]) 
			}
		}
	}
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
func WrapHalo(hps []HaloProfiles, hd *io.SheetHeader) {
	tw := float32(hd.TotalWidth)
	newC := &geom.Vec{}
	for i := range hps {
		h := &hps[i]
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
