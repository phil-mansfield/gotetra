package los

import (
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
	// This can only be parallelized if we sychronize afterwards. Might not be
	// worth it.
	for i := range buf.xs {
		for j := 0; j < 3; j++ {
			if buf.xs[i][j] < hd.Origin[j] {
				buf.xs[i][j] += tw
			}
		}
	}

	// Remember: Grid -> All particles; Segment -> Particles that can be turned
	// into tetrahedra.
	n := hd.SegmentWidth*hd.SegmentWidth
	tFactor := float64(tw*tw*tw) / float64(n * 6)
	idxBuf := new(rGeom.TetraIdxs)
	for segIdx := int64(0); segIdx < n; segIdx++ {
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

// UnpackTetrahedra converts the raw position data in a sheet segment into
// tetrahedra.
func unpackTetrahedra(
	xs []rGeom.Vec, hd *io.SheetHeader, tsBuf []geom.Tetra,
) {
	n := hd.SegmentWidth*hd.SegmentWidth*hd.SegmentWidth
	idxBuf := new(rGeom.TetraIdxs)
	for writeIdx := int64(0); writeIdx < n; writeIdx++ {
		x, y, z := coords(writeIdx, hd.SegmentWidth)
		for dir := int64(0); dir < 6; dir++ {
			tIdx := 6 * writeIdx + dir
			idxBuf.InitCartesian(x, y, z, hd.GridWidth, int(dir))
			unpackTetra(idxBuf, xs, &tsBuf[tIdx])
			tsBuf[tIdx].Orient(+1)
		}
	}

}

// TetraDensity writes the density of a slice of tetrahedra to a slice of
// floats.
func tetrahedraDensity(hd *io.SheetHeader, ts []geom.Tetra, rhos []float64) {
	tw := hd.TotalWidth
	tCount := float64(hd.Count * 6)
	for i := range ts {
		rhos[i] = (tw * tw * tw) / (ts[i].Volume() * tCount)
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

// WrapXs updates a slice of vectors so that they will be as close to thee given
// sheet as periodic boundary conditions will allow.
func wrapXs(xs []rGeom.Vec, hd *io.SheetHeader) {
	tw := float32(hd.TotalWidth)
	for i := range xs {
		for j := 0; j < 3; j++ {
			if xs[i][j] < hd.Origin[j] { xs[i][j] += tw }
		}
	}
}
