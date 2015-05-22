package density

import (
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/phil-mansfield/gotetra/rand"	
)

type Interpolator interface {
	Interpolate(
		buf Buffer, xs, vs []geom.Vec,
		ptVal float64, weights Buffer,
		low, high, jump int,
	)

	// The bounding box around the box being written to.
	DomainCellBounds() *geom.CellBounds
	// The bounding box around the box being written from.
	BufferCellBounds() *geom.CellBounds

	Cells() int
}

type ngp struct { }

// The ordering of these fields makes no goddamned sense.
type mcarlo struct {
	subIntr Interpolator
	segWidth int64
	points, cells int

	gen *rand.Generator

	skip int64
	// Buffers
	idxBuf geom.TetraIdxs
	tet geom.Tetra
	vtet geom.Tetra

	unitBufs [][]geom.Vec
	vecBuf []geom.Vec
}

func MonteCarlo(
	segWidth int64,
	points, cells int,
	skip int64,
	unitBufs [][]geom.Vec,
	subIntr Interpolator,
) Interpolator {
	mc := &mcarlo{
		subIntr, segWidth, points, cells,
		rand.NewTimeSeed(rand.Xorshift), skip,
		geom.TetraIdxs{}, geom.Tetra{}, geom.Tetra{},
		unitBufs, make([]geom.Vec, points),
	}

	return mc
}

func (intr *mcarlo) DomainCellBounds() *geom.CellBounds {
	return intr.subIntr.DomainCellBounds()
}

func (intr *mcarlo) BufferCellBounds() *geom.CellBounds {
	return intr.subIntr.BufferCellBounds()
}

func (intr *mcarlo) Cells() int {
	return intr.subIntr.Cells()
}

func (intr *mcarlo) Interpolate(
	buf Buffer, xs, vs []geom.Vec,
	ptVal float64, weights Buffer,
	low, high, jump int,
) {
	// weights will contain the velocity if needed and will be a NilBuffer
	// otherwise.

	segWidth := intr.segWidth
	gridWidth := segWidth + 1
	idxWidth := intr.segWidth / intr.skip

	ptVal = ptVal / float64(intr.points) / 6.0 *
		float64(intr.skip * intr.skip * intr.skip)

	tetCb := &geom.CellBounds{}

	relCb := &geom.CellBounds{}
	relCb.Width = intr.DomainCellBounds().Width
	relCb.Origin[0], relCb.Origin[1], relCb.Origin[2] =
		cbSubtr(intr.DomainCellBounds(), intr.BufferCellBounds())

	jump64 := int64(jump)

	// I hate this so much. I hate it so much:
	//
	// (Also, I have no idea how it works.)
	neverMod := *intr.subIntr.DomainCellBounds() !=
		*intr.subIntr.BufferCellBounds() ||
		(intr.Cells() / 2 > relCb.Width[0] &&
		intr.Cells() / 2 > relCb.Width[1] &&
		intr.Cells() / 2 > relCb.Width[2])
	mods := 0

	reqVel := buf.Quantity().RequiresVelocity()

	for idx := int64(low); idx < int64(high); idx += jump64 {
		x, y, z := coords(idx, idxWidth)
		gridIdx := index(x * intr.skip, y * intr.skip, z * intr.skip, gridWidth)

		for dir := 0; dir < 6; dir++ {
			intr.idxBuf.Init(gridIdx, gridWidth, intr.skip, dir)
			
			intr.tet.Init(
				&xs[intr.idxBuf[0]],
				&xs[intr.idxBuf[1]],
				&xs[intr.idxBuf[2]],
				&xs[intr.idxBuf[3]],
			)

			if reqVel {
				intr.vtet.Init(
					&vs[intr.idxBuf[0]],
					&vs[intr.idxBuf[1]],
					&vs[intr.idxBuf[2]],
					&vs[intr.idxBuf[3]],
				)
			}

			intr.tet.CellBoundsAt(1.0, tetCb)

			if !tetCb.Intersect(relCb, intr.cells) { continue }

			bufIdx := intr.gen.UniformInt(0, len(intr.unitBufs))
			intr.tet.DistributeTetra(
				intr.unitBufs[bufIdx],
				intr.vecBuf,
			)

			// Lol, whatever.
			if !neverMod {
				for j := 0; j < 3; j++ {
					if coordNeg(&intr.tet, j) {
						mods++
						modCoord(intr.vecBuf, j, float32(intr.Cells()))
					}
				}
			}

			var bweights Buffer
			switch b := buf.(type) {
			case *densityBuffer, *gradientBuffer:
				bweights = NilBuffer
			case *velocityBuffer:
				wbuf, ok := b.weights.VectorBuffer()
				if !ok { panic("buf is non-vector when vector is required.") }
				intr.vtet.DistributeTetra64(intr.unitBufs[bufIdx], wbuf)
				bweights = b.weights
			case *curlBuffer:
				wbuf, ok := b.weights.VectorBuffer()
				if !ok { panic("buf is non-vector when vector is required.") }
				intr.vtet.DistributeTetra64(intr.unitBufs[bufIdx], wbuf)
				bweights = b.weights
			case *divergenceBuffer:
				wbuf, ok := b.weights.VectorBuffer()
				if !ok { panic("buf is non-vector when vector is required.") }
				intr.vtet.DistributeTetra64(intr.unitBufs[bufIdx], wbuf)
				bweights = b.weights
			}

			intr.subIntr.Interpolate(
				buf, intr.vecBuf, nil, ptVal, bweights, 0, intr.points, 1,
			)
		}
	}
}

func coordNeg(tet *geom.Tetra, j int) bool {
	for i := 0; i < 4; i++ {
		if tet.Corners[i][j] < 0 { return true }
	}
	return false
}

func modCoord(buf []geom.Vec, j int, width float32) {
	for i := range buf {
		if buf[i][j] < 0 { buf[i][j] += width }
	}
}

func index(x, y, z, cells int64) int64 {
	return x + y * cells + z * cells * cells
}

func coords(idx, cells int64) (x, y, z int64) {
	x = idx % cells
	y = (idx % (cells * cells)) / cells
	z = idx / (cells * cells)
	return x, y, z
}

func cbSubtr(cb1, cb2 *geom.CellBounds) (i, j, k int) {
	i = cb1.Origin[0] - cb2.Origin[0]
	j = cb1.Origin[1] - cb2.Origin[1]
	k = cb1.Origin[2] - cb2.Origin[2]
	return i, j, k
}
