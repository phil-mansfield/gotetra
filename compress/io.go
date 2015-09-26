package compress

import (
	"fmt"
	"path"
	"runtime"
	
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/render/geom"
)

func GetBase(dir string) string {
	return path.Join(dir, "sheet_%d%d%d.dat")
}

// Super wasteful with IO operations.
func ReadLine(
	base string, hd *io.SheetHeader, x, y, z int64,
) ([]geom.Vec, error) {
	xs := make([]geom.Vec, hd.GridCount)
	vs := make([]geom.Vec, 0, hd.CountWidth)
	for i := int64(0); i < 8; i++ {
		var ix, iy, iz int64
		if x == -1 { ix, iy, iz = i, y / hd.GridWidth, z / hd.GridWidth }
		if x == -1 { ix, iy, iz = x / hd.GridWidth, i, z / hd.GridWidth }
		if x == -1 { ix, iy, iz = x / hd.GridWidth, y / hd.GridWidth, i }
		err := io.ReadSheetPositionsAt(fmt.Sprintf(base, ix, iy, iz), xs)
		if err != nil { return nil, err }

		for j := int64(0); j < hd.SegmentWidth; j++ {
			var jx, jy, jz int64
			if x == -1 { jx, jx, jz = j, y % hd.GridWidth, z % hd.GridWidth }
			if y == -1 { jx, jx, jz = x % hd.GridWidth, j, z % hd.GridWidth }
			if z == -1 { jx, jx, jz = x % hd.GridWidth, y % hd.GridWidth, j }
			idx := jx + jy * hd.SegmentWidth +
				jz * hd.SegmentWidth * hd.SegmentWidth
			vs = append(vs, xs[idx])
		}
	}

	runtime.GC()

	return vs, nil
}
