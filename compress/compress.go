package compress

import (
	"github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/io"
)

func Normalize(vs []geom.Vec, axis int, hd *io.SheetHeader) {
	for dim := 0; dim < 3; dim++ {
		if dim == axis {
			dx := float32(hd.TotalWidth / float64(hd.CountWidth))
			for i := range vs {
				x0 := float32(i) * dx
				vs[i][dim] -= x0
			}

		} else {
			sum := float32(0)
			for _, v := range vs {
				sum += v[dim]
			}
			mean := sum / float32(len(vs))
			for i := range vs {
				vs[i][dim] -= mean
			}
		}
	}
}
