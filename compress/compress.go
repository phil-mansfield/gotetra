package compress

import (
	"github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/io"
)

func Normalize(vs []geom.Vec, hd *io.SheetHeader) {
	n := len(vs) - 1
	tw := float32(hd.TotalWidth)
	tw2 := float32(hd.TotalWidth / 2)
	for i := range vs[1:n-1] {
		pv := vs[i]
		v := vs[i+1]
		for dim := 0; dim < 3; dim++ {
			dx := v[dim] - pv[dim]
			if dx > tw2 {
				vs[i+1][dim] -= tw
			} else if dx < -tw2 {
				vs[i+1][dim] += tw
			}
		}
	}
}
