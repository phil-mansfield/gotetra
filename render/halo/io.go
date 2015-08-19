package halo

import (
	"sort"
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/table"
)

// halos allows for arrays of halo properties to be sorted simultaneously.
type halos struct {
	xs, ys, zs, ms, rs []float64
}

func (hs *halos) Len() int { return len(hs.rs) }
func (hs *halos) Less(i, j int) bool { return hs.rs[i] < hs.rs[j] }
func (hs *halos) Swap(i, j int) {
	hs.rs[i], hs.rs[j] = hs.rs[j], hs.rs[i]
	hs.ms[i], hs.ms[j] = hs.ms[j], hs.ms[i]
	hs.xs[i], hs.xs[j] = hs.xs[j], hs.xs[i]
	hs.ys[i], hs.ys[j] = hs.ys[j], hs.ys[i]
	hs.zs[i], hs.zs[j] = hs.zs[j], hs.zs[i]
}

// ReadRockstar reads halo information from the given Rockstar catalog, sorted
// from largest to smallest.
func ReadRockstar(
	file string, rType Radius, cosmo *io.CosmologyHeader,
) (xs, ys, zs, ms, rs []float64, err error) {
	rCol := rType.RockstarColumn()
	xCol, yCol, zCol := 17, 18, 19
	
	colIdxs := []int{ xCol, yCol, zCol, rCol }
	cols, err := table.ReadTable(file, colIdxs, nil)
	if err != nil { return nil, nil, nil, nil, nil, err }
	
	xs, ys, zs = cols[0], cols[1], cols[2]
	if rType.RockstarMass() {
		ms = cols[3]
		rs = make([]float64, len(ms))
		rType.Radius(cosmo, ms, rs)
	} else {
		rs = cols[3]
		ms = make([]float64, len(rs))
		for i := range rs { rs[i] /= 1000 } // kpc -> Mpc
		rType.Mass(cosmo, rs, ms)
	}

	sort.Sort(sort.Reverse(&halos{ xs, ys, zs, ms, rs }))
	return xs, ys, zs, ms, rs, nil
}
