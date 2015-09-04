package halo

import (
	"sort"

	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/table"
)

// halos allows for arrays of halo properties to be sorted simultaneously.
type halos struct {
	rids []int
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
	hs.rids[i], hs.rids[j] = hs.rids[j], hs.rids[i]
}

// ReadRockstar reads halo information from the given Rockstar catalog, sorted
// from largest to smallest.
func ReadRockstar(
	file string, rType Radius, cosmo *io.CosmologyHeader,
) (rids []int, xs, ys, zs, ms, rs []float64, err error) {
	rCol := rType.RockstarColumn()
	idCol, xCol, yCol, zCol := 1, 17, 18, 19
	
	colIdxs := []int{ idCol, xCol, yCol, zCol, rCol }
	cols, err := table.ReadTable(file, colIdxs, nil)
	if err != nil { return nil, nil, nil, nil, nil, nil, err }
	
	ids := cols[0]
	xs, ys, zs = cols[1], cols[2], cols[3]
	if rType.RockstarMass() {
		ms = cols[4]
		rs = make([]float64, len(ms))
		rType.Radius(cosmo, ms, rs)
	} else {
		rs = cols[4]
		ms = make([]float64, len(rs))
		for i := range rs { rs[i] /= 1000 } // kpc -> Mpc
		rType.Mass(cosmo, rs, ms)
	}

	rids = make([]int, len(ids))
	for i := range rids { rids[i] = int(ids[i]) }

	sort.Sort(sort.Reverse(&halos{ rids, xs, ys, zs, ms, rs }))
	return rids, xs, ys, zs, ms, rs, nil
}
