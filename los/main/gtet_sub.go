package main

import (
	"fmt"
	"log"
	"math"

	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/los/geom"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
)

const (
	finderCells = 150
	overlapMult = 3
)

type Params struct {
}

func main() {
	p := parseCmd()
	ids, snaps, _, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }

	disps, ms, err := subDisplacements(ids, snaps, p)
	if err != nil { log.Fatal(err.Error()) }
	rs := dispToRad(disps)

	for i := range rs { rs[i] = append(rs[i], ms[i]...) }
	util.PrintRows(ids, snaps, rs)
}

func parseCmd() *Params { return &Params{} }

func dispToRad(disps [][]geom.Vec) [][]float64 {
	rs := make([][]float64, len(disps))
	for i, ds := range disps {
		rs[i] = make([]float64, len(ds))

		for j := range ds {
			sum := float32(0)
			for _, comp := range ds[j] { sum += comp*comp }
			rs[i][j] = math.Sqrt(float64(sum))
		}
	}

	return rs
}

func subDisplacements(ids, snaps []int, p *Params) (
	subDisps [][]geom.Vec, subMs [][]float64, err error,
) {
	snapGroups := make(map[int][]int)
	groupIdxs := make(map[int][]int)
	for i, id := range ids {
		snap := snaps[i]
		snapGroups[snap] = append(snapGroups[snap], id)
		groupIdxs[snap] = append(groupIdxs[snap], i)
	}

	subDisps = make([][]geom.Vec, len(ids))
	subMs = make([][]float64, len(ids))

	for snap, group := range snapGroups {
		// Read position data
		hd, err := util.ReadSnapHeader(snap)
		if err != nil { return nil, nil, err }
		sids, err := util.ReadSortedRockstarIDs(snap, -1, halo.M200b)
		if err != nil { return nil, nil, err }
		vals, err := util.ReadRockstar(
			snap, sids, halo.X, halo.Y, halo.Z, halo.Rad200b, halo.M200c,
		)
		xs, ys, zs, rs, ms := vals[0], vals[1], vals[2], vals[3], vals[4]

		// Find subhalos
		g := halo.NewGrid(finderCells, hd.TotalWidth, len(xs))
		g.Insert(xs, ys, zs)
		sf := halo.NewSubhaloFinder(g)
		sf.FindSubhalos(xs, ys, zs, rs, overlapMult)

		f := util.NewIntFinder(sids)
		// Convert subhalo ids to dispacement vectors
		for i, id := range group {
			origIdx := groupIdxs[snap][i]
			sIdx, ok := f.Find(id)
			if !ok { return nil, nil, fmt.Errorf("Could not find ID %d", id) }

			subs := sf.Subhalos(sIdx)
			subDisps[origIdx] = make([]geom.Vec, len(subs))
			subMs[origIdx] = make([]float64, len(subs))

			x, y, z := xs[sIdx], ys[sIdx], zs[sIdx]

			for j := range subs {
				sx, sy, sz := xs[subs[j]], ys[subs[j]], zs[subs[j]]
				subDisps[origIdx][j] = geom.Vec{
					float32(sx -x), float32(sy - y), float32(sz - z),
				}
				subMs[origIdx][j] = ms[subs[j]]
			}
		}
	}

	return subDisps, subMs, nil
}
