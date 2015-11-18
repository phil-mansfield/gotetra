package gtet_util

import (
	"fmt"
	
	"github.com/phil-mansfield/gotetra/render/halo"
)

var (
	FinderCells = 150
	OverlapMult = 3.0
)

func SubIDs(ids, snaps []int) (
	sIDs, sSnaps, hIDs []int, err error,
) {
	snapGroups := make(map[int][]int)
	groupIdxs := make(map[int][]int)
	for i, id := range ids {
		snap := snaps[i]
		snapGroups[snap] = append(snapGroups[snap], id)
		groupIdxs[snap] = append(groupIdxs[snap], i)
	}

	subIDs := make([][]int, len(ids))
	
	for snap, group := range snapGroups {
		// Read position data
		hd, err := ReadSnapHeader(snap)
		if err != nil { return nil, nil, nil, err }
		allIDs, err := ReadSortedRockstarIDs(snap, -1, halo.M200b)
		if err != nil { return nil, nil, nil, err }
		vals, err := ReadRockstar(
			snap, allIDs, halo.X, halo.Y, halo.Z, halo.Rad200b,
		)
		xs, ys, zs, rs := vals[0], vals[1], vals[2], vals[3]
		
		// Find subhalos
		g := halo.NewGrid(FinderCells, hd.TotalWidth, len(xs))
		g.Insert(xs, ys, zs)
		sf := halo.NewSubhaloFinder(g)
		sf.FindSubhalos(xs, ys, zs, rs, OverlapMult)

		f := NewIntFinder(allIDs)
		// Convert subhalo ids to dispacement vectors
		for i, id := range group {
			origIdx := groupIdxs[snap][i]
			sIdx, ok := f.Find(id)
			if !ok {
				return nil, nil, nil, fmt.Errorf("Could not find ID %d", id)
			}
			subMassIDs := sf.Subhalos(sIdx)
			subIDs[origIdx] = make([]int, len(subMassIDs))
			for i := range subMassIDs {
				subIDs[origIdx][i] = allIDs[subMassIDs[i]]
			}
		}
	}

	sIDs, sSnaps, hIDs = []int{}, []int{}, []int{}
	for i := range subIDs {
		sIDs = append(sIDs, ids[i])
		sSnaps = append(sSnaps, snaps[i])
		hIDs = append(hIDs, ids[i])
		for _, id := range subIDs[i] {
			sIDs = append(sIDs, id)
			sSnaps = append(sSnaps, snaps[i])
			hIDs = append(hIDs, ids[i])
		}
	}
	
	return sIDs, sSnaps, hIDs, nil
}

