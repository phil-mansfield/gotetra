package tree

import (
	"fmt"
	
	ct "github.com/phil-mansfield/consistent_trees"
)

const SCALE_FACTOR_MUL = 100000

func HaloHistories(
	files []string, roots []int,
) (ids [][]int, snaps [][]int, err error) {

	ids, snaps = make([][]int, len(roots)), make([][]int, len(roots))

	for _, file := range files {
		ct.ReadTree(file)
		for i, id := range roots {
			if ids[i] != nil { continue }
			if idHs, idSnaps, ok := findHalo(id); ok {
				ids[i] = make([]int, len(idHs))
				for j, h := range idHs { ids[i][j] = h.ID() }
				snaps[i] = idSnaps
			}
		}
		ct.DeleteTree()
	}

	for i, idSnaps := range snaps {
		if idSnaps == nil {
			return nil, nil, fmt.Errorf(
				"Halo %d not found in given files.", roots[i],
			)
		}
	}

	return ids, snaps, nil
}

func HaloSnaps(files []string, ids []int) (snaps []int, err error) {
	snaps = make([]int, len(ids))
	for i := range snaps { snaps[i] = -1 }

	numLists := ct.GetHaloTree().NumLists()

	for _, file := range files {
		ct.ReadTree(file)
		for i, id := range ids {
			if snaps[i] != -1 { continue }
			if h, ok := ct.LookupHaloInList(ct.GetAllHalos(), id); ok {
				snaps[i] = numLists - ct.LookupIndex(h.Scale()) + 1
			}
		}
		ct.DeleteTree()
	}

	for i, snap := range snaps {
		if snap == -1 {
			return nil, fmt.Errorf(
				"Halo %d not found in given files.", ids[i],
			)
		}
	}

	return snaps, nil
}

func findHalo(id int) ([]ct.Halo, []int, bool) {
	list, ok := ct.FindClosestScale(1)
	if !ok { panic("HaloTree doesn't contain scale-1 HaloList. Impossible.") }
	if h, ok := ct.LookupHaloInList(list, id); !ok {
		return nil, nil, false
	} else {
		hs, steps := []ct.Halo{ h }, []int{ ct.LookupIndex(h.Scale()) }
		for {
			h, ok = h.Prog()
			if !ok { break }
			hs, steps = append(hs, h), append(steps, ct.LookupIndex(h.Scale()))
		}
		return hs, steps, true
	}
}
