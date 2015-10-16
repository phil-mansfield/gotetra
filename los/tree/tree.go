package tree

import (
	"fmt"
	
	ct "github.com/phil-mansfield/consistent_trees"
)

const SCALE_FACTOR_MUL = 100000

func HaloHistories(
	files []string, roots []int, snapOffset int,
) (ids [][]int, snaps [][]int, err error) {
	if len(roots) == 0 { return [][]int{}, [][]int{}, nil }

	ids, snaps = make([][]int, len(roots)), make([][]int, len(roots))

	foundCount := 0
	for _, file := range files {
		ct.ReadTree(file)
		var ok bool
		for i, id := range roots {
			if ids[i] != nil { continue }
			if ids[i], snaps[i], ok = findHistory(id); ok {
				foundCount++
			}
		}
		ct.DeleteTree()
		if foundCount == len(roots) { break }
	}

	for i, idSnaps := range snaps {
		if idSnaps == nil {
			return nil, nil, fmt.Errorf(
				"Halo %d not found in given files.", roots[i],
			)
		}
	}

	for i := range snaps {
		for j := range snaps[i] { snaps[i][j] += snapOffset }
	}
	return ids, snaps, nil
}

func HaloSnaps(files []string, ids []int) (snaps []int, err error) {
	if len(ids) == 0 { return []int{}, nil }

	snaps = make([]int, len(ids))
	for i := range snaps { snaps[i] = -1 }

	foundCount := 0
	for _, file := range files {
		ct.ReadTree(file)
		for i, id := range ids {
			if snaps[i] != -1 { continue }
			if _, snap, ok := findHalo(id); ok {
				snaps[i] = snap
				foundCount++
			}
		}
		ct.DeleteTree()
		if foundCount == len(ids) { break }
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

func findHalo(id int) (ct.Halo, int, bool) {
	tree := ct.GetHaloTree()
	for i := 0; i < tree.NumLists(); i++ {
		list := tree.HaloLists(i)
		h, ok := ct.LookupHaloInList(list, id)
		if ok {
			return h, tree.NumLists() - i, true
		}
	}
	return ct.Halo{}, 0, false
}

func findHistory(id int) (ids, snaps []int, ok bool) {
	h, snap, ok := findHalo(id)
	if !ok { return nil, nil, false }
	desc, descSnaps := descTree(h)
	prog, progSnaps := progTree(h)

	ids = combine(reverse(prog), []int{id}, desc)
	snaps = combine(reverse(progSnaps), []int{snap}, descSnaps)
	return ids, snaps, true
}

func descTree(h ct.Halo) (ids, snaps []int) {
	ids, snaps = []int{}, []int{}
	var ok bool
	numLists := ct.GetHaloTree().NumLists()
	for {
		h, ok = h.Desc()
		if !ok { break }
		ids = append(ids, h.ID())
		snaps = append(snaps, numLists - ct.LookupIndex(h.Scale()))
	}
	return ids, snaps
}

func progTree(h ct.Halo) (ids, snaps []int) {
	ids, snaps = []int{}, []int{}
	var ok bool
	numLists := ct.GetHaloTree().NumLists()
	for {
		h, ok = h.Prog()
		if !ok { break }
		ids = append(ids, h.ID())
		snaps = append(snaps, numLists - ct.LookupIndex(h.Scale()))
	}
	return ids, snaps
}

func reverse(xs []int) []int {
	out := make([]int, len(xs))
	for i := range xs {
		out[i] = xs[len(xs) - 1 - i]
	}
	return out
}

func combine(slices ...[]int) []int {
	out := []int{}
	for _, slice := range slices {
		out = append(out, slice...)
	}
	return out
}
