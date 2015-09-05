package tree

import (
	"fmt"
	
	ct "github.com/phil-mansfield/consistent_trees"
)

const SCALE_FACTOR_MUL = 100000

func HaloHistories(files []string, roots []int) (ids [][]int, steps [][]int) {
	ids, steps = make([][]int, len(roots)), make([][]int, len(roots))

	for _, file := range files {
		fmt.Println(file)
		ct.ReadTree(file)
		for i, id := range roots {
			if ids[i] != nil { continue }
			if idHs, idSteps, ok := findHalo(id); ok {
				ids[i] = make([]int, len(idHs))
				for j, h := range idHs { ids[i][j] = h.ID() }
				steps[i] = idSteps
			}
		}
		ct.DeleteTree()
	}

	for i, idSteps := range steps {
		if idSteps == nil {
			panic(fmt.Sprintf("Halo %d not found in given files.", roots[i]))
		}
	}

	return ids, steps
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
