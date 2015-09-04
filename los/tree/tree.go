package tree

import (
	"fmt"
	
	ct "github.com/phil-mansfield/consistent_trees"
)

const SCALE_FACTOR_MUL = 100000

func HaloHistories(files []string, ids []int) (hs [][]ct.Halo, steps [][]int) {
	hs, steps = make([][]ct.Halo, len(ids)), make([][]int, len(ids))

	for _, file := range files {
		ct.ReadTree(file)
		for i, id := range ids {
			if hs[i] != nil { continue }
			if idHs, idSteps, ok := findHalo(id); ok {
				hs[i], steps[i] = idHs, idSteps
			}
		}
		ct.DeleteTree()
	}

	for i, idSteps := range steps {
		if idSteps == nil {
			panic(fmt.Sprintf("Halo %d not dound in given files.", ids[i]))
		}
	}

	return hs, steps
}

func scaleStep(scale float64) int {
	return ct.GetHaloTree().ScaleFactorConv(int(scale * SCALE_FACTOR_MUL))
}

func findHalo(id int) ([]ct.Halo, []int, bool) {
	list, ok := ct.FindClosestScale(1)
	if !ok { panic("HaloTree doesn't contain scale-1 HaloList. Impossible.") }
	if h, ok := ct.LookupHaloInList(list, id); ok {
		return nil, nil, false
	} else {
		hs, steps := []ct.Halo{ h }, []int{ scaleStep(h.Scale()) }
		for {
			h, ok = h.Prog()
			if !ok { break }
			hs, steps = append(hs, h), append(steps, scaleStep(h.Scale()))
		}
		return hs, steps, true
	}
}
