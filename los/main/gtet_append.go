package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"path"
	"sort"
	"strings"
	
	"github.com/phil-mansfield/gotetra/render/halo"
	
	"github.com/phil-mansfield/gotetra/los/tree"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
)

const (
	GammaDZ = 0.5
)

var valMap = map[string]halo.Val {
	"scale":halo.Scale,
	"id":halo.ID,
    "descscale":halo.DescScale,
    "descid":halo.DescID,
    "numprog":halo.NumProg,
    "pid":halo.PID,
    "upid":halo.UPID,
    "descpid":halo.DescPID,
    "phantom":halo.Phantom,
    "sammvir":halo.SAMMVir,
    "mvir":halo.MVir,
    "rvir":halo.RVir,
    "rs":halo.Rs,
    "vrms":halo.Vrms,
    "mmp":halo.MMP,
    "scaleoflastmmp":halo.ScaleOfLastMMP,
    "vmax":halo.VMax,
    "x":halo.X,
    "y":halo.Y,
    "z":halo.Z,
    "vx":halo.Vx,
    "vy":halo.Vy,
    "vz":halo.Vz,
    "jx":halo.Jx,
    "jy":halo.Jy,
    "jz":halo.Jz,
    "spin":halo.Spin,
    "breadthfirstid":halo.BreadthFirstID,
    "depthfirstid":halo.DepthFirstID,
    "treerootid":halo.TreeRootID,
    "orighaloid":halo.OrigHaloID,
    "snapnum":halo.SnapNum,
    "nextcoprogenitordepthfirstid":halo.NextCoprogenitorDepthFirstID,
    "lastprogenitordepthfirstid":halo.LastProgenitorDepthFirstID,
    "rsklypin":halo.RsKylpin,
    "mvirall":halo.MVirAll,
    "m200b":halo.M200b,
    "m200c":halo.M200c,
    "m500c":halo.M500c,
    "m2500c":halo.M2500c,
    "xoff":halo.XOff,
    "voff":halo.Voff,
    "spinbullock":halo.SpinBullock,
    "btoa":halo.BToA,
    "ctoa":halo.CToA,
    "ax":halo.Ax,
    "ay":halo.Ay,
    "az":halo.Az,
    "btoa500c":halo.BToA500c,
    "ctoa500c":halo.CToA500c,
    "ax500c":halo.Ax500c,
    "ay500c":halo.Ay500c,
    "az500c":halo.Az500c,
    "tu":halo.TU,
    "macc":halo.MAcc,
    "mpeak":halo.MPeak,
    "vacc":halo.VAcc,
    "vpeak":halo.VPeak,
    "halfmassscale":halo.HalfmassScale,
    "accrateinst":halo.AccRateInst,
    "accrate100myr":halo.AccRate100Myr,
    "accratetdyn":halo.AccRateTdyn,
    "rad200b":halo.Rad200b,
    "rad200c":halo.Rad200c,
    "rad500c":halo.Rad500c,
    "rad2500c":halo.Rad2500c,
}

type GammaFlag int
const (
	DK14Gamma GammaFlag = iota
)

var gammaMap = map[string]GammaFlag{
	"dk14gamma": DK14Gamma,
}

func (gf GammaFlag) EvalAt(
	snapTrees [][]int, mvirs, scales [][]float64,
	evalSnaps []int, out []float64,
)  {
	for ti := range snapTrees {
		out[ti] = gf.Eval(
			snapTrees[ti], mvirs[ti], scales[ti], evalSnaps[ti],
		)
	}
}

func (gf GammaFlag) Eval(
	snaps []int, mvirs, scales []float64, evalSnap int,
) float64 {
	n := len(snaps)
	
	switch gf {
	case DK14Gamma:
		evalIdx := sort.SearchInts(snaps, evalSnap)
		if evalIdx < 0 || evalIdx >= n {
			panic(fmt.Sprintf("Snap %d not found", evalSnap))
		}
		evalZ := (1 / scales[evalIdx]) - 1
		backIdx := sort.SearchFloat64s(scales, 1 / (1 + evalZ + GammaDZ))
		return (math.Log(mvirs[evalIdx]) - math.Log(mvirs[backIdx])) /
			(math.Log(scales[evalIdx]) - math.Log(scales[backIdx]))
	}
	panic("Impossible :3")
}

type AppendFlag struct {
	ValFlag halo.Val
	GammaFlag GammaFlag
	IsGamma bool
}

func main() {	
	flags, err := parseCmd()
	if err != nil { log.Fatal(err.Error()) }
	
	ids, snaps, inVals, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }
	log.Println("gtet_append")
	vals, err := readVals(ids, snaps, flags)
	if err != nil { log.Fatal(err.Error()) }

	printVals(ids, snaps, inVals, vals)
}


func parseCmd() ([]AppendFlag, error) {
	flag.Parse()
	args := flag.Args()
	flags := make([]AppendFlag, len(args))
	for i, arg := range args {
		val, ok := valMap[strings.ToLower(arg)]
		if !ok {
			gamma, ok := gammaMap[strings.ToLower(arg)]
			if !ok {
				return nil, fmt.Errorf(
					"Flag %d, %s, not recognized.", i+1, arg,
				)
			} else {
				flags[i].GammaFlag = gamma
				flags[i].IsGamma = true
			}
		} else {
			flags[i].ValFlag = val
		}
	}
	return flags, nil
}

func readVals(ids, snaps []int, aFlags []AppendFlag) ([][]float64, error) {
	var (
		idTrees, snapTrees [][]int
		scales, mVirs [][]float64
	)
	if requiresTrees(aFlags) {
		trees, err := treeFiles()
		if err != nil { return nil, err }
		snapOffset, err := util.SnapOffset()
		if err != nil { return nil, err }
		idTrees, snapTrees, err = tree.HaloHistories(trees, ids, snapOffset)
		if err != nil { return nil, err }
		scales, mVirs, err = massHistories(idTrees, snapTrees)
		if err != nil { return nil, err }
	}
	
	snapBins, idxBins := binBySnap(snaps, ids)
	vals := make([][]float64, len(ids))

	sortedSnaps := []int{}
	for snap := range snapBins {
		sortedSnaps = append(sortedSnaps, snap)
	}
	sort.Ints(sortedSnaps)

	// Now comes the part where we gaze whistfully into the distance while
	// wishing Go had algebraic data types.
	
	valFlags := make([]halo.Val, len(aFlags))
	for i := range valFlags {
		if aFlags[i].IsGamma {
			valFlags[i] = halo.Scale
		} else {
			valFlags[i] = aFlags[i].ValFlag
		}
	}
	
	for _, snap := range sortedSnaps {
		if snap == 82 || snap == 84 { continue }
		idSet := snapBins[snap]
		idxSet := idxBins[snap]

		var (
			snapVals [][]float64
			err error
		)

		if snap == -1 { // Handle blank halos.
			snapVals = make([][]float64, len(idSet))
			for i := range snapVals {
				snapVals[i] = make([]float64, len(aFlags))
			}
		} else {
			snapVals, err = util.ReadRockstar(snap, idSet, valFlags...)			
			if err != nil { return nil, err }
			snapVals = flipAxis(snapVals)
		}

		for i := range idSet {
			vals[idxSet[i]] = snapVals[i]
		}
	}

	for fi := range aFlags {
		if !aFlags[fi].IsGamma { continue }
		for ti := range snapTrees {
			vals[ti][fi] = aFlags[fi].GammaFlag.Eval(
				snapTrees[ti], mVirs[ti], scales[ti], snaps[ti],
			)
		}

	}
	
	return vals, nil
}

func requiresTrees(aFlags []AppendFlag) bool {
	for _, flag := range aFlags {
		if flag.IsGamma { return true }
	}
	return false
}

// I hate this function.
func massHistories(
	idTrees, snapTrees [][]int,
) (scales, mvirs [][]float64, err error) {
	
	minSnap, maxSnap := snapTrees[0][0], snapTrees[0][0]
	for _, tree := range snapTrees {
		for _, snap := range tree {
			if snap > maxSnap {
				maxSnap = snap
			} else if snap < minSnap {
				minSnap = snap
			}
		}
	}

	// So much sadness... The point of this is so that I only have to read each
	// catalog file once.
	rids := make([]int, len(idTrees))
	mvirs = make([][]float64, len(idTrees))
	scales = make([][]float64, len(idTrees))
	for i := range mvirs {
		mvirs[i] = make([]float64, len(idTrees[i]))
		scales[i] = make([]float64, len(idTrees[i]))
	}
	
	for snap := minSnap; snap <= maxSnap; snap++ {
		if snap == 82 || snap == 84 { continue }
		defaultID := -1
		for i := range rids { rids[i] = -1 }
		for ti, tree := range snapTrees {
			for si, tSnap := range tree {
				if snap == tSnap {
					defaultID = idTrees[ti][si]
					rids[ti] = defaultID
					break
				}
			}
		}
		
		// Don't read catalogs that no tree actually requires.
		if defaultID == -1 { continue }

		// Trees that don't exist at this snapshot will read some dummy value.
		// This is just to reduce the ugliness of this already ugly code.
		for i := range rids {
			if rids[i] == -1 { rids[i] = defaultID }
		}

		vals, err := util.ReadRockstar(snap, rids, halo.Scale, halo.MVir)
		if err != nil { return nil, nil, err }
		sScales, sMVirs := vals[0], vals[1]

		for ti, tree := range snapTrees {
			for si, tSnap := range tree {
				if snap == tSnap {
					mvirs[ti][si] = sMVirs[ti]
					scales[ti][si] = sScales[ti]
					break
				}
			}
		}
	}

	return scales, mvirs, nil
}

func binBySnap(snaps, ids []int) (snapBins, idxBins map[int][]int) {
	snapBins = make(map[int][]int)
	idxBins = make(map[int][]int)
	for i, snap := range snaps {
		id := ids[i]
		snapBins[snap] = append(snapBins[snap], id)
		idxBins[snap] = append(idxBins[snap], i)
	}
	return snapBins, idxBins
}

func printVals(ids, snaps []int, inVals, vals [][]float64) {
	for i := range inVals {
		inVals[i] = append(inVals[i], vals[i]...)
	}
	vals = inVals

	util.PrintRows(ids, snaps, vals)
}

func flipAxis(vals [][]float64) [][]float64 {
    out := make([][]float64, len(vals[0]))
    for i := range out { out[i] = make([]float64, len(vals)) }
    for i := range out {
        for j := range vals {
            out[i][j] = vals[j][i]
        }
    }
    return out
}

func intr(xs []float64) []interface{} {
	is := make([]interface{}, len(xs))
	for i, x := range xs { is[i] = x }
	return is
}

func treeFiles() ([]string, error) {
	treeDir, err := util.TreeDir()
	if err != nil { return nil, err }
	infos, err := ioutil.ReadDir(treeDir)
	if err != nil { return nil, err }

	names := []string{}
	for _, info := range infos {
		name := info.Name()
		n := len(name)
		// This is pretty hacky.
		if n > 4 && name[:5] == "tree_" && name[n-4:] == ".dat" {
			names = append(names, path.Join(treeDir, name))
		}
	}
	return names, nil
}
