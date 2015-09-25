package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
	
	"github.com/phil-mansfield/gotetra/render/halo"

	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
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

func main() {
	valFlags, err := parseCmd()
	if err != nil { log.Fatal(err.Error()) }
	if len(valFlags) == 0 {
		err = printStdin()
		if err != nil { log.Fatal(err.Error()) }
	}
	
	ids, snaps, inVals, err := parseStdin()
	if err != nil { log.Fatal(err.Error()) }
	vals, err := readVals(ids, snaps, valFlags)
	if err != nil { log.Fatal(err.Error()) }

	printVals(ids, snaps, inVals, vals)
}


func parseCmd() ([]halo.Val, error) {
	flag.Parse()
	args := flag.Args()
	vals := make([]halo.Val, len(args))
	var ok bool
	for i, arg := range args {
		vals[i], ok = valMap[strings.ToLower(arg)]
		if !ok {
			return nil, fmt.Errorf("Flag %d, %s, not recognized.", i+1, arg)
		}
	}
	return vals, nil
}

func printStdin() error {
	lines, err := stdinLines()
	if err != nil { return err }
	fmt.Println(strings.Join(lines, "\n"))
	return nil
}

func stdinLines() ([]string, error) {
	bs, err := ioutil.ReadAll(os.Stdin)
	if err != nil {
		return nil, fmt.Errorf(
			"Error reading stdin: %s.", err.Error(),
		)
	}

	text := string(bs)
	return strings.Split(text, "\n"), nil
}

func parseStdin() (ids, snaps []int, inVals [][]float64, err error) {
	ids, snaps, inVals = []int{}, []int{}, [][]float64{}
	lines, err := stdinLines()
	if err != nil { return nil, nil, nil, err }
	for i, line := range lines {
		rawTokens := strings.Split(line, " ")
		tokens := make([]string, 0, len(rawTokens))
		for _, tok := range rawTokens {
			if len(tok) != 0 { tokens = append(tokens, tok) }
		}

		var (
			id, snap int
			vals []float64
			err error
		)
		switch {
		case len(tokens) == 0:
			continue
		case len(tokens) == 1:
			if tokens[0] == "" { continue }
			return nil, nil, nil, fmt.Errorf(
				"Line %d of stdin has 1 token, but 2 are required.", i + 1,
			)
		case len(tokens) >= 2:
			id, err = strconv.Atoi(tokens[0])
			if err != nil {
				return nil, nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[0],
				)
			} 
			snap, err = strconv.Atoi(tokens[1]) 
			if err != nil {
				return nil, nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[1],
				)
			}
			
			vals = make([]float64, len(tokens) - 2) 
			for i := range vals {
				vals[i], err = strconv.ParseFloat(tokens[i + 2], 64)
			}
		}

		ids = append(ids, id)
		snaps = append(snaps, snap)
		inVals = append(inVals, vals)
	}

	return ids, snaps, inVals, nil
}

func readVals(ids, snaps []int, valFlags []halo.Val) ([][]float64, error) {
	snapBins, idxBins := binBySnap(snaps, ids)
	vals := make([][]float64, len(ids))

	sortedSnaps := []int{}
	for snap := range snapBins {
		sortedSnaps = append(sortedSnaps, snap)
	}
	sort.Ints(sortedSnaps)

	for _, snap := range sortedSnaps {
		idSet := snapBins[snap]
		idxSet := idxBins[snap]

		var (
			snapVals [][]float64
			err error
		)
		if snap == -1 { // Handle blank halos.
			snapVals = make([][]float64, len(idSet))
			for i := range snapVals {
				snapVals[i] = make([]float64, len(valFlags))
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
	return vals, nil
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

	idWidth, snapWidth := 0, 0
	valWidths := make([]int, len(vals[0]))
	for i := range ids {
		iWidth := len(fmt.Sprintf("%d", ids[i]))
		sWidth := len(fmt.Sprintf("%d", snaps[i]))
		if iWidth > idWidth { idWidth = iWidth }
		if sWidth > snapWidth { snapWidth = sWidth }
		for j := range vals[i] {
			width := len(fmt.Sprintf("%.10g", vals[i][j]))
			if width > valWidths[j] { valWidths[j] = width }
		}
	}

	rowFmt := fmt.Sprintf("%%%dd %%%dd", idWidth, snapWidth)
	valFmts := make([]string, len(vals[0]))
	for i := range valFmts {
		rowFmt += fmt.Sprintf(" %%%d.10g", valWidths[i])
	}
	rowFmt += "\n"
	for i := range ids {
		args := append([]interface{}{ids[i], snaps[i]}, intr(vals[i])...)
		fmt.Printf(rowFmt, args...)
	}
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

