package main

import (
	"fmt"
	"flag"
	"io/ioutil"
	"log"
	"path"
	"strconv"
	"strings"

	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/los/tree"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
)

type IDType int
const (
	Rockstar IDType = iota
	M500c
	M200c
	M200m

    finderCells = 150
    overlapMult = 3
)

var (
	ms, rs, xs, ys, zs []float64
	rids []int
)

func main() {
	log.Println("gtet_id")
	
	// Parse command line
	var (
		idTypeStr string
		idStart, idEnd, snap, mult int
		allowSubhalos bool
	)

	flag.StringVar(&idTypeStr, "IDType", "Rockstar",
		"[Rockstar | M500c | M200c | M200m]")
	flag.IntVar(&mult, "Mult", 1, "Number of times to print each ID.")
	flag.IntVar(&snap, "Snap", -1, "The snapshot number.")
	flag.IntVar(&idStart, "IDStart", -1, "ID start range (inclusive).")
	flag.IntVar(&idEnd, "IDEnd", -1, "ID stop range (exclusive)")
	flag.BoolVar(&allowSubhalos, "AllowSubhalos", false,
		"Allow subhalo ids to be passed through.")
	flag.Parse()

	if idStart > idEnd {
		log.Fatalf("IDStart, %d, is larger than IDEnd, %d.", idStart, idEnd)
	} else if idStart != idEnd && idStart < 0 {
		log.Fatalf("Non-positive IDStart %d.")
	}

	if mult <= 0 { log.Fatal("Mult must be positive.") }

	idType, err := parseIDType(idTypeStr)
	if err != nil { err.Error() }
	if idType != Rockstar && snap == -1 {
		log.Fatalf("Must set the Snap flag if using a non-default IDType.")
	}

	snapNum, err := util.SnapNum()
	if err != nil {
		log.Fatalf(
			"Error encountered when finding rockstar directory: %s",err.Error(),
		)
	} else if (snap < 1 && idType != Rockstar) || snap > snapNum {
		log.Fatalf("Snap %d is out of bounds for %d snaps.", snap, snapNum)
	}

	// Get IDs and snapshots

	rawIds := getIDs(idStart, idEnd, flag.Args())
	if len(rawIds) == 0 { return }

	var ids, snaps []int
	switch idType {
	case Rockstar:
		if snap != -1 {
			snaps = make([]int, len(rawIds))
			for i := range snaps { snaps[i] = snap }
		} else {
			snaps, err = findSnaps(rawIds)
		}
		ids = rawIds
		if err != nil { log.Fatalf(err.Error()) }
	case M200m:
		snaps = make([]int, len(rawIds))
		for i := range snaps { snaps[i] = snap }
		ids, err = convertSortedIDs(rawIds, snap)
		if err != nil { log.Fatal(err.Error()) }
	default:
		log.Fatal("Unsupported IDType for now. Sorry :3")
	}

	// Tag subhalos, if neccessary.
	var isSub []bool
	if allowSubhalos {
		isSub = make([]bool, len(ids))
	} else {
		isSub, err = findSubs(ids, snaps)
		if err != nil { log.Fatal(err.Error()) }
	}

	// Output
	printIds(ids, snaps, isSub, mult)
}

func parseIDType(str string) (IDType, error) {
	switch strings.ToLower(str) {
	case "rockstar": return Rockstar, nil
	case "m500c": return M500c, nil
	case "m200c": return M200c, nil
	case "m200m": return M200m, nil
	}
	return -1, fmt.Errorf("IDType '%s' unrecognized", str)
}

func getIDs(idStart, idEnd int, args []string) []int {
	ids := make([]int, 0, idEnd - idStart + len(args))
	for i := idStart; i < idEnd; i++ { ids = append(ids, i) }
	for _, str := range args {
		i, err := strconv.Atoi(str)
		if err != nil {
			log.Fatalf("Could not parse arg %d: '%s' is not an int.", i+1, str)
		}
		ids = append(ids, i)
	}
	return ids
}

func getSnapHaloList(i int) (name string, err error) {	
	rockstarDir, err := util.RockstarDir()
	if err != nil { return "", err }
	infos, err := ioutil.ReadDir(rockstarDir)
	if err != nil { return "", err }
	return path.Join(rockstarDir, infos[i - 1].Name()), nil
}

func findSnaps(ids []int) ([]int, error) {
	treeDir, err := util.TreeDir()
	if err != nil { return nil, err }
	infos, err := ioutil.ReadDir(treeDir)
	if err != nil { return nil, err }

	names := []string{}
	for _, info := range infos {
		name := info.Name()
		n := len(name)
		if n > 4 && name[:5] == "tree_" && name[n-4:] == ".dat" {
			names = append(names, path.Join(treeDir, name))
		}
	}

	return tree.HaloSnaps(names, ids)
}

func readHeader(snap int) (*io.SheetHeader, error) {
	gtetFmt, err := util.GtetFmt()
	if err != nil { return nil, err }

	gtetDir := fmt.Sprintf(gtetFmt, snap)
	gtetFiles, err := util.DirContents(gtetDir)
	if err != nil { return nil, err }

	hd := &io.SheetHeader{}
	err = io.ReadSheetHeaderAt(gtetFiles[0], hd)
	if err != nil { return nil, err }
	return hd, nil
}

func convertSortedIDs(
	rawIDs []int, snap int,
) ([]int, error) {
	maxID := 0
	for _, id := range rawIDs {
		if id > maxID { maxID = id }
	}

	rids, err := util.ReadSortedRockstarIDs(snap, maxID, halo.M200b)
	if err != nil { return nil, err }

	ids := make([]int, len(rawIDs))
	for i := range ids { ids[i] = rids[rawIDs[i]] }
	return ids, nil
}

func findSubs(rawIDs, snaps []int) ([]bool, error) {
	isSub := make([]bool, len(rawIDs))

	// Group by snapshot.
	snapGroups := make(map[int][]int)
	groupIdxs := make(map[int][]int)
	for i, id := range rawIDs {
		snap := snaps[i]
		snapGroups[snap] = append(snapGroups[snap], id)
		groupIdxs[snap] = append(groupIdxs[snap], i)
	}

	// Load each snapshot.
	for snap, group := range snapGroups {
		hd, err := readHeader(snap)
		if err != nil { return nil, err }

		rids, err := util.ReadSortedRockstarIDs(snap, -1, halo.M200b)
		if err != nil { return nil, err }
		vals, err := util.ReadRockstar(
			snap, rids, halo.X, halo.Y, halo.Z, halo.Rad200b,
		)
		xs, ys, zs, rs := vals[0], vals[1], vals[2], vals[3]

		g := halo.NewGrid(finderCells, hd.TotalWidth, len(xs))
		g.Insert(xs, ys, zs)
		sf := halo.NewSubhaloFinder(g)
		sf.FindSubhalos(xs, ys, zs, rs, overlapMult)

		log.Println("gtet_id: search start")
		for i, id := range group {
			origIdx := groupIdxs[snap][i]
			for j, checkID := range rids {
				if checkID == id {
					isSub[origIdx] = sf.HostCount(j) > 0
					break
				} else if j == len(rids) - 1 {
					return nil, fmt.Errorf("ID %d not in halo list.", id)
				}
			}
		}
		log.Println("gtet_id: search start")
	}
	return isSub, nil
}

func printIds(ids []int, snaps []int, isSub []bool, mult int) {
	// Find the maximum width of each column.
	idWidth, snapWidth := 2, 2
	for i := range ids {
		if isSub[i] { continue }
		iWidth := len(fmt.Sprintf("%d", ids[i]))
		sWidth := len(fmt.Sprintf("%d", snaps[i]))
		if iWidth > idWidth { idWidth = iWidth }
		if sWidth > snapWidth { snapWidth = sWidth }
	}

	rowFmt := fmt.Sprintf("%%%dd %%%dd\n", idWidth, snapWidth)
	for i := range ids {
		if isSub[i] { continue }
		for j := 0; j < mult; j++ {
			fmt.Printf(rowFmt, ids[i], snaps[i])
		}
		if mult > 1 { fmt.Printf(rowFmt, -1, -1) }
	}
}
