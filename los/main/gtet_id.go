package main

import (
	"fmt"
	"flag"
	"io/ioutil"
	"log"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/los/tree"
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
		log.Fatalf("IDStart, %d, is larger than IDEnd, %d.")
	} else if idStart != idEnd && idStart <= 0 {
		log.Fatalf("Non-positive IDStart %d.")
	}

	if mult <= 0 { log.Fatal("Mult must be positive.") }

	idType, err := parseIDType(idTypeStr)
	if err != nil { err.Error() }
	if idType != Rockstar && snap == -1 {
		log.Fatalf("Must set the Snap flag if using a non-default IDType.")
	}

	snapNum, err := snapNum()
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
		if err != nil { err.Error() }
	default:
		log.Fatal("Unsupported IDType for now. Sorry :3")
	}

	// Tag subhalos, if neccessary.
	var isSub []bool
	if allowSubhalos {
		isSub = make([]bool, len(ids))
	} else {
		isSub, err = findSubs(rawIds, snaps)
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

func snapNum() (int, error) {
	rockstarDir := os.Getenv("GTET_ROCKSTAR_DIR")
	infos, err := ioutil.ReadDir(rockstarDir)
	return len(infos), err
}

func getSnapHaloList(i int) (name string, err error) {
	rockstarDir := os.Getenv("GTET_ROCKSTAR_DIR")
	infos, err := ioutil.ReadDir(rockstarDir)
	if err != nil { return "", err }
	return path.Join(rockstarDir, infos[i - 1].Name()), nil
}

func findSnaps(ids []int) ([]int, error) {
	treeDir := os.Getenv("GTET_TREE_DIR")
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


// This funciton has side effects :(
func convertSortedIDs(
	rawIDs []int, snap int,
) ([]int, error) {
	list, err := getSnapHaloList(snap)
	if err != nil { return nil, err }

	gtetFmt := os.Getenv("GTET_FMT")
	gtetName := fmt.Sprintf(gtetFmt, snap, 0, 0, 0)
	hd := &io.SheetHeader{}
	err = io.ReadSheetHeaderAt(gtetName, hd)
	if err != nil { return nil, err }
	cosmo := &hd.Cosmo

	rids, xs, ys, zs, ms, rs, err = halo.ReadRockstar(list, halo.R200m, cosmo)
	if err != nil { return nil, err }
	
	ids := make([]int, len(rawIDs))
	for i := range ids { ids[i] = rids[rawIDs[i]] }
	return ids, nil
}

func findSubs(rawIDs, snaps []int) ([]bool, error) {
	isSub := make([]bool, len(rawIDs))

	// Handle the case where we're indexing by sorted mass ID.
	if xs != nil {
		gtetFmt := os.Getenv("GTET_FMT")
		gtetName := fmt.Sprintf(gtetFmt, snaps[0], 0, 0, 0)
		hd := &io.SheetHeader{}
		err := io.ReadSheetHeaderAt(gtetName, hd)
		if err != nil { return nil, err }

		g := halo.NewGrid(finderCells, hd.TotalWidth, len(xs))
		g.Insert(xs, ys, zs)
		sf := halo.NewSubhaloFinder(g)
		sf.FindSubhalos(xs, ys, zs, rs, overlapMult)

		for i, id := range rawIDs {
			isSub[i] = sf.HostCount(id) > 0
		}
		return isSub, nil
	}

	// Handle the case where all we have are halo IDs
	
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
		gtetFmt := os.Getenv("GTET_FMT")
		gtetName := fmt.Sprintf(gtetFmt, snap, 0, 0, 0)
		hd := &io.SheetHeader{}
		err := io.ReadSheetHeaderAt(gtetName, hd)
		if err != nil { return nil, err }

		list, err := getSnapHaloList(snap)
		rids, xs, ys, zs, ms, rs, err = halo.ReadRockstar(
			list, halo.R200m, &hd.Cosmo,
		)

		g := halo.NewGrid(finderCells, hd.TotalWidth, len(xs))
		g.Insert(xs, ys, zs)
		sf := halo.NewSubhaloFinder(g)
		sf.FindSubhalos(xs, ys, zs, rs, overlapMult)
		
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
		for i := 0; i < mult; i++ {
			fmt.Printf(rowFmt, ids[i], snaps[i])
		}
		if mult > 1 { fmt.Printf(rowFmt, -1, -1) }
	}
}
