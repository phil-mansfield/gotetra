package gtet_util

import (
	"encoding/binary"
	"fmt"
	"os"
	"path"
	"sort"
	
	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/render/io"
)

const (
	rockstarMemoDir = "rockstar"
	rockstarMemoFile = "halo_%d.dat"
	rockstarShortMemoFile = "halo_short_%d.dat"

	RockstarShortMemoNum = 10 * 1000
)

func ReadSortedRockstarIDs(snap, maxID int, flag halo.Val) ([]int, error) {
	memoDir, err := MemoDir()
	if err != nil { return nil, err }
	dir := path.Join(memoDir, rockstarMemoDir)
	if !PathExists(dir) { os.Mkdir(dir, 0777) }

	var vals [][]float64
	if maxID >= RockstarShortMemoNum || maxID == -1 {
		file := path.Join(dir, fmt.Sprintf(rockstarMemoFile, snap))
		vals, err = readRockstar(
			file, -1, snap, nil, halo.ID, flag,
		)
		if err != nil { return nil, err }
	} else {
		file := path.Join(dir, fmt.Sprintf(rockstarShortMemoFile, snap))
		vals, err = readRockstar(
			file, RockstarShortMemoNum, snap, nil, halo.ID, flag,
		)
		if err != nil { return nil, err }
	}
	
	fids, ms := vals[0], vals[1]
	ids := make([]int, len(fids))
	for i := range ids { ids[i] = int(fids[i]) }

	if len(ids) < maxID {
		return nil, fmt.Errorf(
			"ID %d too large for snapshot %d", maxID, snap,
		)
	}

	sortRockstar(ids, ms)
	if maxID == -1 { return ids, nil }
	return ids[:maxID+1], nil
}

type massSet struct {
	ids []int
	ms []float64
}

func (set massSet) Len() int { return len(set.ids) }
// We're reverse sorting.
func (set massSet) Less(i, j int) bool { return set.ms[i] > set.ms[j] }
func (set massSet) Swap(i, j int) {
	set.ms[i], set.ms[j] = set.ms[j], set.ms[i]
	set.ids[i], set.ids[j] = set.ids[j], set.ids[i]
}

func sortRockstar(ids []int, ms []float64) {
	set := massSet{ ids, ms }
	sort.Sort(set)
}

// This function does fairly large heap allocations even when it doesn't need
// to. Consider passing it a buffer.
func ReadRockstar(
	snap int, ids []int, valFlags ...halo.Val,
) ([][]float64, error) {
	// Find binFile.
	memoDir, err := MemoDir()
	if err != nil { return nil, err }
	dir := path.Join(memoDir, rockstarMemoDir)
	if !PathExists(dir) { os.Mkdir(dir, 0777) }

	binFile := path.Join(dir, fmt.Sprintf(rockstarMemoFile, snap))
	shortBinFile := path.Join(dir, fmt.Sprintf(rockstarShortMemoFile, snap))

	// This wastes a read the first time it's called. You need to decide if you
	// care. (Answer: probably.)
	vals, err := readRockstar(
		shortBinFile, RockstarShortMemoNum, snap, ids, valFlags...,
	)
	if err == nil { return vals, err }
	return readRockstar(binFile, -1, snap, ids, valFlags...)
}

func readRockstar(
	binFile string, n, snap int, ids []int, valFlags ...halo.Val,
) ([][]float64, error) {
	// If binFile doesn't exist, create it.
	if !PathExists(binFile) {
		rockstarDir, err := RockstarDir()
		if err != nil { return nil, err }
		hlists, err := DirContents(rockstarDir)
		snapNum, err := SnapNum()
		if err != nil { return nil, err }
		negSnap := snapNum - snap
		snapIdx := len(hlists) - 1 - negSnap
		
		if err != nil { return nil, err }
		if n == -1 {
			err = halo.RockstarConvert(hlists[snapIdx], binFile)
			if err != nil { return nil, err }
		} else {
			err = halo.RockstarConvertTopN(hlists[snapIdx], binFile, n)
			if err != nil { return nil, err}
		}
	}

	// Get cosmo header.
	gtetFmt, err := GtetFmt()
	if err != nil { return nil, err }
	files, err := DirContents(fmt.Sprintf(gtetFmt, snap))
	hd := &io.SheetHeader{}
	err = io.ReadSheetHeaderAt(files[0], hd)

	rids, rvals, err := halo.ReadBinaryRockstarVals(
		binFile, &hd.Cosmo, valFlags...,
	)	
	if err != nil { return nil, err }
	
	// Select out only the IDs we want.
	if ids == nil { return rvals, nil }
	vals := make([][]float64, len(rvals))

	for i := range vals { vals[i] = make([]float64, len(ids)) }
	f := NewIntFinder(rids)
	for i, id := range ids {
		line, ok := f.Find(id)
		if !ok { return nil, fmt.Errorf("Could not find ID %d", id) }
		for vi := range vals { vals[vi][i] = rvals[vi][line] }
	}
	
	return vals, nil
}

// A quick generic wrapper for doing those one-to-one mappings I need to do so
// often. Written like this so the backend can be swapped out easily.
type IntFinder struct {
	m map[int]int
}

func NewIntFinder(rids []int) IntFinder {
	f := IntFinder{}
	f.m = make(map[int]int)
	for i, rid := range rids { f.m[rid] = i }
	return f
}

func (f IntFinder) Find(rid int) (int, bool) {
	line, ok := f.m[rid]
	return line, ok
}

func readHeadersFromSheet(snap int) ([]io.SheetHeader, []string, error) {
	gtetFmt, err := GtetFmt()
	if err != nil { return nil, nil, err }
	dir := fmt.Sprintf(gtetFmt, snap)
	files, err := DirContents(dir)
	if err != nil { return nil, nil, err }

	hds := make([]io.SheetHeader, len(files))
	for i := range files {
		err = io.ReadSheetHeaderAt(files[i], &hds[i])
		if err != nil { return nil, nil, err }
	}
	return hds, files, nil
}


func ReadHeaders(snap int) ([]io.SheetHeader, []string, error) {
	memoDir, err := MemoDir()
	if err != nil { return nil, nil, err }
	if _, err := os.Stat(memoDir); err != nil {
		return nil, nil, err
	}

	memoFile := path.Join(memoDir, fmt.Sprintf("hd_snap%d.dat", snap))

	if _, err := os.Stat(memoFile); err != nil {
		// File not written yet.
		hds, files, err := readHeadersFromSheet(snap)
		if err != nil { return nil, nil, err }
		
        f, err := os.Create(memoFile)
        if err != nil { return nil, nil, err }
        defer f.Close()
        binary.Write(f, binary.LittleEndian, hds)

		return hds, files, nil
	} else {
		// File exists: read from it instead.

		f, err := os.Open(memoFile)
        if err != nil { return nil, nil, err }
        defer f.Close()
		
		n, err := SheetNum(snap)
		if err != nil { return nil, nil, err }
		hds := make([]io.SheetHeader, n)
        binary.Read(f, binary.LittleEndian, hds) 

		gtetFmt, err := GtetFmt()
		if err != nil { return nil, nil, err }
		dir := fmt.Sprintf(gtetFmt, snap)
		files, err := DirContents(dir)
		if err != nil { return nil, nil, err }

		return hds, files, nil
	}
}
