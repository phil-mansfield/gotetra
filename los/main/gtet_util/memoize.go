package gtet_util

import (
	"fmt"
	"os"
	"path"
	
	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/render/io"
)

const (
	rockstarMemoDir = "rockstar"
	rockstarMemoFile = "halo_%d.dat"
)

// This function does fairly large heap allocations even when it doesn't need
// to. Consider passing it a buffer.
func ReadRockstar(snap int, ids []int, valFlags ...halo.Val,) ([][]float64,error) {
	// Find binFile.
	memoDir, err := MemoDir()
	if err != nil { return nil, err }
	dir := path.Join(memoDir, rockstarMemoDir)
	if !PathExists(dir) { os.Mkdir(dir, os.ModeDir) }
	binFile := path.Join(dir, fmt.Sprintf(rockstarMemoFile, snap))

	// If binFile doesn't exist, create it.
	if !PathExists(binFile) {
		rockstarDir, err := RockstarDir()
		if err != nil { return nil, err }
		hlists, err := DirContents(rockstarDir)
		if err != nil { return nil, err }
		asciiFile := path.Join(rockstarDir, hlists[snap - 1])
		halo.RockstarConvert(asciiFile, binFile)
	}

	// Get cosmo header.
	gtetFmt, err := GtetFmt()
	if err != nil { return nil, err }
	files, err := DirContents(fmt.Sprintf(gtetFmt, snap))
	hd := &io.SheetHeader{}
	err = io.ReadSheetHeaderAt(files[0], hd)

	// THe actual read happens here :)
	rids, rvals, err := halo.ReadBinaryRockstarVals(
		binFile, &hd.Cosmo, valFlags...,
	)
	if err != nil { return nil, err }
	
	// Select out only the IDs we want.
	vals := make([][]float64, len(rvals))
	found := make([]bool, len(ids))
	for i := range vals { vals[i] = make([]float64, len(ids)) }
	// I think this looping strategy has better cache properties. But it will
	// End up doing more match checks than it should.
	for i, rid := range rids {
		for j, id := range ids {
			// IDs match. Copy over the values.
			if id == rid {
				for vi := range vals { vals[vi][j] = rvals[vi][i] }
				found[j] = true
				continue
			}
		}
	}

	for i := range found {
		if !found[i] {
			return nil, fmt.Errorf("")
		}
	}

	return vals, err
}
