package gtet_util

import (
	"path"
	"fmt"
	"io/ioutil"
	"os"
)

// RockstarDir returns the environment variable $GTET_ROCKSTAR_DIR and
// returns an error if it has not been set.
func RockstarDir() (string, error) {
	str := os.Getenv("GTET_ROCKSTAR_DIR")
	if str == "" { return "", fmt.Errorf("GTET_ROCKSTAR_DIR not set.") }
	return str, nil
}

// TreeDir returns the environment variable $GTET_TREE_DIR and
// returns an error if it has not been set.
func TreeDir() (string, error) {
	str := os.Getenv("GTET_TREE_DIR")
	if str == "" { return "", fmt.Errorf("GTET_TREE_DIR not set.") }
	return str, nil
}

// MemoDir returns the environment variable $GTET_MEMO_DIR and
// returns an error if it has not been set.
func MemoDir() (string, error) {
	str := os.Getenv("GTET_MEMO_DIR")
	if str == "" { return "", fmt.Errorf("GTET_MEMO_DIR not set.") }
	return str, nil
}

// GtetFmt returns the environment variable $GTET_FMT_DIR and
// returns an error if it has not been set.
func GtetFmt() (string, error) {
	str := os.Getenv("GTET_FMT")
	if str == "" { return "", fmt.Errorf("GTET_FMT not set.") }
	// In theory, I could manually check that the format is valid.
	return str, nil
}

// SheetNum returns the number of sheet segments used by given snapshot.
// An error is returned on an I/O error.
func SheetNum(snap int) (int, error) {
    gtetFmt, err := GtetFmt()
    if err != nil { return 0, err }
    dir := fmt.Sprintf(gtetFmt, snap)
    files, err := DirContents(dir)
    if err != nil { return 0, err }
    return len(files), nil
}

// SnapNum returns the largest snapshot which can be indexed by $GTET_FMT.
// An error is returned on an I/O error.
func SnapNum() (int, error) {
	gtetFmt, err := GtetFmt()
	if err != nil { return 0, err }
	fmtDir := path.Dir(gtetFmt)
	fmtName := path.Base(gtetFmt)
	infos, err := ioutil.ReadDir(fmtDir)

	max := -1
	for i := range infos {
		n := 0
		_, err := fmt.Sscanf(infos[i].Name(), fmtName, &n)
		if err != nil { continue }
		if max < n { max = n }
	}
	
	if max == -1 { return -1, fmt.Errorf("No gtet files in $(GTET_FMT)") }
	return max, nil
}

// SnapOffset returns difference between the number of snapshots in $GTET_FMT
// and the number of halo catalogs in GTET_TREE_DIR. Useful for when only
// the last couple of snapshots are the only ones converted into sheet segments.
// An error is returned on I/O error.
func SnapOffset() (int, error) {
	snapNum, err := SnapNum()
	if err != nil { return -1, err }
	rockstarDir, err := RockstarDir()
	if err != nil { return -1, err }
	hlists, err := DirContents(rockstarDir)
	if err != nil { return -1, err }
	if snapNum - len(hlists) < 0 {
		return -1, fmt.Errorf("Fewer particle snapshots than halo lists.")
	}
	return snapNum - len(hlists), nil
}
