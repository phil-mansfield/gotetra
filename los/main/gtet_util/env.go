package gtet_util

import (
	"fmt"
	"os"
)

func RockstarDir() (string, error) {
	str := os.Getenv("GTET_ROCKSTAR_DIR")
	if str == "" { return "", fmt.Errorf("GTET_ROCKSTAR_DIR not set.") }
	return str, nil
}

func TreeDir() (string, error) {
	str := os.Getenv("GTET_TREE_DIR")
	if str == "" { return "", fmt.Errorf("GTET_TREE_DIR not set.") }
	return str, nil
}

func MemoDir() (string, error) {
	str := os.Getenv("GTET_MEMO_DIR")
	if str == "" { return "", fmt.Errorf("GTET_MEMO_DIR not set.") }
	return str, nil
}

func GtetFmt() (string, error) {
	str := os.Getenv("GTET_FMT")
	if str == "" { return "", fmt.Errorf("GTET_FMT not set.") }
	// In theory, I could manually check that the format is valid.
	return str, nil
}

func SheetNum(snap int) (int, error) {
    gtetFmt, err := GtetFmt()
    if err != nil { return 0, err }
    dir := fmt.Sprintf(gtetFmt, snap)
    files, err := DirContents(dir)
    if err != nil { return 0, err }
    return len(files), nil
}
