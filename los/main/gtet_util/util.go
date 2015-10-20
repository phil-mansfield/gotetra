package gtet_util

import (
	"fmt"
	"io/ioutil"
	"os"
	"path"

	"github.com/phil-mansfield/gotetra/render/io"
)

func PathExists(path string) bool {
	_, err := os.Stat(path)
	return err == nil
}

func DirContents(dir string) ([]string, error) {
	infos, err := ioutil.ReadDir(dir)
	if err != nil { return nil, err }
	
	files := make([]string, len(infos))
	for i := range infos {
		files[i] = path.Join(dir, infos[i].Name())
	}

	return files, nil
}

func Filter(xs []int, oks []bool) []int {
	n := 0
	for _, ok := range oks {
		if ok { n++ }
	}

	out := make([]int, 0, n)
	for i, x := range xs {
		if oks[i] { out = append(out, x) }
	}

	return out
}

func ReadSnapHeader(snap int) (*io.SheetHeader, error) {
	gtetFmt, err := GtetFmt()
	if err != nil { return nil, err }

	gtetDir := fmt.Sprintf(gtetFmt, snap)
	gtetFiles, err := DirContents(gtetDir)
	if err != nil { return nil, err }

	hd := &io.SheetHeader{}
	err = io.ReadSheetHeaderAt(gtetFiles[0], hd)
	if err != nil { return nil, err }
	return hd, nil
}
