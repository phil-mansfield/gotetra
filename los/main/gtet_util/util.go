package gtet_util

import (
	"io/ioutil"
	"os"
	"path"
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
