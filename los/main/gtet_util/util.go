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
