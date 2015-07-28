package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path"
	"sort"

	"github.com/phil-mansfield/table"
	"github.com/phil-mansfield/gotetra/render/io"	
	"github.com/phil-mansfield/gotetra/render/halo"

//	"github.com/phil-mansfield/gotetra/los"
)

const (
	rType = halo.R200m
)

func main() {
	if len(os.Args) != 2 {
		log.Fatalf("Usage: $ %s input_dir halo_file", os.Args[0])
	}

	dirName := os.Args[1]
	haloFileName := os.Args[2]

	files, err := fileNames(dirName)
	if err != nil { log.Fatal(err.Error()) }
	hds := make([]io.SheetHeader, len(files))
	for i := range files {
		err = io.ReadSheetHeaderAt(files[i], &hds[i])
		if err != nil { log.Fatal(err.Error()) }
	}

	xs, ys, zs, ms, rs, err := readHalos(haloFileName, &hds[0].Cosmo)
	if err != nil { log.Fatal(err.Error()) }
	
	fmt.Println("xs", xs[:5])
	fmt.Println("ys", ys[:5])
	fmt.Println("zs", zs[:5])
	fmt.Println("ms", ms[:5])
	fmt.Println("rs", rs[:5])
	fmt.Println("# of headers:", len(hds))
}

func fileNames(dirName string) ([]string, error) {
	infos, err := ioutil.ReadDir(dirName)
	if err != nil { return nil, err }

	files := make([]string, len(infos))
	for i := range infos {
		files[i] = path.Join(dirName, infos[i].Name())
	}
	return files, nil
}

type halos struct {
	xs, ys, zs, ms, rs []float64
}

func (hs *halos) Len() int { return len(hs.rs) }
func (hs *halos) Less(i, j int) bool { return hs.rs[i] < hs.rs[j] }
func (hs *halos) Swap(i, j int) {
	hs.rs[i], hs.rs[j] = hs.rs[j], hs.rs[i]
	hs.ms[i], hs.ms[j] = hs.ms[j], hs.ms[i]
	hs.xs[i], hs.xs[j] = hs.xs[j], hs.xs[i]
	hs.ys[i], hs.ys[j] = hs.ys[j], hs.ys[i]
	hs.zs[i], hs.zs[j] = hs.zs[j], hs.zs[i]
}

func readHalos(
	file string, cosmo *io.CosmologyHeader,
) (xs, ys, zs, ms, rs []float64, err error) {
	rCol := rType.RockstarColumn()
	xCol, yCol, zCol := 17, 18, 19
	
	colIdxs := []int{ xCol, yCol, zCol, rCol }
	cols, err := table.ReadTable(file, colIdxs, nil)
	if err != nil { return nil, nil, nil, nil, nil, err }
	
	xs, ys, zs = cols[0], cols[1], cols[2]
	if rType.RockstarMass() {
		ms = cols[3]
		rs = make([]float64, len(ms))
		rType.Radius(cosmo, ms, rs)
	} else {
		rs = cols[3]
		ms = make([]float64, len(rs))
		for i := range rs { rs[i] /= 1000 } // kpc -> Mpc
		rType.Mass(cosmo, rs, ms)
	}

	sort.Sort(sort.Reverse(&halos{ xs, ys, zs, ms, rs }))
	return xs, ys, zs, ms, rs, nil
}
