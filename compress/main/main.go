package main

import (
	"fmt"
	"log"

	"github.com/phil-mansfield/gotetra/compress"
	"github.com/phil-mansfield/gotetra/render/io"
)

const (
	dir = "/project/surph/mansfield/data/sheet_segments/Box_L0063_N1024_G0008_CBol/snapdir_100/"
)

func main() {
	base := compress.GetBase(dir)
	hd := &io.SheetHeader{}
	err := io.ReadSheetHeaderAt(fmt.Sprintf(base, 0, 0, 0), hd)
	if err != nil { log.Fatal(err.Error()) }

	vs, err := compress.ReadLine(base, hd, -1, 900, 900)
	if err != nil { log.Fatal(err.Error()) }
	compress.Normalize(vs, hd)

	dx := (hd.TotalWidth / float64(hd.CountWidth))
	for i, v := range vs {
		x0 := dx * float64(i)
		_, _ = x0, v
		fmt.Printf("%.5g %.5g %.5g %.5g\n", x0, v[0], v[1], v[2])
	}
}
