package main

import (
	"fmt"
	"log"

	"github.com/phil-mansfield/gotetra/compress"
	"github.com/phil-mansfield/gotetra/render/io"
)

const (
	dir = "/project/surph/phil-mansfield/data/sheet_segments"
)

func main() {
	base := compress.GetBase(dir)
	hd := &io.SheetHeader{}
	err := io.ReadSheetHeaderAt(fmt.Sprintf(base, 0, 0, 0), hd)
	if err != nil { log.Fatal(err.Error()) }

	vs, err := compress.ReadLine(base, hd, 100, 100, -1)
	if err != nil { log.Fatal(err.Error()) }
	compress.Normalize(vs, 2, hd)

	dx := (hd.TotalWidth / float64(hd.CountWidth))
	for i := range vs {
		x0 := dx * float64(i)
		fmt.Printf("%.5g %.5g %.5g %.5g\n", x0, vs[0], vs[1], vs[2])
	}
}
