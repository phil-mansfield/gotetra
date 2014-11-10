package main

import (
	"fmt"
	"os"
	"path"

	"github.com/phil-mansfield/gotetra/catalog"
)

func main() {
	for i, arg := range os.Args {
		if i == 0 {
			fmt.Printf("# %20s %3s %3s %3s %10s\n",
				"File Name", "X", "Y", "Z", "Count")
			continue
		}

		h := catalog.ReadHeader(arg)
		name := path.Base(arg)

		x := h.Idx % h.GridWidth
		y := (h.Idx / h.GridWidth) % h.GridWidth
		z := h.Idx / (h.GridWidth * h.GridWidth)

		fmt.Printf("  %20s %3d %3d %3d %10d\n", name, x, y, z, h.Count)
	}
}
