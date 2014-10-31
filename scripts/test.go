package main

import (
	"fmt"

	"encoding/binary"
	"os"
	"time"

	"github.com/phil-mansfield/gotetra"
)

func main() {
	t0 := time.Now().UnixNano()

	_, ps :=  gotetra.ReadGadget(os.Args[1], binary.LittleEndian)

	t1 := time.Now().UnixNano()

	m := make(map[int64]int)
	for i, p := range ps {
		m[p.Id] = i
	}

	t2 := time.Now().UnixNano()
	

	fmt.Printf("Reading: %.0f ms, Mapping: %.0f ms\n",
		float64(t1 - t0) / 1e6, float64(t2 - t1) / 1e6)
}
