package main

import (
	"fmt"

	"encoding/binary"
	"os"

	"github.com/phil-mansfield/gotetra"
)

func main() {
	ps, hd :=  gotetra.ReadGadget(os.Args[1], binary.LittleEndian)
	fmt.Println(hd)
	fmt.Println(ps[0])
	fmt.Println(ps[1])
}
