package main

import (
	"encoding/binary"
	"os"

	"github.com/phil-mansfield/gotetra"
)

func main() {
	gotetra.ReadGadget(os.Args[1], binary.BigEndian)
}
