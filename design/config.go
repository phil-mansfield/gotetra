package main

import (
	"log"
	"os"

	"code.google.com/p/gcfg"
)

type PhaseBoundsConfig struct {
	XOrigin []float64
}

type PhaseBoundsWrapper struct {
	PhaseBounds PhaseBoundsConfig
}

func main() {
	bc := PhaseBoundsWrapper{}
	if err := gcfg.ReadFileInto(&bc, os.Args[1]); err != nil {
		log.Fatal(err.Error())
	}

	log.Println(bc)
}
