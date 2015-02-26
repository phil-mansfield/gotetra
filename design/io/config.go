package main

import (
	"fmt"
	"log"
	"os"

	"code.google.com/p/gcfg"
)

type BallConfig struct {
	// Required
	X, Y, Z, Radius float64

	// Optional
	RadiusMultiplier float64
	Name string
}

func (ball *BallConfig) CheckInit(name string) error {
	if ball.Radius == 0 {
		return fmt.Errorf(
			"Need to specify a positive radius for Ball '%s'", name,
		)
	}

	ball.Name = name
	if ball.RadiusMultiplier == 0 {
		ball.RadiusMultiplier = 1
	}

	return nil
}

type BoxConfig struct {
	// Required
	X, Y, Z float64
	XWidth, YWidth, ZWidth float64

	// Optional
	Name string
}

func (box *BoxConfig) CheckInit(name string) error {
	if box.XWidth == 0 {
		return fmt.Errorf(
			"Need to specify a positive XWidth for Box '%s'", name,
		)
	} else if box.YWidth == 0 {
		return fmt.Errorf(
			"Need to specify a positive YWidth for Box '%s'", name,
		)
	} else if box.ZWidth == 0 {
		return fmt.Errorf(
			"Need to specify a positive ZWidth for Box '%s'", name,
		)
	}

	box.Name = name

	return nil
}

type BoundsConfig struct {
	Ball map[string]*BallConfig
	Box  map[string]*BoxConfig
}

func main() {
	if len(os.Args) != 2 {
		log.Fatal("Expects exactly one argument.")
	}

	bc := BoundsConfig{}

	err := gcfg.ReadFileInto(&bc, os.Args[1])
	if err != nil {
		log.Fatal(err.Error())
	}

	fmt.Println(bc)

	for name, ball := range bc.Ball {
		if err := ball.CheckInit(name); err != nil {
			log.Fatal(err.Error())
		}
	}

	for name, box := range bc.Box {
		if err := box.CheckInit(name); err != nil {
			log.Fatal(err.Error())
		}
	}
}
