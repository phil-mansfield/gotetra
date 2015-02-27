package io

import (
	"fmt"

	"code.google.com/p/gcfg"
)

type BallConfig struct {
	// Required
	X, Y, Z, Radius float64

	// Optional
	RadiusMultiplier float64
	Name string
}

func (ball *BallConfig) CheckInit(name string, totalWidth float64) error {
	if ball.Radius == 0 {
		return fmt.Errorf(
			"Need to specify a positive radius for Ball '%s'.", name,
		)
	}

	if ball.X >= totalWidth || ball.X < 0 {
		return fmt.Errorf(
			"X center of Ball '%s' must be in range [0, %g), but is %g",
			name, totalWidth, ball.X,
		)
	} else if ball.Y >= totalWidth || ball.Y < 0 {
		return fmt.Errorf(
			"Y center of Ball '%s' must be in range [0, %g), but is %g",
			name, totalWidth, ball.Y,
		)
	} else if ball.Z >= totalWidth || ball.Z < 0 {
		return fmt.Errorf(
			"Z center of Ball '%s' must be in range [0, %g), but is %g",
			name, totalWidth, ball.Z,
		)
	}

	ball.Name = name
	if ball.RadiusMultiplier == 0 {
		ball.RadiusMultiplier = 1
	} else if ball.RadiusMultiplier < 0 {
		return fmt.Errorf(
			"Ball '%s' given a negative radius multiplier, %g.",
			name, ball.RadiusMultiplier,
		)
	}

	return nil
}

func (ball *BallConfig) Box(totalWidth float64) *BoxConfig {
	box := &BoxConfig{}
	rad := ball.Radius * ball.RadiusMultiplier

	box.XWidth, box.YWidth, box.ZWidth = 2 * rad, 2 * rad, 2 * rad

	if ball.X > rad {
		box.X = ball.X - rad
	} else {
		box.X = ball.X - rad + totalWidth
	}

	if ball.Y > rad {
		box.Y = ball.Y - rad
	} else {
		box.Y = ball.Y - rad + totalWidth
	}

	if ball.Z > rad {
		box.Z = ball.Z - rad
	} else {
		box.Z = ball.Z - rad + totalWidth
	}

	box.Name = ball.Name
	
	return box
}

type BoxConfig struct {
	// Required
	X, Y, Z float64
	XWidth, YWidth, ZWidth float64

	// Optional
	Name string
}

func (box *BoxConfig) CheckInit(name string, totalWidth float64) error {
	if box.XWidth <= 0 {
		return fmt.Errorf(
			"Need to specify a positive XWidth for Box '%s'", name,
		)
	} else if box.YWidth <= 0 {
		return fmt.Errorf(
			"Need to specify a positive YWidth for Box '%s'", name,
		)
	} else if box.ZWidth <= 0 {
		return fmt.Errorf(
			"Need to specify a positive ZWidth for Box '%s'", name,
		)
	}

	if box.X >= totalWidth || box.X < 0 {
		return fmt.Errorf(
			"X origin of Box '%s' must be in range [0, %g), but is %g",
			name, totalWidth, box.X,
		)
	} else if box.Y >= totalWidth || box.Y < 0 {
		return fmt.Errorf(
			"Y origin of Box '%s' must be in range [0, %g), but is %g",
			name, totalWidth, box.Y,
		)
	} else if box.Z >= totalWidth || box.Z < 0 {
		return fmt.Errorf(
			"Z origin of Box '%s' must be in range [0, %g), but is %g",
			name, totalWidth, box.Z,
		)
	}

	box.Name = name

	return nil
}

type BoundsConfig struct {
	Ball map[string]*BallConfig
	Box  map[string]*BoxConfig
}

func ReadBoundsConfig(fname string, totalWidth float64) ([]BoxConfig, error) {
	bc := BoundsConfig{}

	if err := gcfg.ReadFileInto(&bc, fname); err != nil {
		return nil, err
	}

	boxes := []BoxConfig{}
	for name, ball := range bc.Ball {
		if err := ball.CheckInit(name, totalWidth); err != nil {
			return nil, err
		}
		boxes = append(boxes, *ball.Box(totalWidth))
	}
	for name, box := range bc.Box {
		if err := box.CheckInit(name, totalWidth); err != nil {
			return nil, err
		}
		boxes = append(boxes, *box)
	}

	return boxes, nil
}
