package density

import (
	"fmt"
)

type Quantity int64
const (
	Density Quantity = iota
	DensityGradient
	Velocity
	VelocityDivergence
	VelocityCurl
	EndQuantity
)

func (q Quantity) String() string {
	if q <= 0 || q >= EndQuantity {
		panic(fmt.Sprintf("Value %d out of range for Quantity type.", q))
	}

	switch q {
	case Density:
		return "Density"
	case DensityGradient:
		return "DensityGradient"
	case Velocity:
		return "Velocity"
	case VelocityDivergence:
		return "VelocityDivergence"
	case VelocityCurl:
		return "VelocityCurl"
	}

	panic("Quantity.String() missing a switch clause.")
}

func QuantityFromString(str string) (q Quantity, ok bool) {
	switch str {
	case "Density":
		return Density, true
	case "DensityGradient":
		return DensityGradient, true
	case "Velocity":
		return Velocity, true
	case "VelocityDivergence":
		return VelocityDivergence, true
	case "VelocityCurl":
		return VelocityCurl, true
	}
	return 0, false
}

func (q Quantity) RequiresVelocity() bool {
	switch q {
	case Density, DensityGradient:
		return false
	case Velocity, VelocityDivergence, VelocityCurl:
		return true
	}
	panic(":3")
}

func (q Quantity) CanProject() bool {
	if q <= 0 || q >= EndQuantity {
		panic(fmt.Sprintf("Value %d out of range for Quantity type.", q))
	}

	switch q {
	case Density, Velocity:
		return true
	case DensityGradient, VelocityCurl, VelocityDivergence:
		return false
	}

	panic("Quantity.String() missing a switch clause.")
}
