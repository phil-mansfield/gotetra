package los

import (
	"testing"
)

func TestInRange(t *testing.T) {
	tw := float32(100)
	r := float32(10)
	table := []struct {
		x, low, width float32
		res bool
	} {
		{50, 40, 20, true},
		{20, 40, 20, false},
		{80, 40, 20, false},
		{45, 40, 20, true},
		{65, 40, 20, true},

		{5, 90, 20, true},
		{15, 90, 20, true},
		{25, 90, 20, false},
	}

	for i, test := range table {
		ir := inRange(test.x, r, test.low, test.width, tw)
		if ir != test.res {
			t.Errorf(
				"%d) inRange(%g, %g, %g, %g, %g) != %v",
				i + 1, test.x, r, test.low, test.width, tw, test.res,
			)
		}
	}
}
