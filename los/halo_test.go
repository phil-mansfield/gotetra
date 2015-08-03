package los

import (
	"runtime"
	"testing"

	"github.com/phil-mansfield/gotetra/los/geom"
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

func BenchmarkHaloProfilesClear(b *testing.B) {
	hp := new(HaloProfiles)
	hp.Init(0, 1, &geom.Vec{0, 0, 0}, 0, 1, 200, 100000)
	for i := 0; i < b.N; i++ { hp.Clear() }
}

func BenchmarkHaloProfilesParallelClear(b *testing.B) {
	hs := make([]HaloProfiles, runtime.NumCPU())
	for i := range hs { 
		hs[i].Init(0, 1, &geom.Vec{0, 0, 0}, 0, 1, 200, 100000)
	}
	for i := 0; i < b.N/len(hs); i++ { ParallelClearHaloProfiles(hs) }
}

func BenchmarkHaloProfilesAdd(b *testing.B) {
	hp1, hp2 := new(HaloProfiles), new(HaloProfiles)
	hp1.Init(0, 10, &geom.Vec{0, 0, 0}, 0, 1, 200, 1000)
	hp2.Init(0, 10, &geom.Vec{0, 0, 0}, 0, 1, 200, 1000)
	for i := 0; i < b.N; i++ { hp1.Add(hp2) }
}
