package rand

import (
	"math/rand"
)

type golangGenerator struct {
	r *rand.Rand
}

func (gen *golangGenerator) Init(seed uint64) {
	src := rand.NewSource(int64(seed))
	gen.r = rand.New(src)
}

func (gen *golangGenerator) Next() float64 {
	return gen.r.Float64()
}

func (gen *golangGenerator) NextSequence(target []float64) {
	for i := range target {
		target[i] = gen.r.Float64()
	}
}
