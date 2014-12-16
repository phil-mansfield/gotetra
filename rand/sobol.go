package rand

import (
	"fmt"
	"log"
	"os"
)

const (
	MaxDim = uint32(6)
	MaxBit = uint32(30)
	maxSeqNum = 1<<MaxBit
	fac = 1.0 / maxSeqNum
)

var (
	initMDeg = [MaxDim]uint32{1,2,3,3,4,4}
	initIp = [MaxDim]uint32{0,1,1,2,1,4}
	initIv = [MaxBit*MaxDim]uint32{
		1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,
	}
)

type SobolSequence struct {
	seqNum uint32
	ix, mdeg, ip [MaxDim]uint32
	iv [MaxBit*MaxDim]uint32
	fac float64

	isInit bool
}
// See Press et al. 2007.
func NewSobolSequence() *SobolSequence {
	seq := &SobolSequence{}
	seq.Init()
	return seq
}

// See Press et al. 2007.
func (seq *SobolSequence) Init() {
	for i := range seq.ix { seq.ix[i] = 0 }
	if seq.isInit { return }
	seq.seqNum = 0

	seq.iv = initIv
	seq.ip = initIp
	seq.mdeg = initMDeg

	for k := uint32(0); k < MaxDim; k++ {
		for j := uint32(0); j < seq.mdeg[k]; j++ {
			seq.iv[MaxDim * j + k] <<= MaxBit - j - 1
		}

		deg := seq.mdeg[k]
		for j := deg; j < MaxBit; j++ {
			ipp := seq.ip[k]
			i := seq.iv[MaxDim * (j - deg) + k]
			i ^= (i >> deg)

			for l := deg - 1; l >= 1; l-- {
				if 1 & ipp == 1 { i ^= seq.iv[MaxDim * (j - 1) + k] }
				ipp >>= 1
			}

			seq.iv[MaxDim * j + k] = i
		}
	}

	seq.isInit = true
}

// See Press et al. 2007.
func (seq *SobolSequence) Next(dim int) []float64 {
	target := make([]float64, dim)
	seq.NextAt(target)
	return target
}

// NextAt is equivelent to Next, except the Sobol sequence is returned in-place.
func (seq *SobolSequence) NextAt(target []float64) {
	dim := uint32(len(target))
	if dim > MaxDim {
		log.Fatalf("Target dim %d is larger than MaxDim %d.\n", dim, MaxDim)
	} else if seq.seqNum >= maxSeqNum {
		log.Fatalf("Exceeded maximum seq num of %d for MaxBit %d.\n",
			maxSeqNum, MaxBit,
		)
	}

	seq.seqNum++

	zeroIdx := uint32(0)
	for zeroIdx = 0; zeroIdx < MaxBit; zeroIdx++ {
		if (seq.seqNum & (1<<zeroIdx)) == 0 { break }
	}

	im := zeroIdx * MaxDim
	for k := uint32(0); k < dim; k++ {
		seq.ix[k] ^= seq.iv[im + k]
		target[k] = float64(seq.ix[k]) * fac
	}
}

func min(x, y int) int {
	if x <= y { return x }
	return y
}

func main() {
	seq := NewSobolSequence()
	xs := []float64{0, 0}
	for i := 0; i < 1<<12; i++ {
		seq.NextAt(xs)
		fmt.Println(xs[0], xs[1])
	}
	os.Exit(0)
}
