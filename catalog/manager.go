package catalog

type cmpIdx struct {
	slice, p int32
}

type ParticleManager struct {
	ps   [][]Particle
	Locs map[int64]cmpIdx
	Size int64
}

func NewParticleManager() *ParticleManager {
	man := &ParticleManager{
		[][]Particle{},
		make(map[int64]cmpIdx),
		0,
	}

	return man
}

func (man *ParticleManager) Add(ps []Particle) {
	man.ps = append(man.ps, ps)

	for i := range ps {
		man.Locs[ps[i].Id] = cmpIdx{int32(len(man.ps) - 1), int32(i)}
	}

	man.Size += int64(len(ps))
}

func (man *ParticleManager) Get(id int64) *Particle {
	idx, ok := man.Locs[id]

	if !ok {
		return nil
	}

	return &man.ps[idx.slice][idx.p]
}
