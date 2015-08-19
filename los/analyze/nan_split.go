package analyze

/* I'm actually pretty darn happy with this implementation. */

import (
	"math"
)

type nanSplitParams struct {
	aux [][]float64
	valSets [][]float64
	auxSets [][][]float64
}
type internalNaNSplitOption func(*nanSplitParams)

// NaNSplitOption is an abstract data type which allows for the customization of
// calls to NaNSplit without cluttering the call signature in the common case.
// This works similarly to kwargs in other languages.
type NaNSplitOption internalNaNSplitOption

// Aux gives additional slices to NaNSplit which will be split at the same
// locations as the main input array.
func Aux(aux ...[]float64) NaNSplitOption {
	return func(p *nanSplitParams) { p.aux = aux }
}

// ValSets supplies NaNSplit with a pre-allocated valSets array.
func ValSets(valSets [][]float64) NaNSplitOption {
	return func(p *nanSplitParams) { p.valSets = valSets }
}

// AuxSets supplies NaNSplit with a pre-allocated valSets array.
func AuxSets(auxSets [][][]float64) NaNSplitOption {
	return func(p *nanSplitParams) { p.auxSets = auxSets }
}

func (p *nanSplitParams) loadOptions(opts []NaNSplitOption) {
	for _, opt := range opts { opt(p) }

	// valSets
	if p.valSets != nil {
		// Clear, so append can be used.
		p.valSets = p.valSets[0:0]
	}

	// auxSets
	if p.auxSets == nil {
		p.auxSets = make([][][]float64, len(p.aux))
		for i := range p.aux {
			p.auxSets[i] = [][]float64{}
		}
	} else {
		if len(p.auxSets) != len(p.aux) {
			panic("Length of supplied auxSets doesn't make length of aux.")
		}
		for i := range p.auxSets {
			// Clear, so append can be used.
			p.auxSets[i] = p.auxSets[i][0:0]
		}
	}
}

// NaNSplit splits a slice into all non-empty sub-slices which are separated by
// NaNs. Additional arrays may be passed to NaNSplit, which will be split at
// the same locations.
func NaNSplit(
	vals []float64, opts ...NaNSplitOption,
) (valSets [][]float64, auxSets [][][]float64) {
	params := new(nanSplitParams)
	params.loadOptions(opts)

	aux := params.aux
	valSets = params.valSets
	auxSets = params.auxSets

	for _, a := range aux {
		if len(vals) != len(a) {
			panic("All slices given to NaNSplit must be the same length.")
		}
	}
	
	rangeStart := 0
	midNaN := true
	for i, val := range vals {
		if math.IsNaN(val) {
			// Encountered end of valid range.
			if !midNaN {
				valSets = append(valSets, vals[rangeStart: i])
				for j := range aux {
					auxSets[j] = append(auxSets[j], aux[j][rangeStart: i])
				}
				midNaN = true
			}
		} else {
			// Encountered start of valid range.
			if midNaN {
				rangeStart = i
				midNaN = false
			}
		}
	}

	// The seqeunce ended on a valid range.
	if !midNaN {
		valSets = append(valSets, vals[rangeStart: len(vals)])
		for j := range aux {
			auxSets[j] = append(auxSets[j], aux[j][rangeStart: len(vals)])
		}
		midNaN = true
	}
	
	return valSets, auxSets
}
