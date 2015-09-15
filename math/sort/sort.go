/*package sort provides functions for sorting and finding the median of
float64 slices without the overhead of Go's interfaces.
*/
package sort

const m = 5

// Reverse reverses a slice in place (and returns it for convenience).
func Reverse(xs []float64) []float64 {
	n1, n2 := len(xs) - 1, len(xs) / 2
	for i := 0; i < n2; i++ {
		xs[i], xs[n1 - i] =  xs[n1 - i], xs[i]
	}
	return xs
}
