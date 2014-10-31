package main

import (
	"encoding/binary"
	"io/ioutil"
	"math"
	"os"
	"path"
	"regexp"
	"strings"
	"strconv"
	"sort"
	
	"github.com/phil-mansfield/gotetra"
)

// allMatches returns a slice of all filenames which path the given filepath.
// The base of the file path may be a regular expression using '*' as a
//  wildcard.
func allMatches(reStr string) []string {
	dir := path.Dir(reStr)
	reBase := path.Base(reStr)
	re := regexp.MustCompile(strings.Replace(reBase, "*", "*.", -1))

	fs, err := ioutil.ReadDir(dir)
	if err != nil { panic(err) }

	matches := make([]string, 0)

	for _, f := range fs {
		if re.MatchString(f.Name()) {
			matches = append(matches, path.Join(dir, f.Name()))
		}
	}

	return matches
}

func posToIdx(h *gotetra.Header, pos []float64, width int) int {
	idxs := make([]int, 3)
	
	for i := 0; i < 3; i++ {
		idxs[i] = int(math.Floor(pos[i] / h.BoxSize))
	}

	return idxs[0] + idxs[1] * width + idxs[2] * width * width
}

type idxPaths struct {
	idxs []int
	paths []string
}
func (ip *idxPaths) Len() int { return len(ip.idxs) }
func (ip *idxPaths) Less(i, j int) bool { return ip.idxs[i] < ip.idxs[j] }
func (ip *idxPaths) Swap(i, j int) {
	ip.idxs[i], ip.idxs[j] = ip.idxs[j], ip.idxs[i]
	ip.paths[i], ip.paths[j] = ip.paths[j], ip.paths[i]
}

func main() {
	reStr := os.Args[1]

	boxCount, err := strconv.Atoi(os.Args[2])
	if err != nil { panic(err) }
	width := int(math.Floor(math.Pow(float64(boxCount), 1 / 3.0)))

	matches := allMatches(reStr)

	ip := &idxPaths{ make([]int, 0), make([]string, 0) }
	for _, path := range matches {
		h, pos := gotetra.ReadOnePositionGadget(path, binary.LittleEndian)
		idx := posToIdx(h, pos, width)
		ip.idxs = append(ip.idxs, idx)
		ip.paths = append(ip.paths, path)
	}

	sort.Sort(ip)
	for i := 0; i < ip.Len(); i++ {
		println(i, ip.idxs[i], ip.paths[i])
	}
}
