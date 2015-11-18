package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"path"
	"strconv"

	"github.com/phil-mansfield/gotetra/los/tree"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
)

var usage = "gtet_tree halo_id_1 [halo]  or  gtet_ids ... | gtet_tree"

func main() {
	// Parse input

	flag.Parse()
	log.Println("gtet_tree")
	
	cmdIDs, err := parseCmdArgs(flag.Args())
	if err != nil { log.Fatal(err.Error()) }
	stdinIDs, _, _, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }

	// Calculate and print trees.
	inputIDs := append(cmdIDs, stdinIDs...)
	if err != nil { log.Fatal(err.Error()) }

	trees, err := treeFiles()
	if err != nil { log.Fatal(err.Error()) }
	snapOffset, err := util.SnapOffset()
	if err != nil { log.Fatal(err.Error()) }
	idSets, snapSets, err := tree.HaloHistories(trees, inputIDs, snapOffset)
	if err != nil { log.Fatal(err.Error()) }

	ids, snaps := []int{}, []int{}
	for i := range idSets {
		ids = append(ids, idSets[i]...)
		snaps = append(snaps, snapSets[i]...)
		// Sentinels:
		if i != len(idSets) - 1 {
			ids = append(ids, -1)
			snaps = append(snaps, -1)
		}
	}

	util.PrintCols(ids, snaps)
}

func parseCmdArgs(args []string) ([]int, error) {
	IDs := make([]int, len(args))
	var err error
	for i := range IDs {
		IDs[i], err = strconv.Atoi(args[i])
		if err != nil {
			return nil, fmt.Errorf(
				"Argument %d of command line args cannot be parsed.",
			)
		}
	}
	return IDs, nil
}

func treeFiles() ([]string, error) {
	treeDir, err := util.TreeDir()
	if err != nil { return nil, err }
	infos, err := ioutil.ReadDir(treeDir)
	if err != nil { return nil, err }

	names := []string{}
	for _, info := range infos {
		name := info.Name()
		n := len(name)
		// This is pretty hacky.
		if n > 4 && name[:5] == "tree_" && name[n-4:] == ".dat" {
			names = append(names, path.Join(treeDir, name))
		}
	}
	return names, nil
}
