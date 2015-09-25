package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/phil-mansfield/gotetra/los/tree"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
)

var usage = "gtet_tree halo_id_1 [halo]  or  gtet_ids ... | gtet_tree"

func main() {
	// Parse input

	flag.Parse()

	cmdIDs, err := parseCmdArgs(flag.Args())
	if err != nil { log.Fatal(err.Error()) }

	stdin, err := readStdin()
	if err != nil { log.Fatal(err.Error()) }
	stdinIDs, err := parseStdinArgs(stdin)
	if err != nil { log.Fatal(err.Error()) }

	// Calculate and print trees.
	inputIDs := append(cmdIDs, stdinIDs...)
	if err != nil { log.Fatal(err.Error()) }

	trees, err := treeFiles()
	if err != nil { log.Fatal(err.Error()) }
	idSets, snapSets, err := tree.HaloHistories(trees, inputIDs)
	if err != nil { log.Fatal(err.Error())}

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

	printIDs(ids, snaps)
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

func parseStdinArgs(args []string) ([]int, error) {
	IDs := make([]int, 0, len(args))
	for i := range args {
		rawTokens := strings.Split(args[i], " ")
		tokens := make([]string, 0, len(rawTokens))
		for _, tok := range rawTokens {
			if len(tok) != 0 { tokens = append(tokens, tok) }
		}
		// This should be impossible, but whatever.
		if len(tokens) == 0 { continue }
		id, err := strconv.Atoi(tokens[0])
		if id == -1 { continue }
		IDs = append(IDs, id)
		if err != nil {
			return nil, fmt.Errorf(
				"Argument %d of command line args cannot be parsed.", i,
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

func readStdin() ([]string, error) {
	bs, err := ioutil.ReadAll(os.Stdin)
	if err != nil {
		return nil, fmt.Errorf("Problem reading stdin: %s", err)
	}
	text := string(bs)
	lines := strings.Split(text, "\n")
	fullLines := make([]string, 0, len(lines))
	for _, line := range lines {
		if len(line) != 0 { fullLines = append(fullLines, line) }
	}
	return fullLines, nil
}

func printIDs(ids []int, snaps []int) {
	// Find the maximum width of each column.
	idWidth, snapWidth := 0, 0
	for i := range ids {
		iWidth := len(fmt.Sprintf("%d", ids[i]))
		sWidth := len(fmt.Sprintf("%d", snaps[i]))
		if iWidth > idWidth { idWidth = iWidth }
		if sWidth > snapWidth { snapWidth = sWidth }
	}

	rowFmt := fmt.Sprintf("%%%dd %%%dd\n", idWidth, snapWidth)
	for i := range ids { fmt.Printf(rowFmt, ids[i], snaps[i]) }
}
