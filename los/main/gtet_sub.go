package main

import (
	"flag"
	"log"

	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
)

type Params struct {
	VMaxCutoff float64
}

func main() {
	log.Println("gtet_sub")
	p := parseCmd()
	_ = p
	
	ids, snaps, _, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }

	sIDs, sSnaps, hIDs, err := util.SubIDs(ids, snaps)
	if err != nil { log.Fatal(err.Error()) }
	fIDs := make([]float64, len(hIDs))
	for i := range fIDs { fIDs[i] = float64(hIDs[i]) }
	util.PrintCols(sIDs, sSnaps, fIDs)
}

func parseCmd() *Params {
	p := &Params{}
	// Good choices: 50, 100, 170, 320 km/s for L = 62.5, 125, 250, 500 Mpc/h.
	flag.Float64Var(&p.VMaxCutoff, "VMaxCutoff", 0,
		"Ignore subhalos with VMax below this number. I'm not going to " +
			"force you to, but you should seriously pay attention to this.")
	flag.Parse()
	return p
}
