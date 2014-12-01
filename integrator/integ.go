/*Everything in this pacakge *will* be done functionally, so help me God.
*/
package main

import (
	"fmt"
	"strings"
)

type Flag uint8

type Expression interface {
	fmt.Stringer
	fmt.GoStringer
	Contains(string) bool
	Subst(string, Expression) Expression
	Simplify(string) Expression
}

type TermMaker interface {
	Expression
	Term(string) Term
}

type Term struct {
	
}

func (t Term) Expression() Expression {
	panic("Not Yet Implemented")
}

type Num float64
type Prod []Expression
type Sum []Expression
type Heav struct { Ex Expression }

type Frac struct {
	Num, Den Expression
}

type Var struct { 
	Name string
	Exp int
}

// Typechecking.

var (
	_ Expression = *new(Num)
	_ Expression = *new(Var)
	_ Expression = *new(Prod)
	_ Expression = *new(Sum)
	_ Expression = *new(Frac)
	_ Expression = *new(Heav)

	_ TermMaker = *new(Num)
	_ TermMaker = *new(Var)
	_ TermMaker = *new(Prod)
)

/*
Copy as neccesary. 
func (n Num) 
func (v Var) 
func (p Prod)
func (s Sum) 
func (f Frac)
func (h Heav)
*/

// String
func (n Num)  String() string {
	return fmt.Sprintf("%g", n)
}

func (f Frac) String() string {
	return fmt.Sprintf("(%s / %s)", f.Num.String(), f.Den.String())
}

func (p Prod) String() string {
	strs := make([]string, len(p))
	for i := range p { strs[i] = p[i].String() }
	return fmt.Sprintf("(%s)", strings.Join(strs, " "))
}

func (s Sum) String() string {
	strs := make([]string, len(s))
	for i := range s { strs[i] = s[i].String() }
	return fmt.Sprintf("(%s)", strings.Join(strs, " + "))
}

func (v Var)  String() string {
	if v.Exp == 1 {
		return v.Name
	} else {
		return fmt.Sprintf("%s^%d", v.Name, v.Exp)
	}
}

func (h Heav) String() string {
	return fmt.Sprintf("H(%s)", h.Ex.String())
}

// GoString

func (n Num)  GoString() string { return n.String() }
func (v Var)  GoString() string { return v.String() }
func (p Prod) GoString() string { return p.String() }
func (s Sum)  GoString() string { return s.String() }
func (f Frac) GoString() string { return f.String() }
func (h Heav) GoString() string { return h.String() }

// Contains

func (n Num)  Contains(name string) bool { return false }
func (v Var)  Contains(name string) bool { return v.Name == name }

func (p Prod) Contains(name string) bool {
	for i := range p {
		if p[i].Contains(name) { return true }
	}
	return false
}

func (s Sum)  Contains(name string) bool {
	for i := range s {
		if s[i].Contains(name) { return true }
	}
	return false
}

func (f Frac) Contains(name string) bool {
	return f.Num.Contains(name) || f.Den.Contains(name)
}

func (h Heav) Contains(name string) bool { return h.Ex.Contains(name) }

// Subst
//
// Could save time by memoizing Contains() calls.

func (n Num)  Subst(name string, ex Expression) Expression { return n }

func (v Var)  Subst(name string, ex Expression) Expression {
	if v.Name == name {
		p := Prod(make([]Expression, v.Exp))
		for i := range p {
			p[i] = ex
		}
		return p
	} else {
		return v
	}
}

func (p Prod) Subst(name string, ex Expression) Expression {
	exs := Prod(make([]Expression, len(p)))
	for i := range p {
		exs[i] = p[i].Subst(name, ex)
	}
	return exs
}

func (s Sum)  Subst(name string, ex Expression) Expression {
	exs := Sum(make([]Expression, len(s)))
	for i := range s {
		exs[i] = s[i].Subst(name, ex)
	}
	return exs
}

func (f Frac) Subst(name string, ex Expression) Expression {
	return Frac{ f.Num.Subst(name, ex), f.Den.Subst(name, ex) }
}

func (h Heav) Subst(name string, ex Expression) Expression {
	return Heav{ h.Ex.Subst(name, ex) }
}

// Simplify

func (n Num)  Simplify(name string) Expression {
	return n
}

func (v Var)  Simplify(name string) Expression {
	if (name == "" || name == v.Name) && v.Exp == 0 {
		return Num(1)
	}
	return v
}

func (p Prod) Simplify(name string) Expression {
	return p.Term(name).Expression()
}

func (s Sum)  Simplify(name string) Expression {
	panic("Not Yet Implemented.")
}

func (f Frac) Simplify(name string) Expression {
	num, den := f.Num.Simplify(name), f.Num.Simplify(name)
	switch den := den.(type) {
	case Num:
		switch num := num.(type) {
		case Num:
			return Num(num / den)
		default:
			return Frac{num, Num(den)}
		}
	default:
		return Frac{num, den}
	}
	panic("Impossible.")
}

func (h Heav) Simplify(name string) Expression {
	return Heav{ h.Ex.Simplify(name) }
}

// Term

func (n Num)  Term(name string) Term {
	panic("Not Yet Implemented")
}

func (v Var)  Term(name string) Term {
	panic("Not Yet Implemented")
}

func (p Prod) Term(name string) Term {
	panic("Not Yet Implemented")
}

func main() {
	var ex Expression = Sum{ Num(1), Num(2), Prod{ Var{"x", 2}, Heav{Var{"y", 1}} } }
	ex = ex.Subst("y", Sum{Num(1), Var{"z", 3}})
	fmt.Println(ex)
	fmt.Println(ex.Contains("z"))
}
