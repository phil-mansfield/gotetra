open Printf

type expression = Var of (string * int)
                  | Num of float
                  | Sum of expression array
                  | Prod of expression array
                  | Frac of (expression * expression)

let eps = 0.0000001

let print_ex (ex : expression) : unit =
  let rec print_ex' (ex' : expression) : unit =
    match ex' with
      Num x -> printf "%g" x;
    | Var (name, exp) -> printf "%s^%d" name exp;
    | Sum exs -> 
      if Array.length exs <= 2 then
        raise (failwith "Impossible Sum encountered.")
      else ();
      printf "(";
      for i = 0 to Array.length exs - 2 do
        print_ex' (exs.(i));
        printf " + ";
      done;
      print_ex' (exs.(Array.length exs - 1));
      printf ")";
    | Prod exs ->
      if Array.length exs <= 2 then
        raise (failwith "Impossible Prod encountered.")
      else ();
      printf "("; 
      for i = 0 to Array.length exs - 2 do
        print_ex' (exs.(i));
        printf " * ";
      done;
      print_ex' (exs.(Array.length exs - 1));
      printf ")";
    | Frac (num, den) ->
      printf "(";
      print_ex' num;
      printf "/";
      print_ex' den;
      printf ")"
        
  in
  print_ex' ex;
  print_newline ();
;; 

let array_filter pred xs =
  Array.of_list (List.filter pred (Array.to_list xs))
;;

let sum_array (ex : expression) : expression array =
  match ex with
    Sum exs -> exs
  | _ -> [| ex |]
;;

let concat_sums (exs : expression array) : expression array =
  Array.concat (Array.to_list (Array.map sum_array exs))
;;

let prod_array (ex : expression) : expression array =
  match ex with 
    Prod exs -> exs
  | Sum _ -> raise (failwith "I found a Sum in prod_array call.")
  | _ -> [| ex |]
;;

let sum_prod_array (ex : expression) : expression array array =
  match ex with
    Sum exs -> Array.map prod_array exs
  | _ -> [| prod_array ex |]

let interleave (xs : 'a array) (ys : 'b array) : ('a * 'b) array = 
  let w, h = Array.length xs, Array.length ys in
  Array.init (w * h) (fun i -> (xs.(i mod w), ys.(i / w)))
;;

(* ex1 and ex2 should alreay be expanded. *)
let rec foil ((ex1, ex2) : expression * expression) : expression =
  match (ex1, ex2) with
    (Sum _, _) | (_, Sum _) -> 
      let p_exs1, p_exs2 = sum_prod_array ex1, sum_prod_array ex2 in
      Sum (Array.map foil (Array.map (fun (x, y) -> (Prod x, Prod y))
                             (interleave p_exs1 p_exs2)))
  | (_, _) -> Prod(Array.append (prod_array ex1) (prod_array ex2))
;;

let is_frac (ex : expression) = 
  match ex with
    Frac _ -> true
  | _ -> false
;;

let split_fracs (ex : expression) :
    (expression array * expression array) =
  let exs = match ex with
      Prod exs -> exs
    | Sum exs -> exs
    | _ -> [| ex |]
  in
  let fracs = array_filter is_frac exs in
  let non_fracs = array_filter (fun ex -> not (is_frac ex)) exs in
  (fracs, non_fracs)
;;

let rec expand (ex : expression) : expression =
  match ex with
    Sum exs ->
      let exs' = Array.map expand exs in
      let (fracs, non_fracs) = split_fracs exs' in
      if Array.length fracs = 0 then
        Sum (concat_sums non_fracs)
      else
      Sum (concat_sums (Array.map expand exs))
  | Prod exs -> 
    let exs' = Array.fold_left (fun x y -> foil (x, y)) (Prod [||]) 
      (Array.map expand exs) in
    let fracs, non_fracs = split_fracs exs' in
    if Array.length fracs = 0 then Prod non_fracs
    else raise (failwith "NYI")
  | Num _ -> ex
  | Var _ -> ex
  | Frac (num, den) ->
    let num', den' = expand num, expand den in
    match (num', den') with
      (Frac (nn, nd), Frac (dn, dd)) ->
        Frac (expand (Prod [| nn; dd |]), expand (Prod [| nd; dn |]))
    | (Frac (nn, nd), d) -> Frac (nn, expand (Prod [| nd; d |]))
    | (n, Frac (dn, dd)) -> Frac (expand (Prod [| n; dd |]), dn)
    | (n, d) -> Frac(n, d)
;;

type term = { mutable num : float; vars : (string, int) Hashtbl.t }

let bindings (tbl : ('a, 'b) Hashtbl.t) : ('a * 'b) array =
  Array.of_list (Hashtbl.fold (fun k v acc -> (k, v) :: acc ) tbl [])
;;

(* does not preserve expand properties, obviously *)
let unpack_term (t : term) : expression = 
  let vars = bindings t.vars in
  match (t.num, Array.length vars) with
    (0.0, _) -> Num 0.0
  | (num, 0) ->
    Num num
  | (1.0, _) ->
    let vars' = (array_filter (fun (_, e) -> e != 0 ) vars) in
    let exs = (Array.map (fun (n, e) -> Var (n, e)) vars') in
    if Array.length exs = 0 then Num 1.0
    else if Array.length exs = 1 then exs.(0)
    else Prod exs
  | (num, _) -> 
    let exs = Array.map (fun (n, e) -> Var (n, e)) vars in
    if Array.length exs = 0 then Num num
    else Prod (Array.append ([| Num num |]) exs)
;;

(* assumes expand has already been called *)
let pack_term (ex : expression) : term =
  match ex with
    Num num -> {num = num; vars = Hashtbl.create 0}
  | Var (name, exp) -> 
    let term = {num = 1.0; vars = Hashtbl.create 1} in
    Hashtbl.add term.vars name exp;
    term
  | Prod exps ->
    let term = {num = 1.0; vars = Hashtbl.create (Array.length exps)} in
    for i=0 to Array.length exps - 1 do
      match exps.(i) with
        Num num -> term.num <- term.num *. num;
      | Var (name, exp) ->
        if Hashtbl.mem term.vars name then
          Hashtbl.replace term.vars name (exp + Hashtbl.find term.vars name)
        else Hashtbl.add term.vars name exp
      | _ -> raise (failwith "I found a compound expression in pack_term.")
    done;
    term
  | Sum _ -> raise (failwith "I cannot pack a Sum expression into a term.")
  | Frac _ -> raise (failwith "I cannot pack a Frac expression into a term.")
;;

(* Could speed this up by using a helper function in add_term. *)
let vars_eq (t1 : term) (t2 : term) : bool =
  if Hashtbl.length t1.vars != Hashtbl.length t2.vars then
    false
  else
    let eq = ref true in
    let b2 = bindings t2.vars in
    for i = 0 to Array.length b2 - 1 do
      let n, e = b2.(i) in
      eq := !eq && (Hashtbl.mem t1.vars n && (Hashtbl.find t1.vars n) = e)
    done;
    !eq
;;

(* I hate you, OCaml. *)
let update_num (t : term) (num : float) : term =
  t.num <- t.num +. num; t
;;

let rec add_term (ts : term list) (t : term) : term list =
  match ts with
    [] -> [t]
  | head :: tail ->
    if vars_eq head t then 
      if (update_num head t.num).num = 0.0 then
        tail
      else
        head :: tail
    else
      head :: add_term tail t
;;

let is_nonzero (ex : expression) : bool =
  match ex with
    Num num -> (num > eps) || (num < -.eps)
  | _ -> true
;;

(* assumes expand has already been called. *)
let rec simplify (ex : expression) : expression =
  match ex with 
  | Num _ -> ex 
  | Var _ -> ex
  | Prod _ -> unpack_term (pack_term ex)
  | Sum exs ->
    let terms = Array.map pack_term exs in
    let c_terms = Array.of_list (Array.fold_left add_term [] terms) in
    let exs = array_filter is_nonzero (Array.map unpack_term c_terms) in
    if Array.length exs = 0 then Num 0.0
    else if Array.length exs = 1 then exs.(0)
    else Sum exs
  | Frac (num, den) -> Frac (simplify num, simplify den)
;;

let rec subst (sub_name : string) (sub_ex : expression) (ex : expression) = 
  match ex with
    Num num -> ex
  | Var (name, exp) ->if name = sub_name then
      Prod (Array.make exp sub_ex) else ex
  | Prod exs -> Prod (Array.map (subst sub_name sub_ex) exs)
  | Sum exs -> Sum (Array.map (subst sub_name sub_ex) exs)
  | Frac (num, den) -> Frac (subst sub_name sub_ex num,
                             subst sub_name sub_ex den)
;;

let rec contains_var (var_name : string) (ex : expression) : bool =
  match ex with
  | Num _ -> false
  | Var (name, _) -> var_name = name
  | Prod exs ->
    Array.fold_left (||) false  (Array.map (contains_var var_name) exs)
  | Sum exs ->
    Array.fold_left (||) false  (Array.map (contains_var var_name) exs)
  | Frac (num, den) ->
    contains_var var_name num || contains_var var_name den
;;

let rec indef_integ (var_name : string) (ex : expression) : expression =
  match ex with
    Num _ ->  Prod [| ex; Var (var_name, 1)|]
  | Var (name, exp) -> 
    if name = var_name then
      Prod [| Num (1.0 /. ((float_of_int exp) +. 1.0));
              Var (name, exp + 1)|]
    else Prod [| Var (var_name, 1); ex |]
  | Prod _ -> 
    let term = pack_term ex in
    if Hashtbl.mem term.vars var_name then
      let exp = Hashtbl.find term.vars var_name in
      Prod [| Frac (Num 1.0, Num (1.0 +.  (float_of_int exp))); ex;
              Var (var_name, 1) |]
    else Prod [| ex; Var (var_name, 1)|]
  | Sum exs -> Sum (Array.map (indef_integ var_name) exs)
  | Frac (num, den) -> 
    if contains_var var_name den then
      raise (failwith "I can't evalue rational integrals.")
    else Frac (indef_integ var_name num, den)
;;

let integ (lo : expression) (hi : expression)
    (var_name : string) (ex : expression) : expression =
  let indef_ex = indef_integ var_name ex  in
  Sum [| subst var_name hi indef_ex; 
         Prod [| Num (-1.0);  subst var_name lo indef_ex|] |]
;;
 

(* So, remember: expand -> simplify, subst whenever, expand -> integrate *)
let main () =
  let ex = Prod [| Num 3.0; Var ("x", 1) |] in
  let indef_ex = simplify (expand (indef_integ "x" ex)) in
  let def_ex = simplify (expand (integ (Num (-2.0)) (Num 2.0) "x" ex)) in
  printf "Start with: ";
  print_ex ex;
  printf "Indefinite integral: ";
  print_ex indef_ex;
  printf "Definite integral: ";
  print_ex def_ex;
;;

main ();;
