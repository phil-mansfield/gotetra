open Printf

(* This file is one giant cesspool. *)

type expression = Var of (string * int)
                  | Num of float
                  | Sum of expression array
                  | Prod of expression array
                  | Frac of (expression * expression)
                  | Heaviside of expression

let eps = 0.0000001

let print_type (ex : expression) : unit =
  print_endline (match ex with
    Num _ -> "Num"
  | Var _ -> "Var"
  | Prod _ -> "Prod"
  | Sum _ -> "Sum"
  | Frac _ -> "Frac"
  | Heaviside _ -> "Heaviside")
;;

let x_max (x : int) (y : int) : int = if x > y then x else y ;;
let xs_max (xs : int array) : int =
  Array.fold_left x_max xs.(0) xs
;;

let rec max_depth (ex: expression) : int = 
  1 + (match ex with
    Num _ -> 0
  | Var _ -> 0
  | Prod xs -> xs_max (Array.map max_depth xs)
  | Sum xs -> xs_max (Array.map max_depth xs)
  | Frac (n, d) -> x_max (max_depth n) (max_depth d)
  | Heaviside x -> max_depth x)
;;

let print_ex (ex : expression) : unit =
  let rec print_ex' (ex' : expression) : unit =
    match ex' with
      Num x -> printf "%.1g" x;
    | Var (name, exp) ->
      if exp != 1 then
        printf "math.Pow(%s, %d)" name exp
      else 
        printf "%s" name
    | Sum exs -> 
      if Array.length exs < 2 then
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
      if Array.length exs < 2 then
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
      printf " / ";
      print_ex' den;
      printf ")"
    | Heaviside ex ->
      printf "H(";
      print_ex' ex;
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

let rec foil ((ex1, ex2) : (expression * expression)) : expression =
  match (ex1, ex2) with
    (Sum _, _) | (_, Sum _) -> 
      let p_exs1, p_exs2 = sum_prod_array ex1, sum_prod_array ex2 in
      Sum (Array.map foil (Array.map (fun (x, y) -> (Prod x, Prod y))
                             (interleave p_exs1 p_exs2)))
  | (_, _) -> 
    Prod(Array.append (prod_array ex1) (prod_array ex2))
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

let get_denom (ex : expression) : expression =
  match ex with
    Frac (_, den) -> den
  | _ -> raise (failwith "I can only get_denom on a Frac.")
;;

let get_num (ex : expression) : expression =
  match ex with
    Frac (num, _) -> num
  | _ -> raise (failwith "I can only get_num on a Frac.")
;;

let total_denom (fracs : expression array) (ignore : int option) =
  match ignore with
    None -> 
      let d_exs = Array.map get_denom fracs in
      if Array.length d_exs = 1 then d_exs.(0) else Prod d_exs
  | Some ig_idx ->
    let d_exs = Array.make (Array.length fracs - 1) (Num 0.0) in
    for r_idx=0 to Array.length fracs - 1 do
      let w_idx = if r_idx > ig_idx then r_idx - 1 else r_idx in
      if r_idx != ig_idx then d_exs.(w_idx) <- get_denom fracs.(r_idx) else ();
    done;
    if Array.length d_exs = 0 then Num 1.0
    else if Array.length d_exs = 1 then d_exs.(0)
    else Prod d_exs
;;

type term = { mutable num : float; vars : (string, int) Hashtbl.t;
              mutable h_exs : expression list }

let bindings (tbl : ('a, 'b) Hashtbl.t) : ('a * 'b) array =
  Array.of_list (Hashtbl.fold (fun k v acc -> (k, v) :: acc ) tbl [])
;;

let unpack_term (t : term) : expression = 
  let vars = bindings t.vars in
  match (t.num, Array.length vars + List.length t.h_exs) with
    (0.0, _) -> Num 0.0
  | (num, 0) ->
    Num num
  | (1.0, _) ->
    let vars' = (array_filter (fun (_, e) -> e != 0 ) vars) in
    let exs = Array.append (Array.map (fun (n, e) -> Var (n, e)) vars')
      (Array.of_list (List.map (fun x -> Heaviside x) t.h_exs)) in
    if Array.length exs = 0 then Num 1.0
    else if Array.length exs = 1 then exs.(0)
    else Prod exs
  | (num, _) -> 
    let vars' = (array_filter (fun (_, e) -> e != 0 ) vars) in
    let exs = Array.append (Array.map (fun (n, e) -> Var (n, e)) vars')
      (Array.of_list (List.map (fun x -> Heaviside x) t.h_exs)) in
    if Array.length exs = 0 then Num num
    else Prod (Array.append ([| Num num |]) exs)
;;

let pack_term (ex : expression) : term =
  match ex with
    Num num -> {num = num; vars = Hashtbl.create 0; h_exs = []}
  | Var (name, exp) -> 
    let term = {num = 1.0; vars = Hashtbl.create 1; h_exs = []} in
    Hashtbl.add term.vars name exp;
    term
  | Prod exps ->
    let term = {num = 1.0; vars = Hashtbl.create (Array.length exps);
                h_exs = []} in
    for i=0 to Array.length exps - 1 do
      match exps.(i) with
        Num num -> term.num <- term.num *. num;
      | Var (name, exp) ->
        if Hashtbl.mem term.vars name then
          Hashtbl.replace term.vars name (exp + Hashtbl.find term.vars name)
        else Hashtbl.add term.vars name exp
      | Heaviside h_ex -> term.h_exs <- h_ex :: term.h_exs
      | Prod _ | Sum _ | Frac _ ->
        raise (failwith "I found a compound expression in pack_term.")
    done;
    term
  | Heaviside h_ex -> {num = 1.0; vars = Hashtbl.create 0; h_exs = [h_ex]}
  | Sum _ -> raise (failwith "I cannot pack a Sum expression into a term.")
  | Frac _ -> raise (failwith "I cannot pack a Frac expression into a term.")
;;

let vars_eq (t1 : term) (t2 : term) : bool =
  if Hashtbl.length t1.vars != Hashtbl.length t2.vars then
    false
  else if List.length t1.h_exs != 0 || List.length t2.h_exs != 0 then
    (* Temporary measure to preserve my sanity. *)
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

let rec simplify (ex : expression) : expression =
  match ex with 
  | Num _ -> ex 
  | Var _ -> ex
  | Prod exs -> unpack_term (pack_term (Prod (Array.map simplify exs)))
  | Sum exs ->
    let exs' = Array.map simplify exs in
    let terms = Array.map pack_term exs' in
    let c_terms = Array.of_list (Array.fold_left add_term [] terms) in
    let exs'' = array_filter is_nonzero (Array.map unpack_term c_terms) in
    if Array.length exs'' = 0 then Num 0.0
    else if Array.length exs'' = 1 then exs''.(0)
    else Sum exs''
  | Frac (num, den) -> 
    (match simplify den with
      Num 1.0 -> simplify num
    | Num den_n ->
      (match simplify num with
        Num num_n -> Num (num_n /. den_n)
      | num' -> Frac (num', Num den_n))
    | den' -> Frac (simplify num, den'))
  | Heaviside h_ex -> 
    (match (simplify h_ex) with 
      Num n -> if n > 0.0 then Num 1.0 else Num 0.0
    | h_ex' -> Heaviside h_ex')
;;

let rec expand (ex : expression) : expression =
  let res = (match ex with
    Sum exs ->
      let exs' = Array.map expand exs in
      let (fracs, non_fracs) = split_fracs (Sum exs') in
      if Array.length fracs = 0 then
        Sum (concat_sums non_fracs)
      else
        let t_denom = expand (total_denom fracs None) in
        let mult_frac idx f =
          Prod [| total_denom fracs (Some idx); get_num f |]
        in
        let mult_nonfrac nf =
          Prod [| t_denom; nf |]
        in
        let f_nums = Array.mapi mult_frac fracs in
        let nf_nums = Array.map mult_nonfrac non_fracs in
        
        Frac (expand (Sum (Array.append f_nums nf_nums)), t_denom)
  | Prod exs ->
    let ex' = Array.fold_left (fun x y -> foil (x, y))
      (Prod [||]) (Array.map expand exs) in
    (match ex' with
      Sum _ -> expand ex'
    | Prod _ -> 
      let fracs, non_fracs = split_fracs ex' in
      if Array.length fracs = 0 then ex'
      else 
        let num =  Prod (Array.append (Array.map get_num fracs) non_fracs) in
        Frac (expand num,expand ( total_denom fracs None))
    | _ -> ex')
  | Num _ -> ex
  | Var _ -> ex
  | Frac (num, den) ->
    let num', den' = expand num, expand den in
    (match (num', den') with
      (Frac (nn, nd), Frac (dn, dd)) ->
        Frac (expand (Prod [| nn; dd |]), expand (Prod [| nd; dn |]))
    | (Frac (nn, nd), d) -> Frac (nn, expand (Prod [| nd; d |]))
    | (n, Frac (dn, dd)) -> Frac (expand (Prod [| n; dd |]), dn)
    | (n, d) -> Frac(n, d))
  | Heaviside h_ex -> Heaviside (expand h_ex))
  in
  simplify res
;;

let is_linear (name : string) (t : term) : bool =
  if not (Hashtbl.mem t.vars name) then false
  else if Hashtbl.find t.vars name > 1 then
    raise (failwith "I can't solve non-linear equation.")
  else true
;;

let strip_var (name : string) (t : term) : unit =
  while Hashtbl.mem t.vars name do
    Hashtbl.remove t.vars name
  done
;;

let rec solve_for (name : string) (ex : expression) : expression = 
  match ex with
    Frac (num, den) -> solve_for name num
  | _ ->
    let terms = match ex with
        Sum exs -> Array.map pack_term exs
      | _ -> [| pack_term ex |] in
    let lin_terms = array_filter (is_linear name) terms in
    let non_lin_terms = array_filter (fun t -> not (is_linear name t)) terms in
    Array.iter (strip_var name) lin_terms;
    
    let den_exs = Array.map unpack_term lin_terms in
    let num_exs = Array.map unpack_term non_lin_terms in
    let frac = match (Array.length num_exs, Array.length den_exs) with
        (_, 0) -> raise (failwith "Cannot solve without target variable.")
      | (0, _) -> Num 0.0
      | (1, 1) -> Frac (num_exs.(0), den_exs.(0))
      | (1, _) -> Frac (num_exs.(0), Sum den_exs)
      | (_, 1) -> Frac (Sum num_exs, den_exs.(0))
      | (_, _) -> Frac (Sum num_exs, Sum den_exs)
    in
    Prod [| Num (-1.0); frac |]
;;

let rec subst (sub_name : string) (sub_ex : expression) (ex : expression) = 
  match ex with
    Num num -> ex
  | Var (name, exp) ->
    if name = sub_name then
      if exp = 1 then sub_ex
      else Prod (Array.make exp sub_ex)
    else ex

  | Prod exs ->
    Prod (Array.map (subst sub_name sub_ex) exs)
  | Sum exs -> Sum (Array.map (subst sub_name sub_ex) exs)
  | Frac (num, den) -> Frac (subst sub_name sub_ex num,
                             subst sub_name sub_ex den)
  | Heaviside h_ex -> Heaviside (subst sub_name sub_ex h_ex)
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
  | Heaviside h_ex -> contains_var var_name h_ex
;;

let rec indef_integ (var_name : string) (ex : expression) : expression =
  match ex with
    Num _ -> Prod [| ex; Var (var_name, 1)|]
  | Var (name, exp) -> 
    if name = var_name then
      Prod [| Num (1.0 /. ((float_of_int exp) +. 1.0));
              Var (name, exp + 1)|]
    else Prod [| Var (var_name, 1); ex |]
  | Prod _ ->
    let term = pack_term ex in
    if List.length term.h_exs = 0 then
      if Hashtbl.mem term.vars var_name then
        let exp = Hashtbl.find term.vars var_name in
        Prod [| Frac (Num 1.0, Num (1.0 +.  (float_of_int exp))); ex;
                Var (var_name, 1) |]
      else Prod [| ex; Var (var_name, 1)|]
    else
      let h_ex, h_exs = List.hd term.h_exs, List.tl term.h_exs in
      term.h_exs <- h_exs;
      let ex' = unpack_term term in
      if contains_var var_name h_ex then
        let h_lim = solve_for var_name h_ex in
        let lo, hi = h_lim, Var(var_name, 1)  in
        let h_ex' = simplify (expand (Sum [|hi; Prod [|Num (-1.0); lo|]|])) in
        Prod [| Heaviside h_ex';
                integ lo hi var_name ex' |]
      else
        Prod [| Heaviside h_ex; indef_integ var_name ex'|]
  | Sum exs -> Sum (Array.map (indef_integ var_name) exs)
  | Frac (num, den) -> 
    if contains_var var_name den then
      raise (failwith "I can't evaluate rational integrals.")
    else Frac (indef_integ var_name num, den)
  | Heaviside _ -> indef_integ var_name (Prod [| Num 1.0; ex |])
and integ (lo : expression) (hi : expression)
    (var_name : string) (ex : expression) : expression =
  let indef_ex = indef_integ var_name ex in
  Sum [| subst var_name hi indef_ex; 
         Prod [| Num (-1.0);  subst var_name lo indef_ex|] |]
;;

let main () =
  (* let l3 = Sum [| Num (1.0);
                  Prod [|Num (-1.0); Var("L1", 1) |];
                  Prod [|Num (-1.0); Var("L2", 1) |] |]in *)
  let x = Var("x", 1)
    (*Sum [| Prod [| Var("x1", 1); Var("L1", 1) |];
      Prod [| Var("x2", 1); Var("L2", 1) |];
      Prod [| Var("x3", 1); l3|]|]*) in
  let y = Var("y", 1)
    (*Sum [| Prod [| Var("x1", 1); Var("L1", 1) |];
      Prod [| Var("x2", 1); Var("L2", 1) |];
      Prod [| Var("x3", 1); l3|]|] *) in

  let xl_diff = Sum [| x; Prod [| Num (-1.0); Var("xl", 1) |] |] in
  let xh_diff = Sum [| x; Prod [| Num (-1.0); Var("xh", 1) |] |] in
  let yl_diff = Sum [| y; Prod [| Num (-1.0); Var("yl", 1) |] |] in
  let yh_diff = Sum [| y; Prod [| Num (-1.0); Var("yh", 1) |] |] in

  let x_diff = Sum [| Heaviside(xl_diff);
                      Prod[| Num (-1.0); Heaviside(xh_diff) |] |] in
  let y_diff = Sum [| Heaviside(yl_diff);
                      Prod[| Num (-1.0); Heaviside(yh_diff) |] |] in
  let area = Prod [| x_diff; y_diff |] in

  let area = subst "x1" (Num 0.0) area in
  let area = subst "y1" (Num 0.0) area in
  let area = subst "x2" (Num 1.0) area in
  let area = subst "y2" (Num 0.0) area in
  let area = subst "x3" (Num 0.0) area in
  let area = subst "y3" (Num 1.0) area in

  let area = subst "xl" (Sum [| Var("xl", 1); Prod[| Var("mx", 1); x |] |]) area in
  let area = subst "yl" (Sum [| Var("yl", 1); Prod[| Var("mx", 1); x |] |]) area in
  let area = subst "xh" (Sum [| Var("xh", 1); Prod[| Var("my", 1); y |] |]) area in
  let area = subst "yh" (Sum [| Var("yh", 1); Prod[| Var("my", 1); y |] |]) area in

  let ex = expand area in
  (* print_ex ex;
  printf "# ^ex"; *)
  let lo1 = Num 0.0 in
  let hi1 = Sum [| Num 1.0; Prod [| Num (-1.0); Var("y", 1) |] |] (*Sum [| Num 1.0; Prod [| Num (-1.0); Var("L2", 1) |] |]*) in
  let lo2 = Num 0.0 in
  let hi2 = Num 1.0 in

  let int_ex = expand (integ lo1 hi1 "x" ex) in
  (* print_ex int_ex;
  printf "# ^int ex"; *)
  let int_ex' = integ lo2 hi2 "y" int_ex in
  print_ex int_ex';
  (* printf "# ^int int ex"; *)
;;

main ();;
