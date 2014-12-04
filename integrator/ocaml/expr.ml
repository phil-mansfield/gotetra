open Printf

type expr = Num of float
            | Var of (string * int)
            | Mult of expr list
            | Sum of expr list
            | Inv of expr
            | Heav of expr
            | Const of expr
;;

type term = Term of float * (string * int) list
            | Compound of expr
       
type var_table = { hs : expr list; 
                   terms : expr list array }
;;

(* This is very slow. *)
let rec to_string (ex : expr) : string = 
  match ex with
    Num n -> sprintf "%g" n
  | Var (v, e) -> if e != 1 then sprintf "%s**%d" v e else v
  | Mult exs -> sprintf "(%s)" (String.concat " * " (List.map to_string exs))
  | Sum exs -> sprintf "(%s)" (String.concat " + " (List.map to_string exs))
  | Inv ex' -> sprintf "1.0 / %s" (to_string ex')
  | Heav ex' -> sprintf "H(%s)" (to_string ex')
  | Const ex' -> sprintf "`%s`" (to_string ex')
;;

let card (ex : expr) : int =
  match ex with
    Const _ -> 0
  | Num _ -> 1
  | Var _ -> 2
  | Mult _ | Sum _ -> 3
  | Inv _ -> 5
  | Heav _ -> 6
;;

let rec list_compare (cmp : 'a -> 'a -> int)
    (xs : 'a list) (ys : 'a list) : int = 
  match (xs, ys) with
    ([], []) -> 0
  | ([], _) -> -1
  | (_, []) -> +1
  | (x::xs', y::ys') ->
    let compare_xy = cmp x y in 
    if compare_xy != 0 then compare_xy else list_compare cmp xs' ys'
;;

let rec compare (ex1 : expr) (ex2: expr) : int = 
  match (ex1, ex2) with
    (Num n1, Num n2) -> int_of_float (n1 -. n2)
  | (Var (v1, e1), Var (v2, e2)) -> 
    let v_comp = String.compare v1 v2 in
    if v_comp != 0 then v_comp else e1 - e2
  | (Mult exs1, Mult exs2) | (Sum exs1, Sum exs2) ->
    let sexs1, sexs2 = List.sort compare exs1, List.sort compare exs2 in
    list_compare compare sexs1 sexs2
  | (Inv (ex1), Inv (ex2)) -> compare ex1 ex2
  | (Heav ex1', Heav ex2') | (Const ex1', Const ex2') -> compare ex1' ex2'
  | (ex1', ex2') -> (card ex1') - (card ex2')  
;;

let equal (ex1 : expr) (ex2 : expr) : bool =
  0 = (compare ex1 ex2)
;;

let is_sum ex = match ex with Sum _ -> true | _ -> false ;;
let is_mult ex = match ex with Mult _ -> true | _ -> false ;;
let is_const ex = match ex with Const _ -> true | _ -> false ;;
let is_leaf ex = match ex with Var _ | Num _ -> true | _ -> false ;;
let is_one ex = match ex with Num n when n = 1.0 -> true | _ -> false ;;

let rec check (ex : expr) : bool =
  match ex with
    Num _ -> true
  | Var (v, e) -> v != "" && e != 0
  | Mult exs -> List.length exs > 1 && not (List.exists is_mult exs)
  | Sum exs -> List.length exs > 1 && not (List.exists is_sum exs)
  | Inv ex' -> check ex'
  | Heav ex' | Const ex' -> check ex'
;;

let is_constant (v : string) (ex : expr) : bool =
  match ex with 
    Num _ -> true
  | Var (v', _) -> v' != v
  | Mult exs | Sum exs -> List.for_all is_const exs
  | Inv ex' -> is_const ex'
  | Heav ex' -> is_const ex'
  | Const _ -> true
;;

let rec flag_constants (v : string) (ex : expr) : expr = 
  let ex' = match ex with 
      Sum exs -> Sum (List.map (flag_constants v) exs)
    | Mult exs -> Mult (List.map (flag_constants v) exs)
    | Heav ex'' -> Heav (flag_constants v ex'')
    | Const _ -> failwith "Flagging an already constant expr."
    | _ -> ex in
  if is_constant v ex' then Const ex' else ex'
;;

let rec unflag_constants (ex : expr) : expr =
  match ex with
    Num _ | Var _-> ex
  | Mult exs -> Mult (List.map unflag_constants exs)
  | Sum exs -> Sum (List.map unflag_constants exs)
  | Inv ex' -> Inv (unflag_constants ex')
  | Heav ex' -> Heav (unflag_constants ex')
  | Const ex' -> unflag_constants ex'
;;

let to_term (ex : expr) : term = 
  match ex with 
    Num n -> Term (n, [])
  | Var (v, e) -> Term (1.0, [(v, e)])
  | Mult exs when not (List.for_all is_leaf exs) -> Compound ex
  | Mult exs -> 
    let fold_term (n, vars) ex =
      match ex with
        Num n' -> (n *. n', vars)
      | Var (v, e) when List.mem_assoc v vars -> 
        (n, (v, e + List.assoc v vars) :: (List.remove_assoc v vars))
      | Var (v, e) -> (n, (v, e) :: vars)
      | _ -> failwith "Impossible"
    in
    let n, vars = List.fold_left fold_term (1.0, []) (List.sort compare exs) in
    Term (n, vars)
  | Inv _ -> Compound ex
  | Sum _ | Heav _  -> Compound ex
  | Const _ -> failwith "Tried to convert Const to term."
;;

let prod_list_to_expr (n : float) (exs : expr list) : expr =
  match (n, exs) with
    (0.0, _) -> Num 0.0
  | (_, []) -> Num n
  | (1.0, [ex]) -> ex
  | (1.0, exs) -> Mult exs
  | (_, exs) -> Mult ((Num n) :: exs)
;;

(* Very slow. *)
let can_add_term (t1 : term) (t2 : term) : bool =
  match (t1, t2) with
    (_, Compound _) | (Compound _, _) -> false (* One day... I promise. *)
  | (Term (_, vars1), Term (_, vars2)) ->
    if List.length vars1 != List.length vars2 then false
    else List.for_all (fun var -> List.mem var vars2) vars1
;;

let add_term (t1 : term) (t2 : term) : term =
  match (t1, t2) with
    Term(n1, var1), Term(n2, var2) -> Term(n1 +. n2, var1)
  | _, _ -> failwith "Tried to add terms that can't be added."
;;

let append_term ((ts, cs) : (term list * term list))
    (t : term) : (term list * term list) =
  match t with
    Term _ ->
      let rec append_term' ts' =
        match ts' with
          [] -> [t]
        | t'' :: ts'' ->
          if can_add_term t'' t then (add_term t t''):: ts''
          else t'' :: append_term' ts''
      in (append_term' ts, cs)
  | Compound _ -> (ts, t :: cs)
;;

let of_term (t : term) : expr =
  match t with 
    Compound ex -> ex
  | Term (n, vars) ->
    prod_list_to_expr n
      (List.fold_left 
         (fun vs (v, e) -> if e = 0 then [Num 1.0] else Var(v, e) :: vs)
         [] vars)
;;

let expand (ex : expr) : expr list =
  match ex with
    Sum exs | Mult exs -> exs
  | _ -> [ex]
;;

let is_inv (ex : expr) : bool =
  match ex with
    Num _ | Var _ | Inv _ -> true
  | _ -> false
;;

let invert (ex : expr) : expr = 
  match ex with
    Num n -> Num (1.0 /. n)
  | Var (v, e) -> Var (v, -e)
  | Inv ex' -> ex'
  | _ -> failwith "Tried to invert non-invertible expr."
;;

let rec simplify (ex : expr) : expr =
  match ex with
    Num _ -> ex
  | Var (_, 0) -> Num 1.0
  | Var _ -> ex
  | Mult exs ->
    let exs' = List.map simplify exs in
    let mults, n_mults = List.partition is_mult exs' in
    let exs'' : expr list = List.concat (n_mults::(List.map expand mults)) in
    let ls, nls = List.partition is_leaf exs'' in
    let m_exs = match of_term (to_term (Mult ls)) with
        Num 0.0 -> [Num 0.0]
      | (Num _ | Var _) as ex' -> ex' :: nls
      | Mult m_exs' -> List.append m_exs' nls
      | _ -> failwith "Impossible of_term value." in
    (match List.filter (fun ex' -> not (is_one ex')) m_exs with
      [] -> Num 1.0
    | [ex'] -> ex'
    | m_exs' -> Mult m_exs')
   
  | Sum exs ->
    let exs' = List.map simplify exs in
    let sums, n_sums = List.partition is_sum exs' in
    let exs'' = List.concat (n_sums :: (List.map expand sums)) in
    let ts = List.map to_term exs'' in
    let (ts', cs') = List.fold_left append_term ([], []) ts in
    let exs''' = List.map of_term (List.append ts' cs') in
    (match exs''' with [ex] -> ex | _ -> Sum exs''')
  | Inv ex' ->
    (match simplify ex' with 
      ex'' when is_inv ex'' -> invert ex''
    | Mult exs ->
      let is, nis = List.partition is_inv exs in
      (match nis with
        [] -> Mult (List.map invert is)
      | [ex] -> Mult (Inv ex :: (List.map invert is))
      | _ -> Mult (Inv (Mult (nis)) :: (List.map invert is)))
    | Const _ -> failwith "Tried to simplify Const."
    | ex''  -> Inv ex'')
  | Heav ex' ->
    (match simplify ex' with
      Num n when n > 0.0 -> Num 1.0
    | Num n -> Num 0.0
    | ex'' -> Heav ex'')
  | Const _ -> failwith "Tried to simplify Const."
;;

let subst (v : string) (s_ex : expr) (ex : expr) : expr =
  match ex with 
    Num _ -> ex
  | Var (v', e') -> Array.to_list (Array.make e' s_ex)
  | Sum exs -> Sum (List.map (subst v s_ex) exs)
  | Mult exs -> Mult (List.map (subst v s_ex) exs)
  | Inv ex' -> Inv (subst v s_ex ex')
  | Heav ex' -> Heav (subst v s_ex ex')
  | Const _ -> failwith "Tried to subst into a Const expr."
;;

let of_var_table (vt : var_table) : expr =
  failwith "NYI"
;;

let to_var_table (ex : expr) : var_table = 
  failwith "NYI"
;;

let main () =
  let ex = in
  print_endline (to_string ex);
  print_endline (to_string (simplify ex));
;;

main ();;
