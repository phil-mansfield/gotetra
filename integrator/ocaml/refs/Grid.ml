open Array
open Signatures

module FlatGrid : GRID = 
struct
  type coord = int * int

  type 'a grid = { width : int; height : int;
                   x0 : int; y0 : int;
                   x1 : int; y1 : int;
                   (* TODO: make this immutable. *)
                   mutable xs : 'a array }

  type iter_order = SequentialSweep | CheckerboardSweep | RandomSweep

  let int_of_coord g (x, y) = (y - g.y0) * g.width + (x - g.x0)
  let coord_of_int g i = (g.x0 + i mod g.width, g.y0 + i / g.width)      

  let width g = g.width
  let height g = g.height

  let get g (x, y) = 
    g.xs.((x - g.x0) + (y - g.y0) * g.width)
  let set g (x, y) elem = g.xs.((x - g.x0) + (y - g.y0) * g.width) <- elem

  let right g (x, y) = ((if x = g.x1 then g.x0 else x + 1), y)
  let left g (x, y) = ((if x = g.x0 then g.x1 else x - 1), y)
  let up g (x, y) = (x, if y = g.y0 then g.y1 else y - 1)
  let down g (x, y) = (x, if y = g.y1 then g.y0 else y + 1)

  let init (lo_x, lo_y) (hi_x, hi_y) (f : coord -> 'a) =
    let g = { width = hi_x - lo_x + 1; height = hi_y - lo_y + 1;
              x0 = lo_x; y0 = lo_y; x1 = hi_x; y1 = hi_y;
              xs = Array.make 0 (f (lo_x, lo_y))} in
    g.xs <- Array.init (g.width * g.height) (fun i -> f (coord_of_int g i));
    g

  let make lo_coord hi_coord elem = init lo_coord hi_coord (fun _ -> elem)

  let iter order f g =
    let idx_f i = f (coord_of_int g i) in
    match order with
      SequentialSweep -> Array.iteri idx_f g.xs
    | CheckerboardSweep -> raise (failwith "Not Yet Implemented")
    | RandomSweep -> raise (failwith "Not Yet Implemented")

  let fold f x g = Array.fold_left f x g.xs

  (* HATERS GONNA HATE *)
  let foldi (f : 'a -> coord -> 'b -> 'a) (x : 'a) (g : 'b grid) : 'a =
    let sum = ref x in
    iter SequentialSweep (fun c elem -> sum := f !sum c elem) g;
    !sum

  let print g print_f =
    let print_cell i elem = 
      if i mod g.width = 0 then print_string "[ ";
      if (i + 1) mod g.width = 0 then (print_f elem; print_endline " ]") 
      else print_f elem in
    iteri print_cell g.xs;
end
