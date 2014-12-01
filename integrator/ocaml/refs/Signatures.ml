(** A repository for all the signatures used by the bellagio package. *)
open Array

module type GRID = 
sig
  (** The type of a grid containing elements of type ['a]. *)
  type 'a grid

  (** The type of an index into the grid.*)
  type coord = int * int

  (** The type representing the order in which iteration is done over the 
      grid.*)
  type iter_order = SequentialSweep | CheckerboardSweep | RandomSweep

  (** [int_of_coord g (x, y)] converts a coordinate [(x, y)] into a integer 
      between 0 and [(Grid.width g) * (Grid.height g) - 1], inclusive. *)
  val int_of_coord : 'a grid -> coord -> int

  (** [coord_of_int g i] converts an index between 0 and
      [(Grid.width g) * (Grid.height g) - 1], inclusive, to a unique coordinate.
      indexing into [g]. *)
  val coord_of_int : 'a grid -> int -> coord

  (** [width c] returns the number of elements in a single row of [g]. *)
  val width : 'a grid -> int

  (** [height g] returns the number of elements in a single column of [g]. *)
  val height : 'a grid -> int
    
  (** [get g (x, y)] returns the element stored at the coordinate [(x, y)] in 
      the grid [g].*)
  val get : 'a grid -> coord -> 'a

  (** [set g (x, y) elem] sets the element stored a the coordinate [(x, y)] in
      [g] to [elem].*)
  val set : 'a grid -> coord -> 'a -> unit

  (** [make (lo_x, lo_y) (hi_x, hi_y) elem] creates a grid with width 
      [hi_x - lo_x], height [hi_y - lo_y], and a bottom-left coordinate
      (lo_x, hi_x). Every element is initialized to [elem]. *)
  val make : coord -> coord -> 'a -> 'a grid

  (** [iter (lo_x, lo_y) (hi_x, hi_y) f] creates a grid with width 
      [hi_x - lo_x], height [hi_y - lo_y], and a bottom-left coordinate
      (lo_x, hi_x). Each element is initialized to [f (x, y)], where [(x, y)],
      is that element's coordinate. *)
  val init : coord -> coord -> (coord -> 'a) -> 'a grid
    
  (** [iter order f g] applies the funciton [f] to all the coordinates and
      elements of [g] in the order specified by [order]. *)
  val iter : iter_order -> (coord -> 'a -> unit) -> 'a grid -> unit

  (** [fold f x g] applies some function [f] to every element in [g] along with
      the result from the previous application of [f]. This is equivelent to the
      way [fold] is structured on literally every other data structure that's
      ever existed. *)
  val fold : ('a -> 'b -> 'a) -> 'a -> 'b grid -> 'a

  (** [foldi f x g] is equivelent to [Grid.fold], except the coordinate of the
      element is also passed to [f].*)
  val foldi : ('a -> coord -> 'b -> 'a) -> 'a  -> 'b grid -> 'a

  (** [print g print_f] applies [print_f] to every element of [g] in order and
      seperates rows by newlines. *)
  val print : 'a grid -> ('a -> unit) -> unit

  (** [right g (x, y)]*)
  val right : 'a grid -> coord -> coord
  val left : 'a grid -> coord -> coord
  val up : 'a grid -> coord -> coord
  val down : 'a grid -> coord -> coord
end

module type GROUPER =
sig
  type group_record
  type node = int
  type edge = node * node
  type group_id = int
    
  val union : node -> edge list -> group_record
    
  val find : group_record -> node -> group_id
  val group_ids : group_record -> group_id list
  
  val group : group_record -> group_id -> node list
  val group_size : group_record -> group_id -> int
    
  val largest_group : group_record -> group_id
end

module type HISTOGRAM =
sig
  type histogram
  
  val init : float -> float -> int -> histogram
  val init_bounded : float -> float -> int -> histogram

  val bin_value : histogram -> int -> float
  val bin_index : histogram -> float -> int

  val add : histogram -> float -> unit
  val bins : histogram -> int list
end

module type MCARLO = 
sig
  type lattice
  type histogram
  type bond_type = NN1 | NN2 | NN3 | DMD | SQR | MAG
  type normalize_type = Random2By2 | InPlace2By2 | InPlace3By3

  val init : int -> lattice
  val site_count : lattice -> int

  val reset : lattice -> unit
  val sweep : lattice -> unit
  val set_temp : lattice -> float -> unit
  val set_bond_types : lattice -> (bond_type * float) array -> unit

  val renormalize : lattice -> normalize_type -> float option -> lattice

  val correlation : lattice -> bond_type -> float

  val energy : lattice -> float
  val magnetization : lattice -> float
  val max_cluster_size : lattice -> int
    
  val create_histograms : lattice -> int -> (histogram * histogram)

  val energy_from_hist : histogram -> float -> float -> float
  val c_from_hist : histogram -> float -> float -> float

  val print : lattice -> unit
end

module type OUT_TABLE =
sig
  type out_table

  val make : string list -> out_table

  val add_row : out_table -> float list -> out_table
  val add_column : out_table -> float list -> out_table

  val write : out_table -> string -> unit
end

module type MCARLO_FUNCTOR = functor (Hist : HISTOGRAM) ->
    MCARLO with type histogram = Hist.histogram
