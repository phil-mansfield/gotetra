module OutTable = AsciiTable.OutTable
module Hist = Histogram.ArrayHistogram
module MC = MCarlo.MakeSwendsenWang(Hist)

let equalize_time = 100

(* Generalize this later. *)
let fixed_histogram_stats dir lat sweep_list =
  let mag_file_name = Printf.sprintf "%s/mag-hist.table" dir in
  let energy_file_name = Printf.sprintf "%s/energy-hist.table" dir in

  let make_name sweeps = Printf.sprintf "%s Sweeps" (string_of_int sweeps) in
  let m_t = OutTable.make ("Magnetism" :: List.map make_name sweep_list) in
  let e_t = OutTable.make ("Energy" :: List.map make_name sweep_list) in
  

  let e_hist, m_hist = MC.create_histograms lat 0 in
  let e_bins, m_bins = Hist.bins e_hist, Hist.bins m_hist in

  let energies = List.mapi (fun i _ -> Hist.bin_value e_hist i) e_bins in
  let mags = List.mapi (fun i _ -> Hist.bin_value m_hist i) m_bins in

  let m_t = OutTable.add_column m_t mags in
  let e_t = OutTable.add_column e_t energies in

  let rec expand_tables sweep_list' e_t m_t =
    match sweep_list' with
      [] -> (e_t, m_t)
    | sweep_head :: sweep_tail ->
      (let e_hist, m_hist = MC.create_histograms lat sweep_head in
       let es = List.map float_of_int (Hist.bins e_hist) in
       let ms = List.map float_of_int (Hist.bins m_hist) in
       expand_tables sweep_tail 
         (OutTable.add_column e_t es)
         (OutTable.add_column m_t ms)
      ) 
  in

  let e_t, m_t = expand_tables sweep_list e_t m_t in
    
  OutTable.write m_t mag_file_name;
  OutTable.write e_t energy_file_name;
;;

let main () =
  if (Array.length Sys.argv) < 4 then
    (Printf.printf "Usage: %s <grid width> <temperature> <dir>\n" Sys.argv.(0);
     exit(1));

  let width = int_of_string Sys.argv.(1) in
  let temp = float_of_string Sys.argv.(2) in
  let lat = MC.init width in
  MC.set_temp lat temp;

  fixed_histogram_stats Sys.argv.(3) lat [1_000; 10_000; 100_000]
;;


main ()
