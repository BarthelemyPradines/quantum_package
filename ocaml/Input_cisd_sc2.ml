open Qptypes;;
open Qputils;;
open Core.Std;;

module Cisd_sc2 : sig
  type t = 
    { n_det_max_cisd_sc2 : Det_number_max.t;
      pt2_max            : PT2_energy.t;
      do_pt2_end         : bool;
    } with sexp
  ;;
  val read : unit -> t
  val to_string : t -> string
  val to_rst : t -> Rst_string.t
end = struct
  type t = 
    { n_det_max_cisd_sc2 : Det_number_max.t;
      pt2_max            : PT2_energy.t;
      do_pt2_end         : bool;
    } with sexp
  ;;

  let get_default = Qpackage.get_ezfio_default "cisd_sc2_selected";;

  let read_n_det_max_cisd_sc2 () = 
    if not (Ezfio.has_cisd_sc2_selected_n_det_max_cisd_sc2 ()) then
       get_default "n_det_max_cisd_sc2"
       |> Int.of_string
       |> Ezfio.set_cisd_sc2_selected_n_det_max_cisd_sc2 
    ;
    Ezfio.get_cisd_sc2_selected_n_det_max_cisd_sc2 ()
    |> Det_number_max.of_int
  ;;


  let read_pt2_max () = 
    if not (Ezfio.has_cisd_sc2_selected_pt2_max ()) then
       get_default "pt2_max"
       |> Float.of_string
       |> Ezfio.set_cisd_sc2_selected_pt2_max
    ;
    Ezfio.get_cisd_sc2_selected_pt2_max ()
    |> PT2_energy.of_float
  ;;


  let read_do_pt2_end () = 
    if not (Ezfio.has_cisd_sc2_selected_do_pt2_end ()) then
       get_default "do_pt2_end"
       |> Bool.of_string
       |> Ezfio.set_cisd_sc2_selected_do_pt2_end
    ;
    Ezfio.get_cisd_sc2_selected_do_pt2_end ()
  ;;


  let read () = 
    { n_det_max_cisd_sc2  = read_n_det_max_cisd_sc2 ();
      pt2_max             = read_pt2_max ();
      do_pt2_end          = read_do_pt2_end ();
    }
  ;;

  let to_string b =
    Printf.sprintf "
n_det_max_cisd_sc2 = %s
pt2_max = %s
do_pt2_end = %s
"
        (Det_number_max.to_string b.n_det_max_cisd_sc2)
        (PT2_energy.to_string b.pt2_max)
        (Bool.to_string b.do_pt2_end)
  ;;

  let to_rst b =
    Printf.sprintf "
Stop when the `n_det` > `n_det_max_cisd_sc2` ::

  n_det_max_cisd_sc2 = %s

Stop when -E(PT2) < `pt2_max` ::

  pt2_max = %s

Compute E(PT2) at the end ::

  do_pt2_end = %s

"
        (Det_number_max.to_string b.n_det_max_cisd_sc2)
        (PT2_energy.to_string b.pt2_max)
        (Bool.to_string b.do_pt2_end)
   |> Rst_string.of_string

end

