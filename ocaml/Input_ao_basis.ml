open Qptypes;;
open Qputils;;
open Core.Std;;

module Ao_basis : sig
  type t = 
    { ao_basis        : string ;
      ao_num          : AO_number.t ;
      ao_prim_num     : AO_prim_number.t array;
      ao_prim_num_max : AO_prim_number.t;
      ao_nucl         : Nucl_number.t array;
      ao_power        : Symmetry.Xyz.t array;
      ao_coef         : AO_coef.t array;
      ao_expo         : AO_expo.t array;
    } with sexp
  ;;
  val read : unit -> t
  val to_string : t -> string
  val to_md5 : t -> MD5.t
  val to_rst : t -> Rst_string.t
end = struct
  type t = 
    { ao_basis        : string ;
      ao_num          : AO_number.t ;
      ao_prim_num     : AO_prim_number.t array;
      ao_prim_num_max : AO_prim_number.t;
      ao_nucl         : Nucl_number.t array;
      ao_power        : Symmetry.Xyz.t array;
      ao_coef         : AO_coef.t array;
      ao_expo         : AO_expo.t array;
    } with sexp
  ;;

  let get_default = Qpackage.get_ezfio_default "ao_basis";;

  let read_ao_basis () = 
    if not (Ezfio.has_ao_basis_ao_basis ()) then
       Ezfio.set_ao_basis_ao_basis ""
    ;
    Ezfio.get_ao_basis_ao_basis ()
  ;;

  let read_ao_num () =
    Ezfio.get_ao_basis_ao_num ()
    |> AO_number.of_int
  ;;

  let read_ao_prim_num () =
    (Ezfio.get_ao_basis_ao_prim_num () ).Ezfio.data
    |> Ezfio.flattened_ezfio_data
    |> Array.map ~f:AO_prim_number.of_int
  ;;

  let read_ao_prim_num_max () =
    (Ezfio.get_ao_basis_ao_prim_num () ).Ezfio.data
    |> Ezfio.flattened_ezfio_data
    |> Array.fold ~f:(fun x y -> if x>y then x else y) ~init:0
    |> AO_prim_number.of_int
  ;;

  let read_ao_nucl () =
    (Ezfio.get_ao_basis_ao_nucl () ).Ezfio.data
    |> Ezfio.flattened_ezfio_data
    |> Array.map ~f:Nucl_number.of_int
  ;;

  let read_ao_power () =
    let x = Ezfio.get_ao_basis_ao_power () in
    let dim = x.Ezfio.dim.(0) in
    let data = Ezfio.flattened_ezfio_data x.Ezfio.data in
    let result = Array.init dim ~f:(fun x -> "") in
    for i=1 to dim
    do
      if (data.(i-1) > 0) then
        result.(i-1) <- result.(i-1)^"x"^(Int.to_string data.(i-1));
      if (data.(dim+i-1) > 0) then
        result.(i-1) <- result.(i-1)^"y"^(Int.to_string data.(dim+i-1));
      if (data.(2*dim+i-1) > 0) then
        result.(i-1) <- result.(i-1)^"z"^(Int.to_string data.(2*dim+i-1));
    done;
    Array.map ~f:Symmetry.Xyz.of_string result
  ;;

  let read_ao_coef () =
    (Ezfio.get_ao_basis_ao_coef () ).Ezfio.data
    |> Ezfio.flattened_ezfio_data
    |> Array.map ~f:AO_coef.of_float
  ;;

  let read_ao_expo () =
    (Ezfio.get_ao_basis_ao_expo () ).Ezfio.data
    |> Ezfio.flattened_ezfio_data
    |> Array.map ~f:AO_expo.of_float
  ;;

  let read () =
    { ao_basis        = read_ao_basis ();
      ao_num          = read_ao_num () ;
      ao_prim_num     = read_ao_prim_num ();
      ao_prim_num_max = read_ao_prim_num_max ();
      ao_nucl         = read_ao_nucl ();
      ao_power        = read_ao_power ();
      ao_coef         = read_ao_coef () ;
      ao_expo         = read_ao_expo () ;
    }
  ;;
    
  let to_long_basis b =
    let ao_num = AO_number.to_int b.ao_num in
    let gto_array = Array.init (AO_number.to_int b.ao_num)
      ~f:(fun i ->
        let s = Symmetry.Xyz.to_symmetry b.ao_power.(i) in
        let ao_prim_num = AO_prim_number.to_int b.ao_prim_num.(i) in
        let prims = List.init ao_prim_num ~f:(fun j ->
          let prim = { Primitive.sym  = s ;
                       Primitive.expo = b.ao_expo.(ao_num*j+i)
                     }
          in
          let coef = b.ao_coef.(ao_num*j+i) in
          (prim,coef)
        ) in
        Gto.of_prim_coef_list prims
      )
    in
    let rec do_work accu sym gto nucl =
      match (sym, gto, nucl) with
        | (s::srest, g::grest, n::nrest) -> 
          do_work ((s,g,n)::accu) srest grest nrest
        | ([],[],[]) -> List.rev accu
        | _ -> assert false
    in
    do_work [] 
      (Array.to_list b.ao_power)
      (Array.to_list gto_array)
      (Array.to_list b.ao_nucl)
  ;;

  let to_basis b =
    to_long_basis b
    |> Long_basis.to_basis
  ;;
      
  let to_rst b =
    let short_basis = to_basis b in
    Printf.sprintf "
Name of the AO basis ::

  ao_basis = %s

Basis set ::
  
%s

" b.ao_basis
    (Basis.to_string short_basis 
     |> String.split ~on:'\n'
     |> List.map ~f:(fun x-> "  "^x)
     |> String.concat ~sep:"\n"
    )
  |> Rst_string.of_string
  ;;

  let to_md5 b =
    let short_basis = to_basis b in
    Basis.to_md5 short_basis
  ;;
    
  let to_string b =
    Printf.sprintf "
ao_basis        = %s
ao_num          = %s
ao_prim_num     = %s
ao_prim_num_max = %s
ao_nucl         = %s
ao_power        = %s
ao_coef         = %s
ao_expo         = %s
md5             = %s
"
    b.ao_basis
    (AO_number.to_string b.ao_num)
    (b.ao_prim_num |> Array.to_list |> List.map
      ~f:(AO_prim_number.to_string) |> String.concat ~sep:", " )
    (AO_prim_number.to_string b.ao_prim_num_max)
    (b.ao_nucl |> Array.to_list |> List.map ~f:Nucl_number.to_string |>
      String.concat ~sep:", ")
    (b.ao_power |> Array.to_list |> List.map ~f:(fun x->
      "("^(Symmetry.Xyz.to_string x)^")" )|> String.concat ~sep:", ")
    (b.ao_coef  |> Array.to_list |> List.map ~f:AO_coef.to_string
      |> String.concat ~sep:", ")
    (b.ao_expo  |> Array.to_list |> List.map ~f:AO_expo.to_string
      |> String.concat ~sep:", ")
    (to_md5 b |> MD5.to_string )

  ;;
end
