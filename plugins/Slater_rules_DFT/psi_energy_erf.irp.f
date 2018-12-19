
BEGIN_PROVIDER [ double precision, psi_energy_erf, (N_states) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes e_0 = <Psi|W_{ee}^{lr}|Psi>/<Psi|Psi>
  !
  END_DOC
  integer :: i
  call u_0_H_u_0_erf(psi_energy_erf,psi_coef,N_det,psi_det,N_int,N_states,psi_det_size)
  do i=N_det+1,N_states
    psi_energy_erf(i) = 0.d0
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_energy_core_and_sr_hartree, (N_states) ]
  implicit none
  BEGIN_DOC
! psi_energy_core                = <Psi| h_{core} + v_{H}^{sr}|Psi>
  END_DOC
  psi_energy_core_and_sr_hartree = psi_energy_core + short_range_Hartree
END_PROVIDER


 BEGIN_PROVIDER [ double precision, psi_energy_core, (N_states) ]
  implicit none
  integer :: i
  integer :: j,k
  double precision :: tmp(mo_tot_num,mo_tot_num),mono_ints(mo_tot_num,mo_tot_num)
  BEGIN_DOC
! psi_energy_core                = <Psi| h_{core} |Psi>
  END_DOC
  psi_energy_core = 0.d0
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    mono_ints(k,j) = mo_nucl_elec_integral(k,j) + mo_kinetic_integral(k,j)
   enddo
  enddo
  do i = 1, N_states
   do j = 1, mo_tot_num
    do k = 1, mo_tot_num
     tmp(k,j) = one_body_dm_mo_alpha_for_dft(k,j,i) + one_body_dm_mo_beta_for_dft(k,j,i)
    enddo
   enddo
   call get_average(mono_ints,tmp,psi_energy_core(i))
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, total_range_separated_electronic_energy, (N_states) ]
  implicit none
  BEGIN_DOC
! Total_range_separated_electronic_energy = <Psi| h_{core} |Psi> + (1/2) <Psi| v_{H}^{sr} |Psi> + <i|W_{ee}^{lr}|i> + E_{x} + E_{c}
  END_DOC
   total_range_separated_electronic_energy = psi_energy_core + short_range_Hartree + psi_energy_erf + energy_x + energy_c 
END_PROVIDER


BEGIN_PROVIDER [ double precision, two_elec_energy_dft, (N_states) ]
  implicit none
  BEGIN_DOC
! two_elec_energy_dft = (1/2) <Psi| v_{H}^{sr} |Psi> + <i|W_{ee}^{lr}|i> 
  END_DOC
   two_elec_energy_dft = short_range_Hartree + psi_energy_erf
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ref_bitmask_energy_erf ]
&BEGIN_PROVIDER [ double precision, bi_elec_ref_bitmask_energy_erf ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Energy of the reference bitmask used in Slater rules
  END_DOC
  
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: i,j
  
  call bitstring_to_list(ref_bitmask(1,1), occ(1,1), i, N_int)
  call bitstring_to_list(ref_bitmask(1,2), occ(1,2), i, N_int)
  
  
  ref_bitmask_energy_erf = 0.d0
  bi_elec_ref_bitmask_energy_erf = 0.d0
  
  do j= 1, elec_alpha_num
    do i = j+1, elec_alpha_num
      bi_elec_ref_bitmask_energy_erf += mo_bielec_integral_erf_jj_anti(occ(i,1),occ(j,1))
      ref_bitmask_energy_erf += mo_bielec_integral_erf_jj_anti(occ(i,1),occ(j,1))
    enddo
  enddo
  
  do j= 1, elec_beta_num
    do i = j+1, elec_beta_num
      bi_elec_ref_bitmask_energy += mo_bielec_integral_erf_jj_anti(occ(i,2),occ(j,2))
      ref_bitmask_energy += mo_bielec_integral_erf_jj_anti(occ(i,2),occ(j,2))
    enddo
    do i= 1, elec_alpha_num
      bi_elec_ref_bitmask_energy += mo_bielec_integral_erf_jj(occ(i,1),occ(j,2))
      ref_bitmask_energy += mo_bielec_integral_erf_jj(occ(i,1),occ(j,2))
    enddo
  enddo
  
END_PROVIDER

