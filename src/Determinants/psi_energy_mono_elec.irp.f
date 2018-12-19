 BEGIN_PROVIDER [ double precision, psi_energy_h_core, (N_states) ]
  implicit none
  integer :: i
  integer :: j,k
  double precision :: tmp(mo_tot_num,mo_tot_num),mono_ints(mo_tot_num,mo_tot_num)
  BEGIN_DOC
! psi_energy_h_core                = <Psi| h_{core} |Psi>
! computed using the one_body_dm_mo_alpha+one_body_dm_mo_beta and mo_mono_elec_integral
  END_DOC
  psi_energy_h_core = 0.d0
  do i = 1, N_states
   do j = 1, mo_tot_num
    do k = 1, mo_tot_num
     psi_energy_h_core(i) += mo_mono_elec_integral(k,j) * (one_body_dm_mo_alpha(k,j,i) + one_body_dm_mo_beta(k,j,i))
    enddo
   enddo
  enddo
 double precision :: accu
 do i = 1, N_states
  do j = 1, mo_tot_num
   accu += one_body_dm_mo_alpha(j,j,i) + one_body_dm_mo_beta(j,j,i) 
  enddo
  accu = (elec_alpha_num + elec_beta_num ) / accu
  psi_energy_h_core(i) = psi_energy_h_core(i) * accu
 enddo
END_PROVIDER

