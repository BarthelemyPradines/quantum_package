BEGIN_PROVIDER [double precision, one_body_dm_mo_alpha_for_dft, (mo_tot_num,mo_tot_num, N_states)]
 implicit none
 BEGIN_DOC
! density used for all DFT calculations based on the density 
 END_DOC
 double precision :: delta_alpha(mo_tot_num,mo_tot_num,N_states)
 if(density_for_dft .EQ. "damping_rs_dft")then
  delta_alpha = one_body_dm_mo_alpha - data_one_body_alpha_dm_mo 
  one_body_dm_mo_alpha_for_dft = data_one_body_alpha_dm_mo + damping_for_rs_dft * delta_alpha
 else if (density_for_dft .EQ. "input_density")then
  one_body_dm_mo_alpha_for_dft = data_one_body_alpha_dm_mo
 else if (density_for_dft .EQ. "WFT")then
  provide mo_coef
  one_body_dm_mo_alpha_for_dft = one_body_dm_mo_alpha
 endif

END_PROVIDER 

BEGIN_PROVIDER [double precision, one_body_dm_mo_beta_for_dft, (mo_tot_num,mo_tot_num, N_states)]
 implicit none
 BEGIN_DOC
! density used for all DFT calculations based on the density 
 END_DOC
 double precision :: delta_beta(mo_tot_num,mo_tot_num,N_states)
 if(density_for_dft .EQ. "damping_rs_dft")then
  delta_beta = one_body_dm_mo_beta - data_one_body_beta_dm_mo 
  one_body_dm_mo_beta_for_dft = data_one_body_beta_dm_mo + damping_for_rs_dft * delta_beta
 else if (density_for_dft .EQ. "input_density")then
  one_body_dm_mo_beta_for_dft = data_one_body_beta_dm_mo
 else if (density_for_dft .EQ. "WFT")then
  provide mo_coef
  one_body_dm_mo_beta_for_dft = one_body_dm_mo_beta
 endif
END_PROVIDER 

BEGIN_PROVIDER [double precision, one_body_dm_mo_for_dft, (mo_tot_num,mo_tot_num, N_states)]
 implicit none
 one_body_dm_mo_for_dft = one_body_dm_mo_beta_for_dft + one_body_dm_mo_alpha_for_dft
END_PROVIDER 

BEGIN_PROVIDER [double precision, one_body_dm_average_mo_for_dft, (mo_tot_num,mo_tot_num)]
 implicit none
 integer :: i
 one_body_dm_average_mo_for_dft = 0.d0
 do i = 1, N_states
  one_body_dm_average_mo_for_dft(:,:) +=  one_body_dm_mo_for_dft(:,:,i) * state_average_weight(i)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, one_body_dm_alpha_ao_for_dft, (ao_num,ao_num,N_states) ]
&BEGIN_PROVIDER [ double precision, one_body_dm_beta_ao_for_dft, (ao_num,ao_num,N_states) ]
 BEGIN_DOC
! one body density matrix on the AO basis based on one_body_dm_mo_alpha_for_dft
 END_DOC
 implicit none
 integer :: i,j,k,l,istate
 double precision :: mo_alpha,mo_beta

 one_body_dm_alpha_ao_for_dft = 0.d0
 one_body_dm_beta_ao_for_dft = 0.d0
 do k = 1, ao_num
  do l = 1, ao_num
   do i = 1, mo_tot_num
    do j = 1, mo_tot_num
     do istate = 1, N_states
      mo_alpha = one_body_dm_mo_alpha_for_dft(j,i,istate)
      mo_beta  = one_body_dm_mo_beta_for_dft(j,i,istate)
      one_body_dm_alpha_ao_for_dft(l,k,istate) += mo_coef(k,i) * mo_coef(l,j) *  mo_alpha
      one_body_dm_beta_ao_for_dft(l,k,istate)  += mo_coef(k,i) * mo_coef(l,j)  *  mo_beta        
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER

