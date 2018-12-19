 BEGIN_PROVIDER [double precision, short_range_Hartree_operator, (mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, short_range_Hartree, (N_states)]
 implicit none
 BEGIN_DOC
! short_range_Hartree_operator(i,j) = \int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}
! short_range_Hartree = 0.5 * \sum_{i,j} \rho_{ij} short_range_Hartree_operator(i,j) 
!                     = 0.5 * \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}
 END_DOC
 integer :: i,j,k,l,m,n,istate
 double precision :: get_mo_bielec_integral,get_mo_bielec_integral_erf
 double precision :: integral, integral_erf, contrib
 double precision :: integrals_array(mo_tot_num,mo_tot_num),integrals_erf_array(mo_tot_num,mo_tot_num)
 short_range_Hartree_operator = 0.d0
 short_range_Hartree = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   if(dabs(one_body_dm_average_mo_for_dft(j,i)).le.1.d-12)cycle
   call get_mo_bielec_integrals_i1j1(i,j,mo_tot_num,integrals_array,mo_integrals_map)
   call get_mo_bielec_integrals_erf_i1j1(i,j,mo_tot_num,integrals_erf_array,mo_integrals_erf_map)
   do istate = 1, N_states
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      integral = integrals_array(l,k)
      integral_erf = integrals_erf_array(l,k)
      contrib = one_body_dm_mo_for_dft(i,j,istate) * (integral  - integral_erf)
      short_range_Hartree_operator(l,k,istate) += contrib 
      short_range_Hartree(istate) += contrib * one_body_dm_mo_for_dft(k,l,istate) 
     enddo
    enddo
   enddo
  enddo
 enddo
 short_range_Hartree = short_range_Hartree * 0.5d0
 print*, 'short_range_Hartree',short_range_Hartree
END_PROVIDER


 BEGIN_PROVIDER [double precision, effective_one_e_potential, (mo_tot_num, mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin, (mo_tot_num, mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, shifted_effective_one_e_potential_without_kin, (mo_tot_num, mo_tot_num,N_states)]
 implicit none
 integer :: i,j,istate
 effective_one_e_potential = 0.d0
 BEGIN_DOC 
! effective_one_e_potential(i,j) = <i| v_{H}^{sr} |j> + <i| h_{core} |j> + <i|v_{xc} |j> 
! Taking the expectation value does not provide any energy
! but effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to be used in any WFT calculation
! shifted_effective_one_e_potential_without_kin  =  effective_one_e_potential_without_kin + shifting_constant on the diagonal 
 END_DOC
 do istate = 1, N_states
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    effective_one_e_potential(i,j,istate) = short_range_Hartree_operator(i,j,istate) + mo_nucl_elec_integral(i,j) + mo_kinetic_integral(i,j)    & 
                                   + 0.5d0 * (potential_x_alpha_mo(i,j,istate) + potential_c_alpha_mo(i,j,istate)                               &
                                   +          potential_x_beta_mo(i,j,istate)  + potential_c_beta_mo(i,j,istate)   )
    effective_one_e_potential_without_kin(i,j,istate) = short_range_Hartree_operator(i,j,istate) + mo_nucl_elec_integral(i,j)                   & 
                                   + 0.5d0 * (potential_x_alpha_mo(i,j,istate) + potential_c_alpha_mo(i,j,istate)                               &
                                   +          potential_x_beta_mo(i,j,istate)  + potential_c_beta_mo(i,j,istate)   )
    shifted_effective_one_e_potential_without_kin(j,i,istate) = effective_one_e_potential_without_kin(j,i,istate)
   enddo
  enddo
  do i = 1, mo_tot_num
   shifted_effective_one_e_potential_without_kin(i,i,istate) += shifting_constant(istate)
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, Fock_matrix_expectation_value]
 implicit none
 call get_average(effective_one_e_potential,one_body_dm_average_mo_for_dft,Fock_matrix_expectation_value)

END_PROVIDER 

 BEGIN_PROVIDER [double precision, Trace_v_xc, (N_states)]
&BEGIN_PROVIDER [double precision, Trace_v_H, (N_states)]
&BEGIN_PROVIDER [double precision, Trace_v_Hxc, (N_states)]
 implicit none
 integer :: i,j,istate
 double precision :: dm
 BEGIN_DOC 
! Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
! Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
! Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij} 
 END_DOC
 do istate = 1, N_states
  Trace_v_xc(istate) = 0.d0
  Trace_v_H(istate) = 0.d0
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
     Trace_v_xc(istate) += (potential_x_alpha_mo(j,i,istate) + potential_c_alpha_mo(j,i,istate)) * one_body_dm_mo_alpha_for_dft(j,i,istate) 
     Trace_v_xc(istate) += (potential_x_beta_mo(j,i,istate)  + potential_c_beta_mo(j,i,istate) ) * one_body_dm_mo_beta_for_dft(j,i,istate)
     dm = one_body_dm_mo_alpha_for_dft(j,i,istate) + one_body_dm_mo_beta_for_dft(j,i,istate)
     Trace_v_H(istate) += dm * short_range_Hartree_operator(j,i,istate)
   enddo
  enddo
  Trace_v_Hxc(istate) = Trace_v_xc(istate) + Trace_v_H(istate)
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, DFT_one_e_energy_potential, (mo_tot_num, mo_tot_num,N_states)]
 implicit none
 integer :: i,j,istate
 BEGIN_DOC 
! one_e_energy_potential(i,j) = <i|h_{core}|j> + \int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}
! If one take the expectation value over Psi, one gets the total one body energy
 END_DOC
 do istate = 1, N_states
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    DFT_one_e_energy_potential(j,i,istate) = mo_nucl_elec_integral(j,i) + mo_kinetic_integral(j,i) + short_range_Hartree_operator(j,i,istate) * 0.5d0
   enddo
  enddo
 enddo

END_PROVIDER 

