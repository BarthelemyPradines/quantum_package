BEGIN_PROVIDER [double precision, electronic_energy_mr_dft, (N_states)]
 implicit none
 BEGIN_DOC
 ! Energy for the multi determinantal DFT calculation
 END_DOC
 
 if(projected_wft_for_dft)then
  print*,'You are using a projected WFT which uses an eigenvalue stored in the EZFIO folder' 
  electronic_energy_mr_dft = data_energy_proj + short_range_Hartree + energy_x + energy_c - Trace_v_Hxc
 else
  print*,'You are using a variational method which uses the wave function stored in the EZFIO folder'
  electronic_energy_mr_dft = total_range_separated_electronic_energy
 endif


END_PROVIDER 

 subroutine print_variational_energy_dft
 implicit none
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
  write(*, '(A22,X,F32.10)') 'mu_erf              = ',mu_erf          
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',electronic_energy_mr_dft+nuclear_repulsion
  print*, ''
  print*, 'Component of the energy ....'
  print*, ''
  write(*, '(A22,X,F16.10)') 'nuclear_repulsion   = ',nuclear_repulsion
  write(*, '(A22,X,F16.10)') 'psi_energy_erf      = ',psi_energy_erf      
  write(*, '(A22,X,F16.10)') 'psi_energy_core     = ',psi_energy_core    
  write(*, '(A22,X,F16.10)') 'short_range_Hartree = ',short_range_Hartree
  write(*, '(A22,X,F16.10)') 'two_elec_energy     = ',two_elec_energy_dft
  write(*, '(A22,X,F16.10)') 'energy_x            = ',energy_x         
  write(*, '(A22,X,F16.10)') 'energy_c            = ',energy_c          
  print*, ''
  print*,  '****************************************'
  print*, ''
  write(*, '(A22,X,F16.10)') 'Approx eigenvalue   = ',Fock_matrix_expectation_value + psi_energy_erf
  write(*, '(A22,X,F16.10)') 'Trace_v_xc          = ',Trace_v_xc
 
  write(*, '(A28,X,F16.10)') 'Variational energy of Psi = ',psi_energy
  write(*, '(A28,X,F16.10)') 'psi_energy_bielec         = ',psi_energy_bielec
  write(*, '(A28,X,F16.10)') 'psi_energy_h_core       = ',psi_dft_energy_kinetic+psi_dft_energy_nuclear_elec
! write(*, '(A28,X,F16.10)') 'corrected Multi-det correl= ',Energy_c_md_on_top(1)
 end

subroutine print_variational_energy_dftTest_julien
implicit none
  print*,  '****************************************'
  print*,  ' Range separated DFT energy new Toulouse method '
  write(*, '(A22,X,F32.10)') 'mu_erf              = ',mu_erf
  print*, ''
  print*, 'Component of the energy ....'
  print*, ''
  write(*, '(A22,X,F16.10)') 'nuclear_repulsion     = ',nuclear_repulsion
  write(*, '(A22,X,F16.10)') 'psi_energy_erf        = ',psi_energy_erf
  write(*, '(A22,X,F16.10)') 'psi_energy_core       = ',psi_energy_core
  write(*, '(A22,X,F16.10)') 'E^LR                  = ',psi_energy_core + psi_energy_erf
  !write(*, '(A22,X,F16.10)') 'E_Hxc_Toulouse_1    = ',energy_Hxc(1)
  !write(*, '(A22,X,F16.10)') 'E_Hxc_Toulouse_2    = ',energy_Hxc_bis(1)
  !write(*, '(A22,X,F16.10)') 'E_Hxc_Toulouse_3    = ',energy_Hxc_ter(1)
  !write(*, '(A22,X,F16.10)') 'E_Hxc_Toulouse_4    = ',energy_Hxc_4(1)
  !write(*, '(A22,X,F16.10)') 'E_Hxc_Toulouse_5    = ',energy_Hxc_5(1)
  !write(*, '(A22,X,F16.10)') 'Tot_Toulouse_1      = ',nuclear_repulsion +  psi_energy_core + psi_energy_erf + energy_Hxc
  !write(*, '(A22,X,F16.10)') 'Tot_Toulouse_2      = ',nuclear_repulsion +  psi_energy_core + psi_energy_erf + energy_Hxc_bis
  !write(*, '(A22,X,F16.10)') 'Tot_Toulouse_3      = ',nuclear_repulsion + psi_energy_core + psi_energy_erf + energy_Hxc_ter
  !write(*, '(A22,X,F16.10)') 'Tot_Toulouse_4      = ',nuclear_repulsion + psi_energy_core + psi_energy_erf + energy_Hxc_4
  !write(*, '(A22,X,F16.10)') 'Tot_Toulouse_5      = ',nuclear_repulsion +psi_energy_core + psi_energy_erf + energy_Hxc_5
  write(*, '(A22,X,F16.10)') 'Expect value V_Sr,    = ',exp_value_V_SR_mu(1)
  write(*, '(A22,X,F16.10)') 'E_x,sr+E_c,sr+E_H,sr  = ',short_range_Hartree+energy_x+energy_c
end

subroutine print_variational_energy_dft_mu_of_r
 implicit none

 print*,  '****************************************'
 print*, 'Functional used   = ',md_correlation_functional
 print*,  '****************************************'
 print*, 'mu_of_r_potential = ',mu_of_r_potential
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 if(md_correlation_functional.EQ."basis_set_LDA")then
   write(*, '(A34,X,F16.10)') 'TOTAL ENERGY ECMD LDA     =',psi_energy+Energy_c_md_mu_of_r_LDA+nuclear_repulsion
   write(*, '(A28,X,F16.10)') 'TOTAL ENERGY EC   LDA     =',psi_energy+energy_c_LDA_mu_of_r+nuclear_repulsion
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print*, ''
   write(*, '(A28,X,F16.10)') 'Variational energy of Psi = ',psi_energy + nuclear_repulsion
   print*, 'Component of the energy ....'
   write(*, '(A28,X,F16.10)') 'psi_energy_bielec         = ',psi_energy_bielec
   write(*, '(A28,X,F16.10)') 'psi_energy_h_core       = ',psi_dft_energy_kinetic+psi_dft_energy_nuclear_elec
   write(*, '(A28,X,F16.10)') 'psi_dft_energy_kinetic        = ',psi_dft_energy_kinetic
   write(*, '(A28,X,F16.10)') 'psi_dft_energy_nuclear_elec   = ',psi_dft_energy_nuclear_elec
   write(*, '(A28,X,F16.10)') 'nuclear_repulsion         = ',nuclear_repulsion
   print*, ''
   write(*, '(A28,X,F16.10)') 'DFT mu(r)     correlation = ',Energy_c_md_mu_of_r_LDA
   write(*, '(A28,X,F16.10)') 'DFT mu(r) NO MD corr      = ',energy_c_LDA_mu_of_r
   print*, ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else if(md_correlation_functional.EQ."basis_set_on_top_PBE")then
   write(*, '(A34,X,F16.10)') 'TOTAL ENERGY ECMD PBE n2 UEG      =',psi_energy+Energy_c_md_on_top_PBE_mu_of_r_UEG+nuclear_repulsion
   write(*, '(A34,X,F16.10)') 'TOTAL ENERGY ECMD PBE n2 NO UEG   =',psi_energy+Energy_c_md_on_top_PBE_mu_of_r+nuclear_repulsion
   write(*, '(A34,X,F16.10)') 'TOTAL ENERGY ECMD LDA             =',psi_energy+Energy_c_md_mu_of_r_LDA+nuclear_repulsion
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print*, ''
   write(*, '(A28,X,F16.10)') 'Variational energy of Psi   = ',psi_energy
   print*, 'Component of the energy ....'
   print*, ''
   write(*, '(A28,X,F16.10)') 'nuclear_repulsion           = ',nuclear_repulsion
   write(*, '(A28,X,F16.10)') 'Energy ECMD UEG        = ',Energy_c_md_on_top_PBE_mu_of_r_UEG
   write(*, '(A28,X,F16.10)') 'Energy ECMD NO UEG     = ',Energy_c_md_on_top_PBE_mu_of_r
   print*,''
   write(*, '(A28,X,F16.10)') 'Energy ECMD LDA        = ',Energy_c_md_mu_of_r_LDA
 endif 
  if(.true.)then
   write(*, '(A28,X,F16.10)') 'mu_average for basis set  = ',mu_average
  endif
 if(mu_of_r_potential.EQ."hf_integral")then
  print*,'integral_f_hf              = ',integral_f_hf
  print*,'HF_alpha_beta_bielec_energy= ',HF_alpha_beta_bielec_energy
  print*,'mu(r) expectation values .......'
  print*,'HF_mu_of_r_bielec_energy   = ',HF_mu_of_r_bielec_energy
  print*,'absolute error             = ',HF_mu_of_r_bielec_energy - HF_alpha_beta_bielec_energy
  print*,'relative error             = ',dabs(HF_mu_of_r_bielec_energy - HF_alpha_beta_bielec_energy)/dabs(HF_alpha_beta_bielec_energy)
 endif


end


subroutine print_contribution_dft_mu_of_r
 implicit none

 print*,  '****************************************'
 print*, 'Functional used   = ',md_correlation_functional
 print*,  '****************************************'
 print*, 'mu_of_r_potential = ',mu_of_r_potential
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 if(md_correlation_functional.EQ."basis_set_LDA")then
   print*, ''
   write(*, '(A28,X,F16.10)') 'DFT mu(r)     correlation = ',Energy_c_md_mu_of_r_LDA
   print*, ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else if(md_correlation_functional.EQ."basis_set_on_top_PBE")then
   write(*, '(A28,X,F16.10)') 'Energy ECMD UEG        = ',Energy_c_md_on_top_PBE_mu_of_r_UEG
   write(*, '(A28,X,F16.10)') 'Energy ECMD NO UEG     = ',Energy_c_md_on_top_PBE_mu_of_r
   print*,''
   write(*, '(A28,X,F16.10)') 'Energy ECMD LDA        = ',Energy_c_md_mu_of_r_LDA
 endif 
  if(.true.)then
   write(*, '(A28,X,F16.10)') 'mu_average for basis set  = ',mu_average
  endif

end



subroutine print_variational_energy_dft_no_ecmd
 implicit none

 print*,  '****************************************'
 print*,  ' Regular range separated DFT energy '
 write(*, '(A22,X,F32.10)') 'mu_erf              = ',mu_erf          
 write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',electronic_energy_mr_dft+nuclear_repulsion
 print*, ''
 print*, 'Component of the energy ....'
 print*, ''
 write(*, '(A22,X,F16.10)') 'nuclear_repulsion   = ',nuclear_repulsion
 write(*, '(A22,X,F16.10)') 'psi_energy_erf      = ',psi_energy_erf      
 write(*, '(A22,X,F16.10)') 'psi_energy_core     = ',psi_energy_core    
 write(*, '(A22,X,F16.10)') 'short_range_Hartree = ',short_range_Hartree
 write(*, '(A22,X,F16.10)') 'two_elec_energy     = ',two_elec_energy_dft
 write(*, '(A22,X,F16.10)') 'energy_x            = ',energy_x         
 write(*, '(A22,X,F16.10)') 'energy_c            = ',energy_c          
 print*, ''
 print*,  '****************************************'
 print*, ''
 write(*, '(A22,X,F16.10)') 'Approx eigenvalue   = ',Fock_matrix_expectation_value + psi_energy_erf
 write(*, '(A22,X,F16.10)') 'Trace_v_xc          = ',Trace_v_xc
 print*, ''
 print*,  '****************************************'


end




subroutine print_projected_energy_dft
 implicit none
 print*,  '****************************************'
 write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',electronic_energy_mr_dft 
 print*, ''
 print*, 'Component of the energy ....'
 print*, ''
 write(*, '(A22,X,F16.10)') 'nuclear_repulsion   = ',nuclear_repulsion
 write(*, '(A22,X,F16.10)') 'projected energy    = ',data_energy_proj + short_range_Hartree + energy_x + energy_c - Trace_v_Hxc
 write(*, '(A22,X,F16.10)') 'short_range_Hartree = ',short_range_Hartree 
 write(*, '(A22,X,F16.10)') 'energy_x            = ',energy_x 
 write(*, '(A22,X,F16.10)') 'energy_c            = ',energy_c 
 write(*, '(A22,X,F16.10)') '- Trace_v_Hxc       = ',- Trace_v_Hxc


end


 BEGIN_PROVIDER [double precision, psi_dft_energy_kinetic, (N_states) ]
&BEGIN_PROVIDER [double precision, psi_dft_energy_nuclear_elec, (N_states) ]
 implicit none
 integer :: i,j,istate
 psi_dft_energy_kinetic = 0.d0
 psi_dft_energy_nuclear_elec = 0.d0 
 do istate = 1, N_states 
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    psi_dft_energy_kinetic(istate)      += ( one_body_dm_mo_alpha(j,i,istate)+one_body_dm_mo_beta(j,i,istate)) * mo_kinetic_integral(j,i) 
    psi_dft_energy_nuclear_elec(istate) += ( one_body_dm_mo_alpha(j,i,istate)+one_body_dm_mo_beta(j,i,istate)) * mo_nucl_elec_integral(j,i) 
   enddo
  enddo
 enddo



END_PROVIDER 
