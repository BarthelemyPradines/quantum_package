program print_ecmd_energy
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
!call ecmd_energy_printer
 call pouet 
end

subroutine pouet
 implicit none
 integer :: i,j,k,l,m
 double precision:: r(3)
  write(*, '(A22,X,F16.10)') 'energy_x            = ',energy_x
  write(*, '(A22,X,F16.10)') 'energy_c            = ',energy_c

  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)

!    write(34,'(100(F16.10,X))')r(:),one_body_dm_mo_alpha_and_grad_at_grid_points(:,l,k,j,1)
    enddo
   enddo
  enddo
!do i = 1, mo_tot_num
! write(33,'(1000(F16.10,X))')  one_body_dm_alpha_mo_for_dft(i,:,1)
!enddo

end



subroutine ecmd_energy_printer
 implicit none
 
 print*,  '****************************************'
 write(*, '(A22,X,F32.10)') 'mu_erf              = ',mu_erf          
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 write(*, '(A22,X,F16.10)') 'EC_MD_LDA           = ',Energy_c_md+psi_energy+nuclear_repulsion
!write(*, '(A22,X,F16.10)') 'EC_MD_ON_TOP        = ',Energy_c_md_on_top(1)+psi_energy+nuclear_repulsion
!write(*, '(A22,X,F16.10)') 'EC_MD_ON_TOP_PBE    = ',Energy_c_md_on_top_PBE(1)+psi_energy+nuclear_repulsion
 write(*, '(A22,X,F16.10)') 'EC_MD_ON_TOP_PBE_cor= ',Energy_c_md_on_top_PBE_mu_corrected_UEG(1)+psi_energy+nuclear_repulsion
 print*, ''
 print*, 'Component of the energy ....'
 print*, ''
 write(*, '(A28,X,F16.10)') 'nuclear_repulsion         = ',nuclear_repulsion
 write(*, '(A28,X,F16.10)') 'Variational energy of Psi = ',psi_energy
 write(*, '(A28,X,F16.10)') 'psi_energy_bielec         = ',psi_energy_bielec
 write(*, '(A28,X,F16.10)') 'psi_energy_monoelec       = ',psi_energy_monoelec
 write(*, '(A28,X,F16.10)') 'LDA Multi-det correlation = ',Energy_c_md
!write(*, '(A28,X,F16.10)') 'on_top Multi-det correl   = ',Energy_c_md_on_top(1)
!write(*, '(A28,X,F16.10)') 'on_top_PBE MD correl      = ',Energy_c_md_on_top_PBE(1)
 write(*, '(A28,X,F16.10)') 'on_top_PBE_cor MD correl  = ',Energy_c_md_on_top_PBE_mu_corrected_UEG(1)

end
