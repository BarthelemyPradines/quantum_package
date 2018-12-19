program projected_operators
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
! call routine_v
! call routine_rho 
! call routine_final
!call print_energy
 call routine_v

end

subroutine print_energy
 implicit none
 print*,'psi_energy_bielec_ab = ',psi_energy_bielec_ab
 print*,'psi_energy_bielec    = ',psi_energy_bielec
 print*,'psi_energy           = ',psi_energy - psi_energy_monoelec

end

subroutine routine_rho
 implicit none
 provide rho2_kl_contracted_transposed

end
subroutine routine_v
 implicit none
 integer :: ipoint,k,l
 double precision :: accu, weight,f,f_HF_aa_integrated,r(3)
 accu = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_functions_at_final_grid_points(ipoint)
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  f = f_HF_aa_integrated(r)
  accu += f * weight 
 enddo
 print*,'accu*0.5          = ',accu*0.5d0
 print*,'psi_energy_bielec = ',psi_energy_bielec
 print*,'psi_energy_bielec_aa',psi_energy_bielec_aa
end
 

subroutine routine_rho2
 implicit none
 integer :: ipoint,k,l
 double precision :: accu, weight
 accu = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_functions_at_final_grid_points(ipoint)
  do l = 1, mo_tot_num
   do k = 1, mo_tot_num
    accu += dabs(rho2_kl_contracted(k,l,ipoint) - rho2_kl_contracted_sequential(k,l,ipoint)) * weight
   enddo
  enddo
 enddo
 print*,'accu = ',accu
end
 



subroutine routine_final
 implicit none
 integer :: ipoint
 double precision :: accu, weight
 accu = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_functions_at_final_grid_points(ipoint)
! accu += dabs(f_psi_ab_old(ipoint) - f_psi_ab(ipoint)) *weight
  accu += (f_psi_ab(ipoint)) *weight
 enddo
 print*,'*******************'
 print*,'*******************'
 print*,'*******************'
 print*,'*******************'
 print*,'accu = ',accu


end
