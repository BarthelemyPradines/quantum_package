BEGIN_PROVIDER [double precision, V_kl_contracted_transposed, (n_points_final_grid,mo_tot_num,mo_tot_num)]
 implicit none
 BEGIN_DOC
! V_kl_contracted_transposed(ipoint,k,l) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l
 integer :: ipoint
 double precision, allocatable :: integrals_array(:,:), mos_array_r(:),r(:)
 ! just not to mess with parallelization
 allocate(integrals_array(mo_tot_num,mo_tot_num))
  k = 1
  l = 1
  call get_mo_bielec_integrals_ij(k,l,mo_tot_num,integrals_array,mo_integrals_map) 
 deallocate(integrals_array)
 double precision :: wall0,wall1
 call wall_time(wall0)

 print*,'Providing  V_kl_contracted_transposed ..... '
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,k,l,i,j,integrals_array) & 
 !$OMP SHARED (mo_tot_num, n_points_final_grid, V_kl_contracted_transposed, mo_integrals_map,final_grid_points,mos_in_r_array)
 allocate(integrals_array(mo_tot_num,mo_tot_num))
 !$OMP DO              
  do l = 1, mo_tot_num ! 2 
   do k = 1, mo_tot_num ! 1 
    call get_mo_bielec_integrals_ij(k,l,mo_tot_num,integrals_array,mo_integrals_map) 
    V_kl_contracted_transposed(k,l,:) = 0.d0
    do ipoint = 1, n_points_final_grid
     do j = 1, mo_tot_num
      do i = 1, mo_tot_num
                                        !1 2                     1 2 
       V_kl_contracted_transposed(ipoint,k,l) += integrals_array(i,j) * mos_in_r_array(j,ipoint) * mos_in_r_array(i,ipoint)
      enddo
     enddo
    enddo
   enddo
  enddo
 !$OMP END DO
 deallocate(integrals_array)
 !$OMP END PARALLEL

 call wall_time(wall1)
 print*,'Time to provide V_kl_contracted_transposed = ',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [double precision, V_kl_contracted, (mo_tot_num,mo_tot_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! V_kl_contracted(k,l,ipoint) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: ipoint,k,l
 do k = 1, mo_tot_num
  do l = 1, mo_tot_num
   do ipoint = 1, n_points_final_grid
    V_kl_contracted(k,l,ipoint) = V_kl_contracted_transposed(ipoint,k,l)
   enddo
  enddo
 enddo
 free V_kl_contracted_transposed 

END_PROVIDER 

BEGIN_PROVIDER [double precision, rho2_kl_contracted_transposed, (n_points_final_grid,mo_tot_num,mo_tot_num)]
 implicit none
 BEGIN_DOC
! rho2_kl_contracted_transposed(k,l,ipoint) = \sum_{ij} rho2_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l
 integer :: ipoint
 double precision, allocatable :: mos_array_r(:),r(:)
 provide two_bod_alpha_beta_mo_physician
 double precision :: wall0,wall1
 print*,'Providing  rho2_kl_contracted_transposed ..... '
 call wall_time(wall0)
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,k,l,i,j) & 
 !$OMP SHARED  (mo_tot_num, n_points_final_grid, rho2_kl_contracted_transposed, final_grid_points,two_bod_alpha_beta_mo_physician,mos_in_r_array )
 !$OMP DO              
 do l = 1, mo_tot_num ! 2 
  do k = 1, mo_tot_num ! 1 
   rho2_kl_contracted_transposed(k,l,:) = 0.d0
   do ipoint = 1, n_points_final_grid
    do j = 1, mo_tot_num
     do i = 1, mo_tot_num
                                          !1 2                                     1 2 1 2 
      rho2_kl_contracted_transposed(ipoint,k,l) += two_bod_alpha_beta_mo_physician(i,j,k,l,1) * mos_in_r_array(j,ipoint) * mos_in_r_array(i,ipoint)
     enddo
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide rho2_kl_contracted_transposed = ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [double precision, rho2_kl_contracted, (mo_tot_num,mo_tot_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! rho2_kl_contracted(k,l,ipoint) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: ipoint,k,l
 do k = 1, mo_tot_num
  do l = 1, mo_tot_num
   do ipoint = 1, n_points_final_grid
    rho2_kl_contracted(k,l,ipoint) = rho2_kl_contracted_transposed(ipoint,k,l)
   enddo
  enddo
 enddo
 free rho2_kl_contracted_transposed 
END_PROVIDER 


BEGIN_PROVIDER [double precision, f_psi_ab, (n_points_final_grid)]
 implicit none
 BEGIN_DOC
! f_PSI_ab(X1,X1) = function f_{\Psi^B}(X_1,X_1) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha alpha spins and a general wave function
! < HF | wee_{\alpha\beta} | HF > = 0.5 * \int (X1,X2) f_HF_ab(X1,X2)
 END_DOC
 integer :: ipoint
 integer :: k,l 
 double precision :: wall0,wall1
 provide two_bod_alpha_beta_mo_physician 
 print*,'Providing  f_psi_ab ..... '
 call wall_time(wall0)
 provide V_kl_contracted rho2_kl_contracted
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,k,l) & 
 !$OMP SHARED  (mo_tot_num, n_points_final_grid, rho2_kl_contracted, V_kl_contracted, f_psi_ab)
 !$OMP DO              
 do ipoint = 1, n_points_final_grid
  f_psi_ab(ipoint) = 0.d0
  do l = 1, mo_tot_num ! 2 
   do k = 1, mo_tot_num ! 1
    f_psi_ab(ipoint) += V_kl_contracted(k,l,ipoint) * rho2_kl_contracted(k,l,ipoint)
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide f_psi_ab = ',wall1 - wall0
END_PROVIDER 


BEGIN_PROVIDER [double precision, f_psi_ab_old, (n_points_final_grid)]
 implicit none
 integer :: ipoint
 double precision :: r(3),coulomb,two_body_dm
 integer :: k,l 
  r = 0.d0
 double precision :: wall0,wall1
 print*,'Providing  f_psi_ab_old ..... '
 provide two_bod_alpha_beta_mo_physician 
 call wall_time(wall0)
 call expectation_value_in_real_space_no_divide(r,r,coulomb,two_body_dm)
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,r,coulomb,two_body_dm) & 
 !$OMP SHARED  (n_points_final_grid, f_psi_ab_old, final_grid_points)
 !$OMP DO              
 do ipoint = 1, n_points_final_grid
  f_psi_ab_old(ipoint) = 0.d0
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call expectation_value_in_real_space_no_divide(r,r,coulomb,two_body_dm)
  f_psi_ab_old(ipoint) = coulomb
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide f_psi_ab_old = ',wall1 - wall0
END_PROVIDER 

