double precision function two_dm_in_r(r1,r2,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 integer :: i,j,k,l
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 two_dm_in_r = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     if(dabs(two_bod_alpha_beta_mo_transposed(l,k,j,i,istate)).gt.1.d-10)then
    !  print*,l,k,j,i
     ! print*,two_bod_alpha_beta_mo_transposed(l,k,j,i,istate)
     endif
     two_dm_in_r += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * mos_array_r1(i) * mos_array_r1(l) * mos_array_r2(k) * mos_array_r2(j)
    enddo
   enddo
  enddo
 enddo
 two_dm_in_r = max(two_dm_in_r,1.d-15)
end

double precision function on_top_two_dm_in_r_mu_corrected(mu,r,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r(3),mu
 double precision :: mos_array_r(mo_tot_num),pi
 integer :: i,j,k,l
 double precision :: accu1,accu2,accu3,accu4,threshold
 threshold = 1.d-10
 pi = 4d0 * datan(1d0)
 call give_all_mos_at_r(r,mos_array_r) 

 on_top_two_dm_in_r_mu_corrected = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     on_top_two_dm_in_r_mu_corrected += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * mos_array_r(i) * mos_array_r(l) * mos_array_r(k) * mos_array_r(j)
    enddo
   enddo
  enddo
 enddo
 on_top_two_dm_in_r_mu_corrected = on_top_two_dm_in_r_mu_corrected / ( 1d0 + 2d0/(dsqrt(pi)*mu) )
 on_top_two_dm_in_r_mu_corrected = max(on_top_two_dm_in_r_mu_corrected,1.d-15)
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function on_top_dm_integral_with_mu_correction(mu,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: mu
 double precision :: two_dm_in_r, pi, r(3)
 double precision :: weight
 integer :: j,k,l
 pi = 4d0 * datan(1d0)
 on_top_dm_integral_with_mu_correction = 0d0

 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    on_top_dm_integral_with_mu_correction += two_dm_in_r(r,r,istate) * weight
   enddo
  enddo
 enddo
 on_top_dm_integral_with_mu_correction = 2d0 * on_top_dm_integral_with_mu_correction / ( 1d0 + 2d0/(dsqrt(pi)*mu) )

end


 BEGIN_PROVIDER [double precision, Energy_c_md_on_top, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density.
  ! Ec_md_on_top = (alpha/mu**3) * int n2(r,r) dr  where alpha = (sqrt(2pi)*(-2+sqrt(2)))/(3mu**3) 
 END_DOC
 implicit none 
 integer :: istate
 double precision :: pi,mu
 double precision :: on_top_dm_integral_with_mu_correction 
 mu = mu_erf
 pi = 4d0 * datan(1d0)
 do istate = 1, N_states
  Energy_c_md_on_top(istate) = ((-2d0+sqrt(2d0))*sqrt(2d0*pi)/(3d0*(mu**3)))*on_top_dm_integral_with_mu_correction(mu,istate)
 enddo
 END_PROVIDER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine give_epsilon_c_md_on_top_PBE(mu,r,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states)
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_c_md_on_top_PBE = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*two_dm_in_r(r,r,istate) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end



subroutine give_epsilon_c_md_on_top_PBE_and_corrected(mu,r,on_top,eps_c_md_on_top_PBE,eps_c_md_on_top_PBE_corrected)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(in)  :: on_top(N_states)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states),eps_c_md_on_top_PBE_corrected(N_states)
  double precision              :: on_top_corrected(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states)
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)
  double precision :: on_top_tmp(N_states)
  on_top_tmp = max(on_top,1.d-15)
  eps_c_md_on_top_PBE = 0d0
  eps_c_md_on_top_PBE_corrected = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*on_top_tmp(istate) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
   on_top_corrected(istate) = on_top_tmp(istate) / ( 1d0 + 2d0/(dsqrt(pi)*mu) )
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*on_top_corrected(istate) )
   eps_c_md_on_top_PBE_corrected(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end




 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision :: eps_c_md_on_top_PBE(N_states)
 double precision :: two_dm_in_r, r(3)
 double precision :: weight,mu
 integer :: j,k,l,istate
 double precision :: wall1,wall0  
!call cpu_time(wall0)
 mu = mu_erf
 Energy_c_md_on_top_PBE = 0d0
  
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    call give_epsilon_c_md_on_top_PBE(mu,r,eps_c_md_on_top_PBE)
    do istate = 1, N_states
     Energy_c_md_on_top_PBE(istate) += eps_c_md_on_top_PBE(istate) * weight
    enddo
   enddo
  enddo
 enddo
!call cpu_time(wall1)
!print*,'cpu time for Energy_c_md_on_top_PBE       '
!print*,wall1 - wall0
 
 END_PROVIDER


 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_cycle, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision :: eps_c_md_on_top_PBE(N_states)
 double precision :: two_dm_in_r, r(3)
 double precision :: weight,mu
 integer :: j,k,l,istate
 double precision :: dm_a,dm_b
 double precision :: wall1,wall0  
 call cpu_time(wall0)
 mu = mu_erf

 Energy_c_md_on_top_PBE_cycle = 0d0
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
    if(weight * one_body_dm_mo_alpha_at_grid_points(l,k,j,1).lt.threshold_grid_dft)cycle
    call give_epsilon_c_md_on_top_PBE(mu,r,eps_c_md_on_top_PBE)
    do istate = 1, N_states
     Energy_c_md_on_top_PBE_cycle(istate) += eps_c_md_on_top_PBE(istate) * weight
    enddo
   enddo
  enddo
 enddo

 call cpu_time(wall1)
 print*,'cpu time for Energy_c_md_on_top_PBE_cycle '
 print*,wall1 - wall0
 
 END_PROVIDER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine give_epsilon_c_md_on_top_PBE_mu_corrected(mu,r,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_c_md_on_top_PBE = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*on_top_two_dm_in_r_mu_corrected(mu,r,istate) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end
 

 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_corrected, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision :: eps_c_md_on_top_PBE(N_states)
 double precision :: two_dm_in_r, r(3)
 double precision :: weight,mu
 integer :: j,k,l,istate
 mu = mu_erf
 Energy_c_md_on_top_PBE_mu_corrected = 0d0
  
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    call give_epsilon_c_md_on_top_PBE_mu_corrected(mu,r,eps_c_md_on_top_PBE)
    do istate = 1, N_states
     Energy_c_md_on_top_PBE_mu_corrected(istate) += eps_c_md_on_top_PBE(istate) * weight
!    print*, Energy_c_md_on_top_PBE_mu_corrected(1)
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


