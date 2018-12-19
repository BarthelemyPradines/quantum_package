subroutine coulomb_operator_in_real_space(r1,r2,coulomb)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb
 integer :: i,j,k,l 
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3
 double precision :: c1,c2,c3
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral
 threshold = 1.d-12

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

!call give_all_mos_at_r_old(r1,mos_array_r1) 
!call give_all_mos_at_r_old(r2,mos_array_r2) 

 coulomb = 0.d0
 do i = 1, mo_tot_num
  a1 = mos_array_r1(i) 
  c1 = dabs(a1)
  if(c1.le.threshold)cycle
  do j = 1, mo_tot_num 
   a2 = a1 * mos_array_r2(j)
   c2 = dabs(a2)
   coulomb += coulomb_integrals_contracted_with_overlap(j,i) * a2 
   enddo
  enddo
 !coulomb = coulomb * 0.5d0 
end

BEGIN_PROVIDER [double precision, coulomb_integrals_contracted_with_overlap, (mo_tot_num,mo_tot_num)]
 implicit none
 BEGIN_DOC
! coulomb_integrals_contracted_with_overlap(i,j) = \sum_{kl} <ij|kl> * (\int_r1 phi_k(r1)) * (\int_r2 phi_l(r2))
 END_DOC
 integer :: i,j,k,l
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 coulomb_integrals_contracted_with_overlap = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num 
   call get_mo_bielec_integrals_ij(i,j,mo_tot_num,integrals_array,mo_integrals_map) 
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      coulomb_integrals_contracted_with_overlap(j,i)+=  integrals_array(k,l) * integrals_mo_over_r(k) * integrals_mo_over_r(l)
     enddo
    enddo
   enddo
  enddo


END_PROVIDER 


subroutine erf_coulomb_operator_in_real_space(r1,r2,coulomb)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb
 integer :: i,j,k,l 
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3
 double precision :: c1,c2,c3
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral_erf
 threshold = 1.d-12

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

!call give_all_mos_at_r_old(r1,mos_array_r1) 
!call give_all_mos_at_r_old(r2,mos_array_r2) 

 coulomb = 0.d0
 do i = 1, mo_tot_num
  a1 = mos_array_r1(i) 
  c1 = dabs(a1)
  if(c1.le.threshold)cycle
  do j = 1, mo_tot_num 
   a2 = a1 * mos_array_r2(j)
   c2 = dabs(a2)
   if(c2.le.threshold)cycle
    do k = 1, mo_tot_num
     a3 = a2 * mos_array_r1(k)
     c3 = dabs(a3)
     if(c3.le.threshold)cycle
     do l = 1, mo_tot_num
      integral = get_mo_bielec_integral_erf(i,j,k,l,mo_integrals_erf_map)
      coulomb += a3 * mos_array_r2(l) * integral
!     coulomb += mos_array_r1(i) * mos_array_r1(k) * mos_array_r2(j) * mos_array_r2(l) * integrals_array(k,l)
     enddo
    enddo
   enddo
  enddo
 coulomb = coulomb * 0.5d0 

end

subroutine nuclear_coulomb_operator_in_real_space(r1,coulomb)
 implicit none
 double precision, intent(in) :: r1(3)
 double precision, intent(out):: coulomb
 integer :: i,j 
 double precision :: mos_array_r1(mo_tot_num)

 call give_all_mos_at_r(r1,mos_array_r1) 
 coulomb = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num 
   coulomb += mo_nucl_elec_integral(j,i) * mos_array_r1(j) * integrals_mo_over_r(i)
  enddo
 enddo

end


subroutine expectation_value_in_real_space_ol(r1,r2,coulomb,two_dm)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb, two_dm
 integer :: i,j,k,l,m,n  
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3,a4
 double precision :: c1,c2,c3,c4
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral,mo_bielec_integral
 double precision :: two_dm_orb
 threshold = 0.d-10

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 two_dm = 0.d0
 coulomb = 0.d0
 do k = 1, mo_tot_num ! electron 1 
  do l = 1, mo_tot_num ! electron 2
   do m = 1, mo_tot_num ! electron 1 
    do n = 1, mo_tot_num ! electron 2 
     if(dabs(two_bod_alpha_beta_mo(m,k,n,l,1)).gt.1.d-10)then
     print*,'m,k,n,l = ',m,k,n,l
     print*,'1,1,2,2'
     print*,two_bod_alpha_beta_mo(m,k,n,l,1)
     endif
     do i = 1, mo_tot_num ! electron 1
      do j = 1, mo_tot_num ! electron 2 
       coulomb += two_bod_alpha_beta_mo(m,k,n,l,1) * mos_array_r1(i) * mos_array_r1(m) * mos_array_r2(j) * mos_array_r2(n) * get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


 do i = 1, mo_tot_num
  do j = 1, mo_tot_num 
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     two_dm += two_bod_alpha_beta_mo_transposed(l,k,j,i,1) * mos_array_r1(i) * mos_array_r1(l)  * mos_array_r2(k) * mos_array_r2(j)  
    enddo
   enddo
  enddo
 enddo
 coulomb = coulomb / two_dm

end

subroutine expectation_value_in_real_space(r1,r2,coulomb,two_body_dm)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb, two_body_dm
 integer :: i,j,k,l,m,n  
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3,a4
 double precision :: c1,c2,c3,c4
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral
 threshold = 0.d0

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 two_body_dm = 0.d0
 coulomb = 0.d0
 double precision :: tmp1
 do i = 1, mo_tot_num ! 1
  do j = 1, mo_tot_num ! 2 
   do m = 1, mo_tot_num ! 1
    do n = 1, mo_tot_num ! 2
     !                                               2 1 2 1
     tmp1 = mos_array_r2(n) * mos_array_r1(m) * mos_array_r2(j) * mos_array_r1(i)
     coulomb     += two_bod_alpha_beta_mo_contracted(n,m,j,i,1) * tmp1 
     two_body_dm += two_bod_alpha_beta_mo_physician (n,m,j,i,1) * tmp1 
    enddo
   enddo
  enddo
 enddo
 if(two_body_dm.gt.1.d-12.and.coulomb.gt.1.d-12)then
  coulomb = coulomb/two_body_dm
 else 
  coulomb = 1.d-10
 endif

end


subroutine expectation_value_in_real_space_no_divide(r1,r2,coulomb,two_body_dm)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb, two_body_dm
 integer :: i,j,k,l,m,n  
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3,a4
 double precision :: c1,c2,c3,c4
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral
 threshold = 0.d0

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 two_body_dm = 0.d0
 coulomb = 0.d0
 double precision :: tmp1
 do i = 1, mo_tot_num ! 1
  do j = 1, mo_tot_num ! 2 
   do m = 1, mo_tot_num ! 1
    do n = 1, mo_tot_num ! 2
     !                                               2 1 2 1
     tmp1 = mos_array_r2(n) * mos_array_r1(m) * mos_array_r2(j) * mos_array_r1(i)
     coulomb     += two_bod_alpha_beta_mo_contracted(n,m,j,i,1) * tmp1 
     two_body_dm += two_bod_alpha_beta_mo_physician (n,m,j,i,1) * tmp1 
    enddo
   enddo
  enddo
 enddo

end


subroutine pure_expectation_value_in_real_space(r1,r2,coulomb)
use map_module
 implicit none
 BEGIN_DOC
 ! f(1,2) of equation (22) of the first paper in JCP
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb
 integer :: i,j,k,l,m,n  
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3,a4
 double precision :: c1,c2,c3,c4
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral
 threshold = 0.d0

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 coulomb = 0.d0
 do i = 1, mo_tot_num ! 1
  do j = 1, mo_tot_num ! 2 
   do m = 1, mo_tot_num ! 1
    do n = 1, mo_tot_num ! 2
     !                                           2 1 2 1
     coulomb += two_bod_alpha_beta_mo_contracted(n,m,j,i,1) * mos_array_r2(n) * mos_array_r1(m) * mos_array_r2(j) * mos_array_r1(i) 
    enddo
   enddo
  enddo
 enddo

end


double precision function integral_pure_expectation_value_in_real_space(r)
 implicit none
! integral over 2 of pure_expectation_value_in_real_space(1,2) 
 double precision, intent(in) :: r(3)
 double precision :: mos_array_r(mo_tot_num)
 integer :: i,j,m,n
 call give_all_mos_at_r(r,mos_array_r)
 integral_pure_expectation_value_in_real_space = 0.d0
 do i = 1, mo_tot_num ! 1
  do j = 1, mo_tot_num ! 2 
   do m = 1, mo_tot_num ! 1
    do n = j,j           ! 2 :: kronecker(n,j) because of orthogonality of MOs
     !                                                            2 1 2 1
     integral_pure_expectation_value_in_real_space += two_bod_alpha_beta_mo_contracted(n,m,j,i,1) * mos_array_r(m) * mos_array_r(i) 
    enddo
   enddo
  enddo
 enddo
end

subroutine numerical_delta_function_coulomb(r1,r2,integral)
 implicit none
 BEGIN_DOC
! compute numerically the following integral: 
! int_r \sum_m \phi_m(r)\phi_m(r2) * 1/|r - r1|
 END_DOC
 double precision,intent(in) :: r1(3),r2(3)
 double precision,intent(out):: integral
 double precision :: r(3),mos_array_r2(mo_tot_num),mos_array_r(mo_tot_num),distance
 integer :: j,k,l ,m
 integral = 0.d0
 call give_all_mos_at_r(r2,mos_array_r2) 
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(1) = grid_points_per_atom(1,l,k,j)
    r(2) = grid_points_per_atom(2,l,k,j)
    r(3) = grid_points_per_atom(3,l,k,j)
    distance  = (r1(1) - r(1))*(r1(1) - r(1))
    distance += (r1(2) - r(2))*(r1(2) - r(2))
    distance += (r1(3) - r(3))*(r1(3) - r(3))
    distance = dsqrt(distance)
    if(dabs(distance).lt.1.d-10)cycle
    call give_all_mos_at_r(r,mos_array_r) 
    do m = 1, mo_tot_num
     integral += mos_array_r(m) * mos_array_r2(m) * 1.d0/(distance) * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
 enddo
end


subroutine local_r12_operator_with_one_e_int(r1,r2,integral)
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral(mo_tot_num)
 integer :: i,j
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: mos_array_ints(mo_tot_num,mo_tot_num)
 BEGIN_DOC 
 ! computes the following integral for all mos 
 ! integral(i) = 1/phi_i(r2) * \sum_{j} phi_j(r2) * \int d_r phi_j(r) phi_i(r1) 1 / |r - r1| 
 END_DOC
 call potential_mono_elec_integral_mo_at_r(r1,mos_array_ints)
 call give_all_mos_at_r(r2,mos_array_r2) 
 do i = 1, mo_tot_num
  integral(i) = 0.d0
  do j = 1, mo_tot_num
   integral(i) += mos_array_r2(j) * mos_array_ints(j,i)
  enddo
  integral(i) = integral(i) / mos_array_r2(i)
 enddo

end

subroutine local_r12_operator_with_one_e_int_on_1s(r1,r2,integral)
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral
 integer :: i,j
 double precision :: mos_array_r2(mo_tot_num)
 BEGIN_DOC 
 ! computes the following integral for all mos 
 ! integral(i) = 1/phi_i(r2) * \sum_{j} phi_j(r2) * \int d_r phi_j(r) phi_i(r1) 1 / |r - r1| 
 END_DOC
 integral = 0.d0
 call give_all_mos_at_r(r2,mos_array_r2) 

 call potential_mono_elec_integral_mo_at_r(r1,mos_array_ints)
 double precision :: mos_array_ints(mo_tot_num,mo_tot_num)
 do j = 1, mo_tot_num
  integral += mos_array_r2(j) * mos_array_ints(j,1)
 enddo

!double precision :: mos_array_ints(mo_tot_num)
!call potential_mono_elec_integral_mo_at_r_on_1s(r1,mos_array_ints)
!do j = 1, mo_tot_num
! integral += mos_array_r2(j) * mos_array_ints(j)
!enddo


 integral = integral / mos_array_r2(1)

end

BEGIN_PROVIDER [double precision, integrals_for_hf_potential_integrated_on_beta, (mo_tot_num,elec_alpha_num)]
 implicit none
 integer :: i,j,k
 double precision :: get_mo_bielec_integral
 integrals_for_hf_potential_integrated_on_beta = 0.d0
 do k = 1, elec_alpha_num ! electron 1 alpha 
  do i = 1, mo_tot_num   ! electron 1 alpha
   do j = 1, elec_beta_num ! electron 2 beta 
    integrals_for_hf_potential_integrated_on_beta(i,k) += get_mo_bielec_integral(i,j,k,j,mo_integrals_map) 
   enddo
  enddo
 enddo
END_PROVIDER 


double precision function mu_coulomb(y,x)
 implicit none
 BEGIN_DOC
! finds the correct mu that such that 
! y = erf(mu_coulomb * x)/x 
 ENd_DOC
 double precision, intent(in) :: y,x
 double precision :: y_tmp,x_tmp,thr,inverse_erf
 thr = 1.d-10
 y_tmp = y * x
 x_tmp =inverse_erf(y_tmp,thr)
 if(x.ne.0.d0)then
  mu_coulomb = x_tmp / x
 else
  mu_coulomb = x_tmp / 1.d-15
 endif
end


subroutine local_r12_operator_on_hf_act(r1,r2,integral_psi)
 implicit none
 BEGIN_DOC
! computes the following ANALYTICAL integral
! \sum_{m,n,i,j} V_{ij}^{mn} phi_i(1) * phi_j(2) * phi_m(1) * phi_n(2) / rho_alpha^{HF}(1)*rho_beta^{HF}(2)
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral_psi
 integer :: i,j,m,n
 integer :: ii,jj,mm,nn
 double precision :: mo_bielec_integral,two_bod
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: get_mo_bielec_integral
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 integral_psi = 0.d0
 two_bod = 0.d0
 do mm = 1, n_act_orb
  m = list_act(mm)
  if(m>elec_alpha_num)cycle
  do nn = 1, n_act_orb
   n = list_act(nn)
   if(n>elec_beta_num)cycle
   two_bod += mos_array_r1(n) * mos_array_r1(n) * mos_array_r2(m) * mos_array_r2(m) 
   do ii = 1, n_act_orb
    i = list_act(ii)
    do jj = 1, n_act_orb
     j = list_act(jj)
     integral_psi +=  integrals_for_hf_potential(j,i,n,m) * mos_array_r1(i) * mos_array_r2(j) * mos_array_r2(n) * mos_array_r1(m) 
    enddo
   enddo
  enddo
 enddo
 integral_psi = integral_psi /  two_bod

end



subroutine integral_of_f_12_hf_over_beta(r1,integral_f)
 implicit none
 BEGIN_DOC
! computes the following ANALYTICAL integral
! integral_(r2,beta) f_HF(r1,r2) with f_HF is the function f of the first paper but defined for a HF two-body tensor 
 END_DOC
 double precision, intent(in) :: r1(3)
 double precision, intent(out):: integral_f
 integer :: i,k
 double precision :: mos_array_r1(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1)

 integral_f = 0.d0
 do k = 1, elec_alpha_num
  do i = 1, mo_tot_num
   integral_f += integrals_for_hf_potential_integrated_on_beta(i,k) * mos_array_r1(i) * mos_array_r1(k)
  enddo
 enddo

end

 BEGIN_PROVIDER [double precision, integral_f_hf ]
 implicit none
 BEGIN_DOC 
! should be the average value of the alpha/beta bielectronic repulsion over the HF wave function
 END_DOC
 double precision :: r(3),integral_f,weight
 integer :: i_point
 integral_f_hf= 0.D0
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call integral_of_f_12_hf_over_beta(r,integral_f)
  weight = final_weight_functions_at_final_grid_points(i_point)
  integral_f_hf += weight * integral_f
 enddo
 END_PROVIDER 


BEGIN_PROVIDER [double precision, hf_bielec_alpha_beta_energy]
 implicit none
 integer :: i,j
 hf_bielec_alpha_beta_energy = 0.d0
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   hf_bielec_alpha_beta_energy += integrals_for_hf_potential(j,i,j,i)
  enddo
 enddo
END_PROVIDER 
