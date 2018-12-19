subroutine f_HF_ab(r1,r2,integral_psi,two_bod)
 implicit none
 BEGIN_DOC
! f_HF_ab(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha beta spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > = 0.5 * \int (X1,X2) f_HF_aa(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral_psi,two_bod
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
 do m = 1, elec_alpha_num
  do n = 1, elec_beta_num
   two_bod += mos_array_r1(n) * mos_array_r1(n) * mos_array_r2(m) * mos_array_r2(m) 
   do i = 1, mo_tot_num
    do j = 1, mo_tot_num
     integral_psi +=  integrals_for_hf_potential(j,i,n,m) * mos_array_r1(i) * mos_array_r2(j) * mos_array_r2(n) * mos_array_r1(m) 
    enddo
   enddo
  enddo
 enddo
end


BEGIN_PROVIDER [double precision, integrals_for_hf_potential, (mo_tot_num,mo_tot_num,elec_alpha_num,elec_alpha_num)]
 implicit none
 BEGIN_DOC
! integrals_for_hf_potential(k,l,i,j) = <ij|kl> where k,l runs over all orbitals and i,j until elec_alpha_num
 END_DOC
 integer :: i,j,m,n
 double precision :: get_mo_bielec_integral
 do m = 1, elec_alpha_num ! electron 1 alpha 
  do n = 1, elec_alpha_num ! electron 2 beta 
   do i = 1, mo_tot_num   ! electron 1 alpha
    do j = 1, mo_tot_num  ! electron 2 beta 
                               !2 1 2 1                           1 2 1 2
     integrals_for_hf_potential(j,i,n,m) = get_mo_bielec_integral(m,n,i,j,mo_integrals_map) 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

