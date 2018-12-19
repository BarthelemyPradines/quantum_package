double precision function f_HF_aa(r1,r2)
 implicit none
 BEGIN_DOC
! f_HF_aa(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha alpha spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > = 0.5 * \int (X1,X2) f_HF_aa(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 f_HF_aa = 0.d0
 do j = 1, elec_alpha_num ! electron 2
  do i = 1, elec_alpha_num ! electron 1 
   do l = 1, mo_tot_num ! electron 2 
    do k = 1, mo_tot_num ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_aa += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r2(j) & 
                * (mos_array_r1(k) * mos_array_r2(l) - mos_array_r1(l) * mos_array_r2(k)) ! direct and exchange term
    enddo
   enddo
  enddo
 enddo
end

double precision function f_HF_aa_integrated(r1)
 implicit none
 BEGIN_DOC
! f_HF_aa_integrated(X_1) = function int(X_2) f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha alpha spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > = 0.5 * \int (X1) f_HF_aa_integrated(X1)
 END_DOC
 double precision, intent(in) :: r1(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1) 
 f_HF_aa_integrated = 0.d0
 do j = 1, elec_alpha_num ! electron 2
  do i = 1, elec_alpha_num ! electron 1 
   do l = j, j ! electron 2 
    do k = 1, mo_tot_num ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_aa_integrated += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r1(k)  
    enddo
   enddo
  enddo
 enddo

 do j = 1, elec_alpha_num ! electron 2
  do i = 1, elec_alpha_num ! electron 1 
   do l = 1, mo_tot_num ! electron 2 
    do k = j, j ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_aa_integrated -= integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r1(l)  
    enddo
   enddo
  enddo
 enddo

end


double precision function f_HF_aa_(r1,r2)
 implicit none
 BEGIN_DOC
! f_HF_aa(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha alpha spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > = 0.5 * \int (X1,X2) f_HF_aa(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 f_HF_aa = 0.d0
 do j = 1, elec_alpha_num ! electron 2
  do i = 1, elec_alpha_num ! electron 1 
   do l = 1, mo_tot_num ! electron 2 
    do k = 1, mo_tot_num ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_aa += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r2(j) & 
                * (mos_array_r1(k) * mos_array_r2(l) - mos_array_r1(l) * mos_array_r2(k)) ! direct and exchange term
    enddo
   enddo
  enddo
 enddo
end

