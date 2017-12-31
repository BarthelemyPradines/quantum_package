
 BEGIN_PROVIDER [double precision, cas_one_body_dm, (n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, diag_cas_two_body_dm, (n_act_orb,n_act_orb,2,2,N_states)]
&BEGIN_PROVIDER [double precision, pseudo_diag_cas_two_body_dm, (n_act_orb,2,n_act_orb,n_act_orb,2,N_states)]
 implicit none
 BEGIN_DOC
 ! diag_cas_two_body_dm(i_a,i_b,ispin,jspin) = \sum_{I} c_I^2 <I| \hat{n}_{a,ispin}  \hat{n}_{b,jspin} |I>
 ! pseudo_diag_cas_two_body_dm(i_a,ispin,i_b,i_c,jspin) = \sum_{I,J} c_I c_J <J| a^{\dagger}_{c,jspin} a_{b,jspin} \hat{n}_{a,ispin} |I>
 END_DOC
 integer :: idet,jdet,j,k,i_a,i_b,a,b,istate,ispin, jspin
 integer :: list_active(N_int*bit_kind_size,2),n_elec_active(2)
 double precision :: accu(N_states)
 cas_one_body_dm = 0.d0
 diag_cas_two_body_dm = 0.d0
 pseudo_diag_cas_two_body_dm = 0.d0
 do idet = 1, N_det_ref
  call bitstring_to_list(psi_active(1,1,idet), list_active(1,1), n_elec_active(1), N_int)
  call bitstring_to_list(psi_active(1,2,idet), list_active(1,2), n_elec_active(2), N_int)
  do istate = 1, N_states
   accu(istate) = psi_ref_coef(idet,istate) * psi_ref_coef(idet,istate)
  enddo
  do ispin = 1, 2
   do j = 1, n_elec_active(ispin)
    i_a = list_act_reverse(list_active(j,ispin))
    do istate = 1, N_states
     cas_one_body_dm(i_a,i_a,ispin,istate) += accu(istate)
    enddo
   enddo
   do jspin = 1, 2
    do j = 1, n_elec_active(ispin)
     i_a = list_act_reverse(list_active(j,ispin))
     do k = 1, n_elec_active(jspin)
      i_b = list_act_reverse(list_active(k,jspin))
      do istate = 1, N_states
       diag_cas_two_body_dm(i_b,i_a,jspin,ispin,istate) += accu(istate)
      enddo
     enddo
    enddo
   enddo
  enddo
  integer           :: degree(N_det_ref)
  integer           :: idx(0:N_det_ref)
  integer :: exc(0:2,2,2)
  double precision :: phase
  integer :: s1,h1,p1
  call get_excitation_degree_vector_mono(psi_ref,psi_ref(1,1,idet),degree,N_int,N_det_ref,idx)
  do jdet = 1, idx(0)
   if(idx(jdet).eq.idet)cycle
   call get_mono_excitation(psi_ref(1,1,idet),psi_ref(1,1,idx(jdet)),exc,phase,N_int)
   if(exc(0,1,1)==1)then
    s1 = 1
    h1 = list_act_reverse(exc(1,1,1))
    p1 = list_act_reverse(exc(1,2,1))
   else
    s1 = 2
    h1 = list_act_reverse(exc(1,1,2))
    p1 = list_act_reverse(exc(1,2,2))
   endif
   do istate = 1, N_states
    cas_one_body_dm(h1,p1,s1,istate) += psi_ref_coef(idet,istate) * psi_ref_coef(idx(jdet),istate) * phase
    do ispin = 1, 2
     do j = 1, n_elec_active(ispin)
      i_a = list_act_reverse(list_active(j,ispin))
      pseudo_diag_cas_two_body_dm(i_a,ispin,h1,p1,s1,istate) += psi_ref_coef(idet,istate) * psi_ref_coef(idx(jdet),istate) * phase
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 


