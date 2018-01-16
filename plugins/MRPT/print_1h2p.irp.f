program print_1h2p
 implicit none
 read_wf = .True.
 touch read_wf
 call routine_1h1p_bis
end


subroutine routine_1h1p_bis
 implicit none
 double precision,allocatable :: matrix_1h1p(:,:,:) 
 allocate (matrix_1h1p(N_det_ref,N_det_ref,N_states))
 integer :: i,j,istate,k
 double precision :: accu_diag_h_apply,accu_of_diag_h_apply
 double precision :: accu_diag_dm(N_states),accu_of_diag_dm(N_states)
 double precision :: accu_bis(N_states)
 double precision :: pt2(N_states)
 integer :: a,b,ispin,jspin,i_a,i_b,i_c,c,i_d,d
 matrix_1h1p = 0.d0
 call H_apply_mrpt_1h1p(matrix_1h1p,N_det_ref)
 do istate = 1, N_states
 accu_diag_h_apply = 0.d0
 accu_of_diag_h_apply = 0.d0
 do i = 1, N_det_ref
! write(*,'(100(F10.5,X))')matrix_1h1p(i,:,istate)
  accu_diag_h_apply += matrix_1h1p(i,i,istate) * psi_ref_coef(i,istate) * psi_ref_coef(i,istate)
  do j = 1, N_det_ref 
   if(i==j)cycle
   accu_of_diag_h_apply+= matrix_1h1p(i,j,istate) * psi_ref_coef(i,istate) * psi_ref_coef(j,istate)
  enddo
 enddo
 print*, 'h_apply_diag   =', accu_diag_h_apply
 print*, 'h_apply_of_di  =', accu_of_diag_h_apply
 print*, 'Total h_apply  =', accu_diag_h_apply + accu_of_diag_h_apply
 enddo
 print*,'******************************************************'

 integer :: occ_act(N_int*bit_kind_size,2),n_elec_act(2)
 accu_diag_dm = 0.d0
 
 do istate = 1, N_states
  do ispin = 1, 2
   do i_a = 1, n_act_orb
     i_b  = i_a
     accu_diag_dm(istate) += effective_active_energies_double_1h1p(i_a,ispin,istate) * cas_one_body_dm(i_b,i_a,ispin,istate)
     do i_b = 1, n_act_orb
      accu_diag_dm(istate) -= effective_coulomb_double_1h1hp(i_b,i_a,ispin,istate) * diag_cas_two_body_exchage_dm(i_b,i_a,ispin,istate)
     enddo
   enddo
  enddo
 enddo
 
 print*, 'accu_diag_dm   =', accu_diag_dm

 print*, 'accu_of_diag_dm=', accu_of_diag_dm
 print*, 'Total dm       =',accu_diag_dm+accu_of_diag_dm



end

subroutine routine_1h1p
 implicit none
 double precision,allocatable :: matrix_1h1p(:,:,:) 
 allocate (matrix_1h1p(N_det_ref,N_det_ref,N_states))
 integer :: i,j,istate,k
 double precision :: accu_diag_h_apply,accu_of_diag_h_apply
 double precision :: accu_diag_dm(N_states),accu_of_diag_dm(N_states)
 double precision :: accu_bis(N_states)
 double precision :: pt2(N_states)
 integer :: a,b,ispin,jspin,i_a,i_b,i_c,c
 matrix_1h1p = 0.d0
 call H_apply_mrpt_1h1p(matrix_1h1p,N_det_ref)
 do istate = 1, N_states
 accu_diag_h_apply = 0.d0
 accu_of_diag_h_apply = 0.d0
 do i = 1, N_det_ref
! write(*,'(100(F10.5,X))')matrix_1h1p(i,:,istate)
  accu_diag_h_apply += matrix_1h1p(i,i,istate) * psi_ref_coef(i,istate) * psi_ref_coef(i,istate)
  do j = 1, N_det_ref 
   if(i==j)cycle
   accu_of_diag_h_apply+= matrix_1h1p(i,j,istate) * psi_ref_coef(i,istate) * psi_ref_coef(j,istate)
  enddo
 enddo
 print*, 'h_apply_diag ', accu_diag_h_apply
 print*, 'h_apply_of_di', accu_of_diag_h_apply
!print*, 'H apply - density matrix = ',accu - pt2(istate)
 enddo

 integer :: occ_act(N_int*bit_kind_size,2),n_elec_act(2)
 accu_diag_dm = scalar_core_inact_contrib_1h1p
 
 do istate = 1, N_states
  do ispin = 1, 2
   do i_a = 1, n_act_orb
    accu_diag_dm(istate) += effective_active_energies_1h1p(i_a,istate) * cas_one_body_dm(i_a,i_a,ispin,istate)
    do jspin = 1,2
     do i_b = 1, n_act_orb
      accu_diag_dm(istate) += effective_coulomb_1h1hp(i_a,i_b,ispin,jspin,istate) * diag_cas_two_body_dm(i_b,i_a,ispin,ispin,istate)
     enddo
    enddo
   enddo
  enddo
 enddo
 print*, '************************'
 print*, 'accu_diag_dm ', accu_diag_dm
 print*, '************************'

 accu_bis = 0.d0
 accu_of_diag_dm = 0.d0
 do istate =1, N_states
  do i_a = 1, n_act_orb
   do i_b = 1, n_act_orb
    if(i_a==i_b)cycle
    do ispin = 1, 2
     accu_of_diag_dm(istate) += effective_Fock_1h1hp(i_a,i_b,ispin,istate) * cas_one_body_dm(i_a,i_b,ispin,istate)
     do jspin = 1, 2
      do i_c = 1, n_act_orb
       accu_of_diag_dm(istate) +=  effective_pseudo_Fock_1h1hp(i_c,jspin,i_b,i_a,ispin,istate) * pseudo_diag_cas_two_body_dm(i_c,jspin,i_b,i_a,ispin,istate)
       
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 print*, 'accu_of_diag ', accu_of_diag_h_apply

 call test_1h1p(pt2)





!accu_bis = 0.d0
!do istate = 1, N_states
! do i = 1, N_det_ref
!  do ispin = 1, 2
!   call bitstring_to_list(psi_active(1,ispin,i), occ_act(1,ispin), n_elec_act(ispin), N_int)
!  enddo
!  do ispin = 1, 2
!   do j = 1, n_elec_act(ispin)
!    a = list_act_reverse(occ_act(j,ispin))
!    print*, 'j,a,ispin',j,a,ispin
!    print*, effective_active_energies_1h1p(a,istate)
!    accu_bis(istate) += effective_active_energies_1h1p(a,istate)
!    do jspin = 1,2
!     do k = 1, n_elec_act(jspin)
!      b = list_act_reverse(occ_act(k,jspin))
!      print*, 'k,b,jspin',k,b,jspin
!      print*, effective_coulomb_1h1hp(a,b,ispin,jspin,istate)

!      accu_bis(istate) += effective_coulomb_1h1hp(a,b,ispin,jspin,istate)
!     enddo
!    enddo
!   enddo
!  enddo
! enddo
!enddo
!do istate = 1, N_states
! print*, 'accu_bis(istate)',accu_bis(istate)
! print*, 'accu_bis+pt2    ',accu_bis(istate)+ pt2(istate)
!enddo

 
!accu_bis = 0.d0
end

subroutine routine_2h1p
 implicit none
 double precision,allocatable :: matrix_1h2p(:,:,:) 
 allocate (matrix_1h2p(N_det_ref,N_det_ref,N_states))
 integer :: i,j,istate
 double precision :: accu
 double precision :: accu_bis(N_states)
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   do istate = 1, N_states
    matrix_1h2p(i,j,istate) = 0.d0
   enddo
  enddo
 enddo
 double precision :: pt2(N_states)
 call contrib_2h1p_dm_based(pt2)
 print*, 'pt2 = ',pt2
 matrix_1h2p = 0.d0
 call H_apply_mrpt_2h1p(matrix_1h2p,N_det_ref)
 do istate = 1, N_states
 accu = 0.d0
 do i = 1, N_det_ref
  do j = 1, N_det_ref 
   accu += matrix_1h2p(i,j,istate) * psi_ref_coef(i,istate) * psi_ref_coef(j,istate)
  enddo
 enddo
 print*, 'second order ', accu
 print*, 'H apply - density matrix = ',accu - pt2(istate)
 enddo
 
 integer :: a,b,ispin
 accu_bis = 0.d0
 do istate = 1, N_states
  accu_bis(istate) += effective_fock_operator_pure_diag_2h1p(istate) 
  do ispin = 1, 2
  do i = 1, n_act_orb
   a = list_act(i)
   do j = 1, n_act_orb
    b = list_act(j)
    accu_bis(istate) += effective_fock_operator_2h1p(i,j,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
   enddo
  enddo
  enddo
  print*, 'accu_bis = ',accu_bis(istate)
  print*, 'H apply - density matrix = ',accu_bis(istate) - pt2(istate)
 enddo
end

subroutine routine_1h2p
 implicit none
 double precision,allocatable :: matrix_1h2p(:,:,:) 
 allocate (matrix_1h2p(N_det_ref,N_det_ref,N_states))
 integer :: i,j,istate
 double precision :: accu
 double precision :: accu_bis(N_states)
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   do istate = 1, N_states
    matrix_1h2p(i,j,istate) = 0.d0
   enddo
  enddo
 enddo

 call contrib_1h2p_dm_based(accu_bis)
 print*, 'accu_bis',accu_bis
 print*, '1h2p    ',contribution_1h2p
 matrix_1h2p = 0.d0
 call H_apply_mrpt_1h2p(matrix_1h2p,N_det_ref)
 do istate = 1, N_states
 accu = 0.d0
 do i = 1, N_det_ref
  do j = 1, N_det_ref 
   accu += matrix_1h2p(i,j,istate) * psi_ref_coef(i,istate) * psi_ref_coef(j,istate)
  enddo
 enddo
 print*, 'second order ', accu
 print*, 'H apply - density matrix = ',accu - accu_bis(istate)
 enddo
 accu = 0.d0
 integer :: ispin,a,b
 do istate = 1, N_states
  do ispin = 1, 2
   do i = 1, n_act_orb
    a = list_act(i)
    do j = 1, n_act_orb 
     b = list_act(j)
     accu += effective_fock_operator_1h2p(i,j,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
    enddo
   enddo
  enddo
  print*, 'accu = ',accu
 print*, 'H apply - density matrix = ',accu - accu_bis(istate)
 enddo


 deallocate (matrix_1h2p)
end
