
use bitmasks
subroutine contrib_1h2p_dm_based(accu)
 double precision, intent(out) :: accu(N_states)
 implicit none
 integer :: i_i,i_r,i_v,i_a,i_b
 integer :: i,r,v,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 accu = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_r = 1, n_virt_orb
   r = list_virt(i_r)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,r,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,r,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do jspin=1, 2
       delta_e(i_a,jspin,istate) = one_anhil(i_a,jspin,istate)                        &
                                 - fock_virt_total_spin_trace(r,istate)               & 
                                 - fock_virt_total_spin_trace(v,istate)               & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
      enddo
     enddo
    enddo
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     do i_b = 1, n_act_orb
      b = list_act(i_b)
      do ispin = 1, 2 ! spin of (i --> r)
       do jspin = 1, 2 ! spin of (a --> v)
        if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
        do istate = 1, N_states
         if(ispin == jspin)then
          accu(istate) += (active_int(i_a,1) - active_int(i_a,2)) * one_body_dm_mo_spin_index(a,b,istate,ispin)   &
                                                          * (active_int(i_b,1) - active_int(i_b,2)) & 
                                                          * delta_e(i_a,jspin,istate)
         else 
          accu(istate) += active_int(i_a,1)  * one_body_dm_mo_spin_index(a,b,istate,ispin) * delta_e(i_a,ispin,istate) & 
                        * active_int(i_b,1) 
         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end


subroutine matrix_1h2p_dm_based(matrix)
 double precision, intent(inout) :: matrix(N_det,N_det,N_states)
 implicit none
 integer :: i_i,i_r,i_v,i_a,i_b
 integer :: i,r,v,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 double precision accu(N_states)
 accu = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_r = 1, n_virt_orb
   r = list_virt(i_r)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,r,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,r,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do jspin=1, 2
       delta_e(i_a,jspin,istate) = one_anhil(i_a,jspin,istate)                        &
                                 - fock_virt_total_spin_trace(r,istate)               & 
                                 - fock_virt_total_spin_trace(v,istate)               & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
      enddo
     enddo
    enddo
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     do i_b = 1, n_act_orb
      b = list_act(i_b)
      do ispin = 1, 2 ! spin of (i --> r)
       do jspin = 1, 2 ! spin of (a --> v)
        if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
        do istate = 1, N_states
         if(ispin == jspin)then
          accu(istate) += (active_int(i_a,1) - active_int(i_a,2)) * one_body_dm_mo_spin_index(a,b,istate,ispin)   &
                                                          * (active_int(i_b,1) - active_int(i_b,2)) & 
                                                          * delta_e(i_a,jspin,istate)
         else 
          accu(istate) += active_int(i_a,1)  * one_body_dm_mo_spin_index(a,b,istate,ispin) * delta_e(i_a,ispin,istate) & 
                        * active_int(i_b,1) 
         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
end

 BEGIN_PROVIDER [logical, list_determinants_single_exc_ok, (2,n_act_orb,n_act_orb)]
&BEGIN_PROVIDER [integer, list_determinants_single_exc, (2,2,n_act_orb,n_act_orb)]
&BEGIN_PROVIDER [double precision, phase_determinants_single_exc, (2,n_act_orb,n_act_orb)]
&BEGIN_PROVIDER [integer(bit_kind), determinants_single_exc, (N_int,2,2,2,n_act_orb,n_act_orb)]
  use bitmasks
 implicit none
 integer :: idet,jdet,aorb,borb,a,b,jspin
 double precision :: hij,phase
 integer           :: degree(N_det)
 integer           :: idx(0:N_det)
 integer :: exc(0:2,2,2)
 list_determinants_single_exc_ok = .False.
 do idet = 1, N_det
  call get_excitation_degree_vector_mono(psi_det,psi_det(1,1,idet),degree,N_int,N_det,idx)
  do jdet = 1, idx(0)
   call get_mono_excitation(psi_det(1,1,idet),psi_det(1,1,idx(jdet)),exc,phase,N_int)
   if (exc(0,1,1) == 1) then
     ! Mono alpha
     aorb  = (exc(1,2,1))   !!! a^{\dagger}_a 
     borb  = (exc(1,1,1))   !!! a_{b}
     jspin = 1
   else
     ! Mono beta
     aorb  = (exc(1,2,2))   !!!  a^{\dagger}_a
     borb  = (exc(1,1,2))   !!!  a_{b}
     jspin = 2
   endif
   a = list_act_reverse(aorb)
   b = list_act_reverse(borb)
   list_determinants_single_exc(1,jspin,a,b) = idet
   list_determinants_single_exc(2,jspin,a,b) = idx(idet)
   list_determinants_single_exc_ok(jspin,a,b) = .True.
  enddo
 enddo

END_PROVIDER 



BEGIN_PROVIDER [double precision, effective_fock_operator_1h2p, (n_act_orb,n_act_orb,2,N_states)]
 implicit none
  use bitmasks
 integer :: i_i,i_r,i_v,i_a,i_b
 integer :: i,r,v,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 effective_fock_operator_1h2p = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_r = 1, n_virt_orb
   r = list_virt(i_r)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,r,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,r,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do jspin=1, 2
       delta_e(i_a,jspin,istate) = one_anhil(i_a,jspin,istate)                        &
                                 - fock_virt_total_spin_trace(r,istate)               & 
                                 - fock_virt_total_spin_trace(v,istate)               & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
      enddo
     enddo
    enddo
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     do i_b = 1, n_act_orb
      b = list_act(i_b)
      do ispin = 1, 2 ! spin of (i --> r)
       do jspin = 1, 2 ! spin of (a --> v)
        if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
        do istate = 1, N_states
         if(ispin == jspin)then
          effective_fock_operator_1h2p(i_a,i_b,ispin,istate) += &
                          (active_int(i_a,1) - active_int(i_a,2))    &
                                                          * (active_int(i_b,1) - active_int(i_b,2)) & 
                                                          * delta_e(i_a,jspin,istate)
         else 
          effective_fock_operator_1h2p(i_a,i_b,ispin,istate) += &
                          active_int(i_a,1)  *  delta_e(i_a,ispin,istate) & 
                        * active_int(i_b,1) 
         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


 double precision :: accu
 accu = 0.d0
 do ispin = 1, 2
  do i_a = 1, n_act_orb
   a = list_act(i_a)
   do i_b = 1, n_act_orb
    b = list_act(i_b)
    accu += effective_fock_operator_1h2p(i_a,i_b,ispin,1) * one_body_dm_mo_spin_index(a,b,istate,ispin)
   enddo
  enddo
 enddo
 print*, 'accu = ',accu
END_PROVIDER 


subroutine contrib_2h1p_dm_based(accu)
 implicit none
 integer :: i_i,i_j,i_v,i_a,i_b
 integer :: i,j,v,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision, intent(out) :: accu(N_states)
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 accu = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_j = 1, n_inact_orb
   j = list_inact(i_j)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,j,v,a,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,j,a,v,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do jspin=1, 2
       delta_e(i_a,jspin,istate) = one_creat(i_a,jspin,istate)  - fock_virt_total_spin_trace(v,istate)  &
                                 + fock_core_inactive_total_spin_trace(i,istate)       &
                                 + fock_core_inactive_total_spin_trace(j,istate)        
       delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
      enddo
     enddo
    enddo
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     do i_b = 1, n_act_orb
      b = list_act(i_b)
      do ispin = 1, 2 ! spin of (i --> v)
       do jspin = 1, 2 ! spin of (j --> a)
        if(ispin == jspin .and. i.le.j)cycle ! condition not to double count 
        do istate = 1, N_states
         if(ispin == jspin)then
          accu(istate) += (active_int(i_a,1) - active_int(i_a,2)) * one_body_dm_dagger_mo_spin_index(a,b,istate,ispin)   &
                                                          * (active_int(i_b,1) - active_int(i_b,2)) & 
                                                          * delta_e(i_a,jspin,istate)
         else 
          accu(istate) += active_int(i_a,1)  * one_body_dm_dagger_mo_spin_index(a,b,istate,ispin) * delta_e(i_a,ispin,istate) & 
                        * active_int(i_b,1) 
         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end


!subroutine contrib_2p_dm_based(accu)
!implicit none
!integer :: i_r,i_v,i_a,i_b,i_c,i_d
!integer :: r,v,a,b,c,d
!integer :: ispin,jspin
!integer :: istate
!double precision, intent(out) :: accu(N_states)
!double precision :: active_int(n_act_orb,n_act_orb,2)
!double precision :: delta_e(n_act_orb,n_act_orb,2,2,N_states)
!double precision :: get_mo_bielec_integral
!accu = 0.d0
!do i_r = 1, n_virt_orb
! r = list_virt(i_r)
! do i_v = 1, n_virt_orb
!   v = list_virt(i_v)
!   do i_a = 1, n_act_orb
!    a = list_act(i_a)
!    do i_b = 1, n_act_orb
!     b = list_act(i_b)
!     active_int(i_a,i_b,1) = get_mo_bielec_integral(a,b,r,v,mo_integrals_map) ! direct
!     active_int(i_a,i_b,2) = get_mo_bielec_integral(a,b,v,r,mo_integrals_map) ! direct
!     do istate = 1, N_states
!      do jspin=1, 2 ! spin of i_a
!       do ispin = 1, 2 ! spin of i_b
!        delta_e(i_a,i_b,jspin,ispin,istate) = two_anhil(i_a,i_b,jspin,ispin,istate)    &
!                                  - fock_virt_total_spin_trace(r,istate)               & 
!                                  - fock_virt_total_spin_trace(v,istate)                 
!        delta_e(i_a,i_b,jspin,ispin,istate) = 1.d0/delta_e(i_a,i_b,jspin,ispin,istate)  
!       enddo
!      enddo
!     enddo
!    enddo
!   enddo
!   ! diagonal terms 
!   do i_a = 1, n_act_orb
!    a = list_act(i_a)
!    do i_b = 1, n_act_orb
!     b = list_act(i_b)
!     do ispin = 1, 2 ! spin of (a --> r)
!      do jspin = 1, 2 ! spin of (b --> v)
!       if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
!       if(ispin == jspin .and. a.le.b)cycle ! condition not to double count 
!       do istate = 1, N_states
!        if(ispin == jspin)then
!         double precision :: contrib_spin
!         if(ispin == 1)then
!          contrib_spin = two_body_dm_aa_diag_act(i_a,i_b)
!         else
!          contrib_spin = two_body_dm_bb_diag_act(i_a,i_b)
!         endif
!         accu(istate) += (active_int(i_a,i_b,1) - active_int(i_a,i_b,2)) * contrib_spin  &
!                       * (active_int(i_a,i_b,1) - active_int(i_a,i_b,2)) & 
!                       * delta_e(i_a,i_b,ispin,jspin,istate)
!        else 
!         accu(istate) += 0.5d0 * active_int(i_a,i_b,1)  * two_body_dm_ab_diag_act(i_a,i_b) * delta_e(i_a,i_b,ispin,jspin,istate) & 
!                               * active_int(i_a,i_b,1) 
!        endif
!       enddo
!      enddo
!     enddo
!    enddo
!   enddo
!  enddo
! enddo


!end

