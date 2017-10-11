subroutine set_intermediate_normalization_lmct_old(norm,i_hole)
 implicit none
 integer, intent(in) :: i_hole
 double precision, intent(out) :: norm(N_states)
 integer :: i,j,degree,index_ref_generators_restart(N_states),k,l
 integer::  number_of_holes,n_h, number_of_particles,n_p
 integer, allocatable :: index_one_hole(:),index_one_hole_one_p(:),index_two_hole_one_p(:),index_two_hole(:)
 integer, allocatable :: index_one_p(:)
 integer :: n_one_hole,n_one_hole_one_p,n_two_hole_one_p,n_two_hole,n_one_p
 logical :: is_the_hole_in_det
 double precision :: inv_coef_ref_generators_restart(N_states),hij,hii,accu
 integer :: index_good_hole(1000)
 integer :: n_good_hole
 logical,allocatable :: is_a_ref_det(:)
 allocate(index_one_hole(n_det),index_one_hole_one_p(n_det),index_two_hole_one_p(N_det),index_two_hole(N_det),index_one_p(N_det),is_a_ref_det(N_det))
 
 n_one_hole = 0
 n_one_hole_one_p = 0
 n_two_hole_one_p = 0
 n_two_hole = 0
 n_one_p = 0
 n_good_hole = 0
 ! Find the one holes and one hole one particle
 is_a_ref_det = .False.


 double precision :: ovlp(N_states)
 integer :: i_good_state(N_states),is_chosen(N_states)
 integer :: n_good_state,iorder(N_States)
 double precision :: overlap_ref(N_states,N_states),norm_overlap(N_states)
 integer, allocatable :: index_ref_gen(:)
 allocate(index_ref_gen(N_det_generators_restart))

 overlap_ref = 0.d0
 norm_overlap = 0.d0

 do i = 1, N_det
   ! Find the reference determinant for intermediate normalization
  do k = 1, N_states
   call get_excitation_degree(ref_generators_restart(1,1,k),psi_det(1,1,i),degree,N_int)   
   if(degree == 0)then
    index_ref_generators_restart(k) = i
    inv_coef_ref_generators_restart(k) = 1.d0/psi_coef(i,k)
   endif
  enddo
  
  ! Find all the determinants present in the reference wave function
  do j = 1, N_det_generators_restart
   call get_excitation_degree(psi_det(1,1,i),psi_det_generators_restart(1,1,j),degree,N_int)  
   if(degree == 0)then
    is_a_ref_det(i) = .True.
    index_ref_gen(j) = i  
    do l = 1, N_states
     norm_overlap(l) += psi_coef(i,l) **2 
     do k = 1, N_states
      overlap_ref(k,l) += psi_coef(i,l) *  psi_coef_generators_restart(j,k)
     enddo
    enddo
    exit
   endif
  enddo
  if(is_a_ref_det(i))cycle
  n_h = number_of_holes(psi_det(1,1,i))
  n_p = number_of_particles(psi_det(1,1,i))
  if(n_h == 1 .and. n_p == 0)then
   if(is_the_hole_in_det(psi_det(1,1,i),1,i_hole).or.is_the_hole_in_det(psi_det(1,1,i),2,i_hole))then  
    n_good_hole +=1
    index_good_hole(n_good_hole) = i
!   print*, 'hole',psi_coef(i,1),psi_coef(i,2)
!   print*, psi_coef(index_ref_generators_restart(1),1),psi_coef(index_ref_generators_restart(2),2)
   else
    do k = 1, N_states
     psi_coef(i,k) = 0.d0
    enddo
   endif
  else
   do k = 1, N_states
    psi_coef(i,k) = 0.d0
   enddo
  endif
 enddo

 n_good_state = 0
 is_chosen = .False.
 do l = 1, N_states
  print*, 'state ',l
  ovlp = 0.d0
  iorder = 0
  do k = 1, N_states
   if(is_chosen(k))cycle
   iorder(k) = k
   overlap_ref(k,l) = overlap_ref(k,l)  / dsqrt(norm_overlap(l)) * 1.d0 / dsqrt(norm_generators_restart(k))
   ovlp(k) = -dabs(overlap_ref(k,l))
  enddo
  call dsort(ovlp,iorder,N_states)
  if(dabs(ovlp(1)).gt.0.75d0)then
   n_good_state +=1 
   i_good_state(l) = iorder(1)
  endif
  print*, 'ovlp(1)',ovlp(1)
  print*, 'iorder(1)',iorder(1)
 enddo




 if(n_good_state .ne. N_states)then
  print*, 'WARNING !!!! '
  print*, 'THE STATES FOUND ARE NOT THE ONE YOU REQUESTED !!!'
  print*, 'SKIPPING THIS OF THE FOBOCI ....'
  do k = 1, N_states
   do j = 1, N_det
    psi_coef(j,k) = 0.d0
   enddo
   do j = 1, N_det_generators_restart
    psi_coef(index_ref_gen(j),k) = psi_coef_generators_restart(j,k)
   enddo
  enddo
 else 
  logical :: need_reordering
  need_reordering = .False.
  do k = 1, N_states
   if(i_good_state(k).ne.k)then
    need_reordering = .True.
   endif
  enddo
  double precision, allocatable :: psi_coef_tmp(:,:)
  if(need_reordering)then
    print*, 'REORDERING THE STATES ..'
    allocate(psi_coef_tmp(N_det, N_states))
    do i = 1, N_states
     do j = 1, N_det
      psi_coef_tmp(j,i) = psi_coef(j,i)
     enddo
    enddo
    do i = 1, N_States
     print*, 'Statae = ',i
     do j = 1, N_det
      psi_coef(j,i) = psi_coef_tmp(j,i_good_state(i))
      ! Find the reference determinant for intermediate normalization
      call get_excitation_degree(ref_generators_restart(1,1,i),psi_det(1,1,j),degree,N_int)   
      if(degree == 0)then
       call debug_det(psi_det(1,1,j),N_int)
       print*, psi_coef(j,i)
       index_ref_generators_restart(i) = i
       inv_coef_ref_generators_restart(i) = 1.d0/psi_coef(j,i)
      endif
     enddo
    enddo
  endif
 endif
 print*,''
 print*,'n_good_hole = ',n_good_hole
 do k = 1,N_states
  print*,'state ',k
  do i = 1, n_good_hole
   print*,'psi_coef(index_good_hole) = ',psi_coef(index_good_hole(i),k) * inv_coef_ref_generators_restart(k)
  enddo
  print*,''
 enddo

 ! Set the wave function to the intermediate normalization
 do k = 1, N_states
  do i = 1, N_det
   psi_coef(i,k) = psi_coef(i,k) * inv_coef_ref_generators_restart(k)
  enddo
 enddo
 double precision :: norm_ref(N_States)
 norm_ref = 0.d0
 norm = 0.d0
 do k = 1,N_states
  print*,'state ',k
  do i = 1, N_det
   if (is_a_ref_det(i))then
    print*,'psi_coef_ref = ',psi_coef(i,k)
!   call debug_det(psi_det(1,1,i),N_int)
    norm_ref(k) += psi_coef(i,k) * psi_coef(i,k)
    cycle
   endif
   norm(k) += psi_coef(i,k) * psi_coef(i,k)
  enddo
  print*,'norm = ',norm(k)
  norm_ref(k) = 1.d0/dsqrt(norm_ref(k))
 enddo
 deallocate(index_one_hole,index_one_hole_one_p,index_two_hole_one_p,index_two_hole,index_one_p,is_a_ref_det)
 soft_touch psi_coef
end


subroutine set_intermediate_normalization_mlct_old(norm,i_particl)
 implicit none
 integer, intent(in) :: i_particl
 double precision, intent(out) :: norm(N_states)
 integer :: i,j,degree,index_ref_generators_restart(N_states),k,l
 integer::  number_of_holes,n_h, number_of_particles,n_p
 integer, allocatable :: index_one_hole(:),index_one_hole_one_p(:),index_two_hole_one_p(:),index_two_hole(:)
 integer, allocatable :: index_one_p(:),index_one_hole_two_p(:)
 integer :: n_one_hole,n_one_hole_one_p,n_two_hole_one_p,n_two_hole,n_one_p,n_one_hole_two_p
 logical :: is_the_particl_in_det
 double precision :: inv_coef_ref_generators_restart(N_states)
 integer          :: exc(0:2,2,2)
 double precision :: phase,hij,hii,accu
 integer :: h1,p1,h2,p2,s1,s2
 integer :: index_good_particl(1000)
 integer :: n_good_particl
 logical,allocatable :: is_a_ref_det(:)
 integer :: i_count
 allocate(index_one_hole(n_det),index_one_hole_one_p(n_det),index_two_hole_one_p(N_det),index_two_hole(N_det),index_one_p(N_det),is_a_ref_det(N_det))
 allocate(index_one_hole_two_p(n_det))

 double precision :: ovlp(N_states)
 integer :: i_good_state(N_states),is_chosen(N_states)
 integer :: n_good_state,iorder(N_States)
 double precision :: overlap_ref(N_states,N_states),norm_overlap(N_states)
 integer, allocatable :: index_ref_gen(:)
 allocate(index_ref_gen(N_det_generators_restart))

 overlap_ref = 0.d0
 norm_overlap = 0.d0

 
 n_one_hole = 0
 n_one_hole_one_p = 0
 n_two_hole_one_p = 0
 n_two_hole = 0
 n_one_p = 0
 n_one_hole_two_p = 0
 n_good_particl = 0
 ! Find the one holes and one hole one particle
 i_count = 0
 is_a_ref_det = .False.
 do i = 1, N_det
  do k =1, N_states
   call get_excitation_degree(ref_generators_restart(1,1,k),psi_det(1,1,i),degree,N_int)
   if(degree == 0)then
    index_ref_generators_restart(k) = i
    inv_coef_ref_generators_restart(k) = 1.d0/psi_coef(i,k)
   endif
  enddo

  ! Find all the determinants present in the reference wave function
  do j = 1, N_det_generators_restart
   call get_excitation_degree(psi_det(1,1,i),psi_det_generators_restart(1,1,j),degree,N_int)  
   if(degree == 0)then
    is_a_ref_det(i) = .True.
    index_ref_gen(j) = i  
    do l = 1, N_states
     norm_overlap(l) += psi_coef(i,l) **2 
     do k = 1, N_states
      overlap_ref(k,l) += psi_coef(i,l) *  psi_coef_generators_restart(j,k)
     enddo
    enddo
    exit
   endif
  enddo
  if(is_a_ref_det(i))cycle


  n_h = number_of_holes(psi_det(1,1,i))
  n_p = number_of_particles(psi_det(1,1,i))
  if(n_h == 0 .and. n_p == 1)then   ! 1p
   if(is_the_particl_in_det(psi_det(1,1,i),1,i_particl).or.is_the_particl_in_det(psi_det(1,1,i),2,i_particl))then  
    n_good_particl += 1
    index_good_particl(n_good_particl) = i
   else
    do k = 1, N_states
     psi_coef(i,k) = 0.d0
    enddo
   endif
  else
   do k = 1, N_states
    psi_coef(i,k) = 0.d0
   enddo
  endif
 enddo

 n_good_state = 0
 is_chosen = .False.
 do l = 1, N_states
  print*, 'state ',l
  ovlp = 0.d0
  iorder = 0
  do k = 1, N_states
   if(is_chosen(k))cycle
   iorder(k) = k
   overlap_ref(k,l) = overlap_ref(k,l)  / dsqrt(norm_overlap(l)) * 1.d0 / dsqrt(norm_generators_restart(k))
   ovlp(k) = -dabs(overlap_ref(k,l))
  enddo
  call dsort(ovlp,iorder,N_states)
  if(dabs(ovlp(1)).gt.0.75d0)then
   n_good_state +=1 
   i_good_state(l) = iorder(1)
  endif
  print*, 'ovlp(1)',ovlp(1)
  print*, 'iorder(1)',iorder(1)
 enddo




 if(n_good_state .ne. N_states)then
  print*, 'WARNING !!!! '
  print*, 'THE STATES FOUND ARE NOT THE ONE YOU REQUESTED !!!'
  print*, 'SKIPPING THIS OF THE FOBOCI ....'
  do k = 1, N_states
   do j = 1, N_det
    psi_coef(j,k) = 0.d0
   enddo
   do j = 1, N_det_generators_restart
    psi_coef(index_ref_gen(j),k) = psi_coef_generators_restart(j,k)
   enddo
  enddo
 else 
  logical :: need_reordering
  need_reordering = .False.
  do k = 1, N_states
   if(i_good_state(k).ne.k)then
    need_reordering = .True.
   endif
  enddo
  double precision, allocatable :: psi_coef_tmp(:,:)
  if(need_reordering)then
    print*, 'REORDERING THE STATES ..'
    allocate(psi_coef_tmp(N_det, N_states))
    do i = 1, N_states
     do j = 1, N_det
      psi_coef_tmp(j,i) = psi_coef(j,i)
     enddo
    enddo
    do i = 1, N_States
     print*, 'Statae = ',i
     do j = 1, N_det
      psi_coef(j,i) = psi_coef_tmp(j,i_good_state(i))
      ! Find the reference determinant for intermediate normalization
      call get_excitation_degree(ref_generators_restart(1,1,i),psi_det(1,1,j),degree,N_int)   
      if(degree == 0)then
       call debug_det(psi_det(1,1,j),N_int)
       print*, psi_coef(j,i)
       index_ref_generators_restart(i) = i
       inv_coef_ref_generators_restart(i) = 1.d0/psi_coef(j,i)
      endif
     enddo
    enddo
  endif
 endif

 norm = 0.d0
 print*,''
 print*,'n_good_particl = ',n_good_particl
 do k = 1, N_states
   print*,'state ',k
   do i = 1, n_good_particl
    print*,'psi_coef(index_good_particl,1) = ',psi_coef(index_good_particl(i),k) * inv_coef_ref_generators_restart(k)
   enddo
   print*,''
 enddo


 double precision :: norm_ref(N_States)
 norm_ref = 0.d0
 do k = 1, N_states
  print*,'state ',k
  do i = 1, N_det
   if (is_a_ref_det(i))then
    print*,'i,psi_coef_ref = ',psi_coef(i,k)
    norm_ref(k) += psi_coef(i,k) * psi_coef(i,k)
    cycle
   endif
   norm(k) += psi_coef(i,k) * psi_coef(i,k)
  enddo
  print*,'norm = ',norm
  norm_ref(k) = 1.d0/dsqrt(norm_ref(k))
 enddo
 ! Set the wave function to the intermediate normalization
 do k = 1, N_states
  do i = 1, N_det
   psi_coef(i,k) = psi_coef(i,k) * inv_coef_ref_generators_restart(k)
  enddo
 enddo

!! Set the wave function to the intermediate normalization
!do k = 1, N_states
! do i = 1, N_det
!  psi_coef(i,k) = psi_coef(i,k) * norm_ref(k)
! enddo
!enddo


 soft_touch psi_coef
 deallocate(index_one_hole,index_one_hole_one_p,index_two_hole_one_p,index_two_hole,index_one_p,is_a_ref_det)
end


subroutine update_density_matrix_osoci
 implicit none
 BEGIN_DOC
 ! one_body_dm_mo_alpha_osoci += Delta rho alpha
 ! one_body_dm_mo_beta_osoci  += Delta rho beta
 END_DOC
 integer :: i,j,k
 do k = 1, N_states
  do i = 1, mo_tot_num
!  print*, one_body_dm_mo_alpha_osoci(i,i,k) ,  one_body_dm_mo_alpha_old(i,i,k) , one_body_dm_mo_alpha_generators_restart(i,i,k) 
   do j = 1, mo_tot_num
    one_body_dm_mo_alpha_osoci(i,j,k) = one_body_dm_mo_alpha_osoci(i,j,k) + (one_body_dm_mo_alpha_old(i,j,k) - one_body_dm_mo_alpha_generators_restart(i,j,k))
    one_body_dm_mo_beta_osoci(i,j,k) = one_body_dm_mo_beta_osoci(i,j,k) + (one_body_dm_mo_beta_old(i,j,k) - one_body_dm_mo_beta_generators_restart(i,j,k))
   enddo
  enddo
 enddo


end

subroutine initialize_density_matrix_osoci
 implicit none
 one_body_dm_mo_alpha_osoci = one_body_dm_mo_alpha_generators_restart
 one_body_dm_mo_beta_osoci  = one_body_dm_mo_beta_generators_restart
end

subroutine rescale_density_matrix_osoci(norm)
 implicit none
 double precision, intent(in) :: norm(N_states)
 integer :: i,j,k
 do k = 1, N_states
  do i = 1, mo_tot_num
   do j = 1,mo_tot_num
    one_body_dm_mo_alpha_osoci(i,j,k) = one_body_dm_mo_alpha_osoci(i,j,k) * norm(k)
    one_body_dm_mo_beta_osoci(j,i,k) = one_body_dm_mo_beta_osoci(j,i,k) * norm(k)
   enddo
  enddo
 enddo
end

subroutine save_osoci_natural_mos(norm_total)

 implicit none
 BEGIN_DOC
 ! Set natural orbitals, obtained by diagonalization of the one-body density matrix in the MO basis
 END_DOC
 character*(64) :: label
 double  precision, intent(in) :: norm_total(N_states)
 double precision, allocatable :: tmp(:,:),tmp_bis(:,:)
 integer, allocatable :: occ(:,:)
 integer       :: n_occ_alpha,i,i_core,j_core,iorb,jorb,j,i_inact,j_inact,i_virt,j_virt,k
 allocate(tmp(size(one_body_dm_mo_alpha_osoci,1),size(one_body_dm_mo_alpha_osoci,2)))
 allocate(tmp_bis(size(one_body_dm_mo_alpha_osoci,1),size(one_body_dm_mo_alpha_osoci,2)))
 allocate (occ(N_int*bit_kind_size,2))

 ! Negation to have the occupied MOs first after the diagonalization

 print*, 'Intermediate normalization DM '
 do k = 1, N_states
  print*,'State ',k
  print*,''
  print*,'Inactive-active Part of the One body DM'
  print*,''
  do i = 1,n_act_orb
   iorb = list_act(i)
   print*,''
   print*,'ACTIVE ORBITAL  ',iorb
   do j = 1, n_inact_orb
    jorb = list_inact(j)
    if(dabs(one_body_dm_mo_alpha_osoci(iorb,jorb,k) + one_body_dm_mo_beta_osoci(iorb,jorb,k)).gt.0.0001d0)then
!   if(jorb == 59)then
     print*,'INACTIVE  '
     print*,'DM ',iorb,jorb,one_body_dm_mo_alpha_osoci(iorb,jorb,k) + one_body_dm_mo_beta_osoci(iorb,jorb,k)
     print*,'DM ',iorb,iorb,one_body_dm_mo_alpha_osoci(iorb,iorb,k) + one_body_dm_mo_beta_osoci(iorb,iorb,k)
     print*,'DM ',jorb,jorb,one_body_dm_mo_alpha_osoci(jorb,jorb,k) + one_body_dm_mo_beta_osoci(jorb,jorb,k)
    endif
   enddo
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    if(dabs(one_body_dm_mo_alpha_osoci(iorb,jorb,k) + one_body_dm_mo_beta_osoci(iorb,jorb,k)).gt.0.0001d0)then
     print*,'VIRT      '
     print*,'DM ',iorb,jorb,one_body_dm_mo_alpha_osoci(iorb,jorb,k) + one_body_dm_mo_beta_osoci(iorb,jorb,k)
     print*,'DM ',iorb,iorb,one_body_dm_mo_alpha_osoci(iorb,iorb,k) + one_body_dm_mo_beta_osoci(iorb,iorb,k)
     print*,'DM ',jorb,jorb,one_body_dm_mo_alpha_osoci(jorb,jorb,k) + one_body_dm_mo_beta_osoci(jorb,jorb,k)
    endif
   enddo
  enddo
 enddo
 
 
  tmp_bis = 0.d0
  print*, 'norm_total(k)',norm_total(:)
  do j = 1, mo_tot_num
   do i = 1, mo_tot_num
    do k = 1, N_states
     tmp_bis(i,j) -=  (one_body_dm_mo_alpha_osoci(i,j,k) * 1.d0/norm_total(k) * state_average_weight(k)+ one_body_dm_mo_beta_osoci(i,j,k) * 1.d0/norm_total(k)* state_average_weight(k))  
    enddo
   enddo
  enddo
 ! Set to Zero the core-inact-act-virt part
 do i = 1, n_core_orb
  i_core = list_core(i)
  tmp_bis(i_core,i_core) = -10.d0
  do j = i+1, n_core_orb
   j_core = list_core(j)
   tmp_bis(i_core,j_core) = 0.d0
   tmp_bis(j_core,i_core) = 0.d0
  enddo
  do j = 1, n_inact_orb
   iorb = list_inact(j)
   tmp_bis(i_core,iorb) = 0.d0
   tmp_bis(iorb,i_core) = 0.d0
  enddo
  do j = 1, n_act_orb
   iorb = list_act(j)
   tmp_bis(i_core,iorb) = 0.d0
   tmp_bis(iorb,i_core) = 0.d0
  enddo
  do j = 1, n_virt_orb
   iorb = list_virt(j)
   tmp_bis(i_core,iorb) = 0.d0
   tmp_bis(iorb,i_core) = 0.d0
  enddo
 enddo
 ! Set to Zero the inact-inact part to avoid arbitrary rotations
 do i = 1, n_inact_orb
  i_inact = list_inact(i)
  do j = i+1, n_inact_orb 
   j_inact = list_inact(j)
   tmp_bis(i_inact,j_inact) = 0.d0
   tmp_bis(j_inact,i_inact) = 0.d0
  enddo
 enddo

 ! Set to Zero the inact-virt part to avoid arbitrary rotations
 do i = 1, n_inact_orb
  i_inact = list_inact(i)
  do j = 1, n_virt_orb 
   j_virt = list_virt(j)
   tmp_bis(i_inact,j_virt) = 0.d0
   tmp_bis(j_virt,i_inact) = 0.d0
  enddo
 enddo

 ! Set to Zero the virt-virt part to avoid arbitrary rotations
 do i = 1, n_virt_orb
  i_virt = list_virt(i)
  do j = i+1, n_virt_orb 
   j_virt = list_virt(j)
   tmp_bis(i_virt,j_virt) = 0.d0
   tmp_bis(j_virt,i_virt) = 0.d0
  enddo
 enddo

 double precision :: accu
 ! Set to Zero the act-act part to avoid arbitrary rotations
 do i = 1,n_act_orb
  iorb = list_act(i)
  do j = i+1,n_act_orb
   jorb = list_act(j)
   tmp_bis(iorb,jorb) = 0.d0
   tmp_bis(jorb,iorb) = 0.d0
  enddo
 enddo

 tmp = tmp_bis

 call bitstring_to_list(reunion_of_bitmask(1,1), occ(1,1), n_occ_alpha, N_int)
 double precision :: maxvaldm,imax,jmax
 maxvaldm = 0.d0
 imax = 1
 jmax = 1
 print*, 'Aproximatly DM '
 print*,''
 print*,'Inactive-active Part of the One body DM'
 print*,''
 do i = 1,n_act_orb
  iorb = list_act(i)
  print*,''
  print*,'ACTIVE ORBITAL  ',iorb
  do j = 1, n_inact_orb
   jorb = list_inact(j)
   if(dabs(tmp(iorb,jorb)).gt.0.0001d0)then
    print*,'INACTIVE  '
    print*,'DM ',iorb,jorb,(tmp(iorb,jorb))
    print*,'DM ',jorb,jorb,(tmp(jorb,jorb))
    print*,'DM ',iorb,iorb,(tmp(iorb,iorb))
   endif
  enddo
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   if(dabs(tmp(iorb,jorb)).gt.0.0001d0)then
    print*,'VIRT      '
    print*,'DM ',iorb,jorb,(tmp(iorb,jorb))
    print*,'DM ',jorb,jorb,(tmp(jorb,jorb))
    print*,'DM ',iorb,iorb,(tmp(iorb,iorb))
   endif
  enddo
 enddo
 do i = 1, mo_tot_num
  do j = i+1, mo_tot_num
   if(dabs(tmp(i,j)).le.threshold_fobo_dm)then
      tmp(i,j) = 0.d0
      tmp(j,i) = 0.d0
   endif
  enddo
 enddo

 label = "Natural"
 
 call mo_as_eigvectors_of_mo_matrix(tmp,size(tmp,1),size(tmp,2),label,1)
 touch mo_coef psi_det psi_coef
 deallocate(tmp,occ)


end

subroutine check_symetry(i_hole,thr,test)
 implicit none
 integer, intent(in) :: i_hole
 double precision, intent(in) :: thr
 logical, intent(out) :: test
 integer :: i,j,k,l
 double precision :: accu
 accu = 0.d0
 do i = 1, n_act_orb
  accu += dabs(mo_mono_elec_integral(i_hole,list_act(i)))
 enddo
 if(accu.gt.thr)then
  test = .True.
 else
  test = .false.
 endif
end

subroutine check_symetry_1h1p(i_hole,i_part,thr,test)
 implicit none
 integer, intent(in) :: i_hole,i_part
 double precision, intent(in) :: thr
 logical, intent(out) :: test
 integer :: i,j,k,l
 double precision :: accu
 accu = dabs(mo_mono_elec_integral(i_hole,i_part))
 if(accu.gt.thr)then
  test = .True.
 else
  test = .false.
 endif
end


 subroutine update_one_body_dm_mo(norm_total)
   implicit none
   double precision, intent(in) :: norm_total(N_states)
   integer :: i,k
   double precision :: accu_tot,accu_sd
   print*,'touched the one_body_dm_mo_beta'
   one_body_dm_mo_alpha_average = 0.d0
   one_body_dm_mo_beta_average  = 0.d0
   do k = 1, N_states
    one_body_dm_mo_alpha_average += one_body_dm_mo_alpha_osoci(:,:,k) * state_average_weight(k) /norm_total(k)
    one_body_dm_mo_beta_average  += one_body_dm_mo_beta_osoci(:,:,k)  * state_average_weight(k) /norm_total(k)
   enddo
   touch one_body_dm_mo_alpha  one_body_dm_mo_beta 
   accu_tot = 0.d0
   accu_sd  = 0.d0
   do i = 1, mo_tot_num
    accu_tot += one_body_dm_mo_alpha_average(i,i) + one_body_dm_mo_beta_average(i,i)
    accu_sd  += one_body_dm_mo_alpha_average(i,i) - one_body_dm_mo_beta_average(i,i)
   enddo
   print*,'accu_tot = ',accu_tot
   print*,'accu_sdt = ',accu_sd 
 end
 
 subroutine provide_properties
   implicit none
   call print_mulliken_sd
   call print_hcc
 end



 subroutine dress_diag_elem_2h1p(dressing_H_mat_elem,ndet,lmct,i_hole)
  use bitmasks
  double precision, intent(inout) :: dressing_H_mat_elem(Ndet)
  integer, intent(in) :: ndet,i_hole
  logical, intent(in) :: lmct
  ! if lmct = .True. ===> LMCT
  ! else             ===> MLCT
  implicit none
  integer :: i
  integer :: n_p,n_h,number_of_holes,number_of_particles
  integer :: exc(0:2,2,2)
  integer :: degree
  double precision :: phase
  integer :: h1,h2,p1,p2,s1,s2
  do i = 1, N_det

   n_h = number_of_holes(psi_det(1,1,i))
   n_p = number_of_particles(psi_det(1,1,i))
   call get_excitation(ref_bitmask,psi_det(1,1,i),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   if (n_h == 0.and.n_p==0)then ! CAS
    dressing_H_mat_elem(i)+= total_corr_e_2h1p
    if(lmct)then
     dressing_H_mat_elem(i) += - corr_energy_2h1p_per_orb_ab(i_hole) - corr_energy_2h1p_per_orb_bb(i_hole)    
    endif
   endif
   if (n_h == 1.and.n_p==0)then ! 1h
     dressing_H_mat_elem(i)+= 0.d0
   else if (n_h == 0.and.n_p==1)then ! 1p
    dressing_H_mat_elem(i)+= total_corr_e_2h1p
    dressing_H_mat_elem(i) += - corr_energy_2h1p_per_orb_ab(p1) - corr_energy_2h1p_per_orb_aa(p1)  
   else if (n_h == 1.and.n_p==1)then ! 1h1p
!   if(degree==1)then
     dressing_H_mat_elem(i)+= total_corr_e_2h1p
     dressing_H_mat_elem(i)+= - corr_energy_2h1p_per_orb_ab(h1) 
!   else 
!    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(h1) &  
!    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(h1) + corr_energy_2h2p_per_orb_bb(h1))
!    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(p2) &  
!    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(p2) + corr_energy_2h2p_per_orb_bb(p2))
!    dressing_H_mat_elem(i) += 0.5d0 * (corr_energy_2h2p_for_1h1p_double(h1,p1))
!   endif
   else if (n_h == 2.and.n_p==1)then ! 2h1p
    dressing_H_mat_elem(i)+= 0.d0
   else if (n_h == 1.and.n_p==2)then ! 1h2p
    dressing_H_mat_elem(i)+= total_corr_e_2h1p
    dressing_H_mat_elem(i) += - corr_energy_2h1p_per_orb_ab(h1) 
   endif
  enddo
 
 end

 subroutine dress_diag_elem_1h2p(dressing_H_mat_elem,ndet,lmct,i_hole)
  use bitmasks
  double precision, intent(inout) :: dressing_H_mat_elem(Ndet)
  integer, intent(in) :: ndet,i_hole
  logical, intent(in) :: lmct
  ! if lmct = .True. ===> LMCT
  ! else             ===> MLCT
  implicit none
  integer :: i
  integer :: n_p,n_h,number_of_holes,number_of_particles
  integer :: exc(0:2,2,2)
  integer :: degree
  double precision :: phase
  integer :: h1,h2,p1,p2,s1,s2
  do i = 1, N_det

   n_h = number_of_holes(psi_det(1,1,i))
   n_p = number_of_particles(psi_det(1,1,i))
   call get_excitation(ref_bitmask,psi_det(1,1,i),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   if (n_h == 0.and.n_p==0)then ! CAS
    dressing_H_mat_elem(i)+= total_corr_e_1h2p
    if(.not.lmct)then
     dressing_H_mat_elem(i) += - corr_energy_1h2p_per_orb_ab(i_hole) - corr_energy_1h2p_per_orb_aa(i_hole)    
    endif
   endif
   if (n_h == 1.and.n_p==0)then ! 1h
     dressing_H_mat_elem(i)+= total_corr_e_1h2p - corr_energy_1h2p_per_orb_ab(h1)
   else if (n_h == 0.and.n_p==1)then ! 1p
    dressing_H_mat_elem(i)+= 0.d0
   else if (n_h == 1.and.n_p==1)then ! 1h1p
    if(degree==1)then
     dressing_H_mat_elem(i)+= total_corr_e_1h2p
     dressing_H_mat_elem(i)+= - corr_energy_1h2p_per_orb_ab(h1) 
    else 
     dressing_H_mat_elem(i) +=0.d0
    endif
!    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(h1) &  
!    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(h1) + corr_energy_2h2p_per_orb_bb(h1))
!    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(p2) &  
!    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(p2) + corr_energy_2h2p_per_orb_bb(p2))
!    dressing_H_mat_elem(i) += 0.5d0 * (corr_energy_2h2p_for_1h1p_double(h1,p1))
!   endif
   else if (n_h == 2.and.n_p==1)then ! 2h1p
    dressing_H_mat_elem(i)+= total_corr_e_1h2p
    dressing_H_mat_elem(i)+= - corr_energy_1h2p_per_orb_ab(h1) - corr_energy_1h2p_per_orb_ab(h1)
   else if (n_h == 1.and.n_p==2)then ! 1h2p
    dressing_H_mat_elem(i) += 0.d0
   endif
  enddo
 
 end

 subroutine dress_diag_elem_2h2p(dressing_H_mat_elem,ndet)
  use bitmasks
  double precision, intent(inout) :: dressing_H_mat_elem(Ndet)
  integer, intent(in) :: ndet
  implicit none
  integer :: i
  integer :: n_p,n_h,number_of_holes,number_of_particles
  integer :: exc(0:2,2,2)
  integer :: degree
  double precision :: phase
  integer :: h1,h2,p1,p2,s1,s2
  do i = 1, N_det
   dressing_H_mat_elem(i)+= total_corr_e_2h2p

   n_h = number_of_holes(psi_det(1,1,i))
   n_p = number_of_particles(psi_det(1,1,i))
   call get_excitation(ref_bitmask,psi_det(1,1,i),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   if (n_h == 1.and.n_p==0)then ! 1h
    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(h1) &  
    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(h1) + corr_energy_2h2p_per_orb_bb(h1))
   else if (n_h == 0.and.n_p==1)then ! 1p
    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(p1) &  
    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(p1) + corr_energy_2h2p_per_orb_bb(p1))
   else if (n_h == 1.and.n_p==1)then ! 1h1p
    if(degree==1)then
     dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(h1) &  
     - 0.5d0 * (corr_energy_2h2p_per_orb_aa(h1) + corr_energy_2h2p_per_orb_bb(h1))
     dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(p1) &  
     - 0.5d0 * (corr_energy_2h2p_per_orb_aa(p1) + corr_energy_2h2p_per_orb_bb(p1))
     dressing_H_mat_elem(i) += 0.5d0 * (corr_energy_2h2p_for_1h1p_a(h1,p1) + corr_energy_2h2p_for_1h1p_b(h1,p1))
    else 
     dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(h1) &  
     - 0.5d0 * (corr_energy_2h2p_per_orb_aa(h1) + corr_energy_2h2p_per_orb_bb(h1))
     dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(p2) &  
     - 0.5d0 * (corr_energy_2h2p_per_orb_aa(p2) + corr_energy_2h2p_per_orb_bb(p2))
     dressing_H_mat_elem(i) += 0.5d0 * (corr_energy_2h2p_for_1h1p_double(h1,p1))
    endif
   else if (n_h == 2.and.n_p==1)then ! 2h1p
    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(h1) - corr_energy_2h2p_per_orb_bb(h1) & 
                              - corr_energy_2h2p_per_orb_ab(h2)                                   & 
                              - 0.5d0 * ( corr_energy_2h2p_per_orb_bb(h2) + corr_energy_2h2p_per_orb_bb(h2)) 
    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(p1)
    if(s1.ne.s2)then
      dressing_H_mat_elem(i) +=  corr_energy_2h2p_ab_2_orb(h1,h2)
    else 
      dressing_H_mat_elem(i) +=  corr_energy_2h2p_bb_2_orb(h1,h2)
    endif
   else if (n_h == 1.and.n_p==2)then ! 1h2p
    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(h1) &  
    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(h1) + corr_energy_2h2p_per_orb_bb(h1))
    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(p1) &  
    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(p1) + corr_energy_2h2p_per_orb_bb(p1))
    dressing_H_mat_elem(i) += - corr_energy_2h2p_per_orb_ab(p2) &  
    - 0.5d0 * (corr_energy_2h2p_per_orb_aa(p2) + corr_energy_2h2p_per_orb_bb(p2))
    if(s1.ne.s2)then
      dressing_H_mat_elem(i) +=  corr_energy_2h2p_ab_2_orb(p1,p2)
    else 
      dressing_H_mat_elem(i) +=  corr_energy_2h2p_bb_2_orb(p1,p2)
    endif
   endif
  enddo
 
 end

 subroutine diag_dressed_2h2p_hamiltonian_and_update_psi_det(i_hole,lmct)
 implicit none
 double precision, allocatable :: dressing_H_mat_elem(:),energies(:)
  integer, intent(in) :: i_hole
  logical, intent(in) :: lmct
  ! if lmct = .True. ===> LMCT
  ! else             ===> MLCT
 integer :: i
 double precision :: hij
 allocate(dressing_H_mat_elem(N_det),energies(N_states_diag))
  print*,''
  print*,'dressing with the 2h2p in a CC logic'
  print*,''
  do i = 1, N_det
   call i_h_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hij)
   dressing_H_mat_elem(i) = hij
  enddo
  call dress_diag_elem_2h2p(dressing_H_mat_elem,N_det)
  call dress_diag_elem_2h1p(dressing_H_mat_elem,N_det,lmct,i_hole)
  call dress_diag_elem_1h2p(dressing_H_mat_elem,N_det,lmct,i_hole)
  call davidson_diag_hjj(psi_det,psi_coef,dressing_H_mat_elem,energies,size(psi_coef,1),N_det,N_states,N_states_diag,N_int,output_determinants)
  do i = 1, 2
   print*,'psi_coef = ',psi_coef(i,1)
  enddo
  
  
 deallocate(dressing_H_mat_elem)



 end
