program pouet
 implicit none
 read_wf = .True.
 touch read_wf
!call test_superci_eigenvectors
!call test_natural_orbitals 
!call test_ref_energy
!call test_guess_superci 
 call test_fock 
!call test_H_superci
!call test_fock_from_act
end 

subroutine test_superci_eigenvectors
 implicit none
 print*, eigenvalues_sci_brutal(1),eigenvalues_sci(1)
end

subroutine test_fock
 implicit none
 integer :: i,j,iorb,jorb
  print*, 'Reference energy = ',reference_energy_superci
! provide diagonal_superci_matrix
  print*, 'eigenvalues',eigenvalues_sci(1),eigenvalues_sci_brutal(1)
! provide diagonal_superci_matrix 
 do i = 1, size_super_ci
!do i = 1, 1
  ! iorb = index_rotation_CI_reverse(i,1)
  ! jorb = index_rotation_CI_reverse(i,2)
  ! print*, superci_matrix(i,i)-total_matrix_superci_brutal(i,i),transformed_occ1_virt1_occ2_virt2(iorb,jorb,iorb,jorb)
  do j = 1, size_super_ci 
!  if(dabs(superci_matrix(i,j)-total_matrix_superci_brutal(i,j)).gt.1.d-10.and.i.ne.1)then
   if(dabs(superci_matrix(i,j)-total_matrix_superci_brutal(i,j)).gt.1.d-10)then
    print*, ' ***** ',i,j
    print*, index_rotation_CI_reverse(i,1),index_rotation_CI_reverse(i,2)
    print*, index_rotation_CI_reverse(j,1),index_rotation_CI_reverse(j,2)
    print*, superci_matrix(i,j),total_matrix_superci_brutal(i,j),dabs(superci_matrix(i,j)-total_matrix_superci_brutal(i,j))
    print*, '***'
   endif
  enddo
 enddo

 

end 


subroutine test_H_superci
 implicit none
 double precision, allocatable :: u1(:),u0(:)
 allocate(u1(size_super_ci),u0(size_super_ci))
 double precision, allocatable :: v1(:),v0(:)
 allocate(v1(size_super_ci),v0(size_super_ci))
 integer :: i,j
 u0 = 0.d0
 u1 = 0.d0
 u0(1) = 1.d0

 v0 = 0.d0
 v1 = 0.d0
 v0(1) = 1.d0

 do i = 1, size_super_ci
  u1(i) += u0(i) * superci_matrix(i,i) 
  do j = i+1, size_super_ci
   u1(i) += u0(j) * superci_matrix(j,i) 
   u1(j) += u0(i) * superci_matrix(j,i)
  enddo
 enddo
 

 call apply_H_superci_to_vector(v0,v1)
 
 do i = 1, size_super_ci
  if(dabs(v1(i) - u1(i)).gt.1.d-12)then
   print*, 'PB !!! '
   print*, i,v1(i),u1(i)
  endif
 enddo

 print*, 'passed the first one '
 pause

 u0 = u1
 u1 = 0.d0
 v0 = v1
 v1 = 0.d0

 do i = 1, size_super_ci
  u1(i) += u0(i) * superci_matrix(i,i) 
  do j = i+1, size_super_ci
   u1(i) += u0(j) * superci_matrix(j,i) 
   u1(j) += u0(i) * superci_matrix(j,i)
  enddo
 enddo
 

 call apply_H_superci_to_vector(v0,v1)
 
 do i = 1, size_super_ci
  if(dabs(v1(i) - u1(i)).gt.1.d-12)then
   print*, 'PB !!! '
   print*, i,v1(i),u1(i)
  endif
 enddo

 deallocate(u1,u0,v0,v1)
end


subroutine test_guess_superci
 implicit none
 provide eigenvalues_sci 
end

subroutine test_ref_energy
 implicit none
 print*,reference_energy_superci,psi_energy(1)+nuclear_repulsion
end 

subroutine test_natural_orbitals
 implicit none
 character*(64) :: label
 label = "Natural"
 integer :: i
 
 print*, ''
 print*, ''
 print*, 'SUPER CI DM'
 print*, ''
 do i = 1, mo_tot_num
  write(*,'(100(F10.7,X))')super_ci_density_matrix_mo(i,:)
 enddo
 call mo_as_svd_vectors_of_mo_matrix(super_ci_density_matrix_mo,size(super_ci_density_matrix_mo,1),mo_tot_num,mo_tot_num,label)


end

subroutine test_fock_from_act
 implicit none
 integer :: i,j,iorb,jorb

 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   if(dabs(Fock_matrix_alpha_from_act_mo_bis(jorb,iorb) - Fock_matrix_beta_from_act_mo(iorb,jorb)).gt.1.d-10)then
    print*, iorb,jorb
    print*, Fock_matrix_alpha_from_act_mo_bis(jorb,iorb),Fock_matrix_beta_from_act_mo(iorb,jorb),dabs(Fock_matrix_alpha_from_act_mo_bis(jorb,iorb) - Fock_matrix_beta_from_act_mo(iorb,jorb))
   endif
  enddo
 enddo



end
