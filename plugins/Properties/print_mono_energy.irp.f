program print_mono_energy
 implicit none
 integer :: i,j
  print*,'i  list(i)   <phi_i | h_core | phi_i> '
 do i = 1, n_act_orb
  write(*,'(2(I3,X),F10.5)')i,list_act(i),mo_mono_elec_integral(list_act(i),list_act(i))
 enddo
 

end
