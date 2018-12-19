 BEGIN_PROVIDER [double precision, big_array_coulomb_integrals, (mo_tot_num,mo_tot_num, mo_tot_num)]
&BEGIN_PROVIDER [double precision, big_array_exchange_integrals,(mo_tot_num,mo_tot_num, mo_tot_num)]
&BEGIN_PROVIDER [double precision, big_array_integrals_anti,(mo_tot_num,mo_tot_num, mo_tot_num)]
 implicit none
 BEGIN_DOC
! big_array_coulomb_integrals (j,i,k) = <ij|kj> = (ik|jj)
! big_array_exchange_integrals(j,i,k) = <ij|jk> = (ij|kj)
! big_array_integrals_anti(j,i,k)     = <ij|kj> - <ij|jk> = (ik|jj) - (ij|kj)
 END_DOC
 integer :: i,j,k,l
 double precision :: get_mo_bielec_integral
 double precision :: integral_direct,integral_exch

 do k = 1, mo_tot_num
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
     l = j
     integral_direct = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
     big_array_coulomb_integrals(j,i,k) = integral_direct
     l = j
     integral_exch = get_mo_bielec_integral(i,j,l,k,mo_integrals_map)
     big_array_exchange_integrals(j,i,k) = integral_exch
     big_array_integrals_anti(j,i,k) = integral_direct - integral_exch
   enddo
  enddo
 enddo


END_PROVIDER 
