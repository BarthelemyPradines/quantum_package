 BEGIN_PROVIDER [ double precision, ao_bi_elec_integral_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_bi_elec_integral_beta ,  (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Alpha Fock matrix in AO basis set
 END_DOC
 
 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer*8                      :: p,q
 double precision               :: integral, c0, c1, c2
 double precision               :: ao_bielec_integral, local_threshold
 double precision, allocatable  :: ao_bi_elec_integral_alpha_tmp(:,:)
 double precision, allocatable  :: ao_bi_elec_integral_beta_tmp(:,:)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: ao_bi_elec_integral_beta_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: ao_bi_elec_integral_alpha_tmp

 ao_bi_elec_integral_alpha = 0.d0
 ao_bi_elec_integral_beta  = 0.d0
 if (do_direct_integrals) then

   !$OMP PARALLEL DEFAULT(NONE)                                      &
       !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,p,q,r,s,i0,j0,k0,l0, &
       !$OMP ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp, c0, c1, c2, &
       !$OMP local_threshold)&
       !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
       !$OMP ao_integrals_map,ao_integrals_threshold, ao_bielec_integral_schwartz, &
       !$OMP ao_overlap_abs, ao_bi_elec_integral_alpha, ao_bi_elec_integral_beta)

   allocate(keys(1), values(1))
   allocate(ao_bi_elec_integral_alpha_tmp(ao_num,ao_num), &
            ao_bi_elec_integral_beta_tmp(ao_num,ao_num))
   ao_bi_elec_integral_alpha_tmp = 0.d0
   ao_bi_elec_integral_beta_tmp  = 0.d0

   q = ao_num*ao_num*ao_num*ao_num
   !$OMP DO SCHEDULE(dynamic)
   do p=1_8,q
           call bielec_integrals_index_reverse(kk,ii,ll,jj,p)
           if ( (kk(1)>ao_num).or. &
                (ii(1)>ao_num).or. &
                (jj(1)>ao_num).or. &
                (ll(1)>ao_num) ) then
                cycle
           endif
           k = kk(1)
           i = ii(1)
           l = ll(1)
           j = jj(1)

           if (ao_overlap_abs(k,l)*ao_overlap_abs(i,j)  &
              < ao_integrals_threshold) then
             cycle
           endif
           local_threshold = ao_bielec_integral_schwartz(k,l)*ao_bielec_integral_schwartz(i,j)
           if (local_threshold < ao_integrals_threshold) then
             cycle
           endif
           i0 = i
           j0 = j
           k0 = k
           l0 = l
           values(1) = 0.d0
           local_threshold = ao_integrals_threshold/local_threshold
           do k2=1,8
             if (kk(k2)==0) then
               cycle
             endif
             i = ii(k2)
             j = jj(k2)
             k = kk(k2)
             l = ll(k2)
             c0 = SCF_density_matrix_ao_alpha(k,l)+SCF_density_matrix_ao_beta(k,l)
             c1 = SCF_density_matrix_ao_alpha(k,i)
             c2 = SCF_density_matrix_ao_beta(k,i)
             if ( dabs(c0)+dabs(c1)+dabs(c2) < local_threshold) then
               cycle
             endif
             if (values(1) == 0.d0) then
               values(1) = ao_bielec_integral(k0,l0,i0,j0)
             endif
             integral = c0 * values(1)
             ao_bi_elec_integral_alpha_tmp(i,j) += integral
             ao_bi_elec_integral_beta_tmp (i,j) += integral
             integral = values(1)
             ao_bi_elec_integral_alpha_tmp(l,j) -= c1 * integral
             ao_bi_elec_integral_beta_tmp (l,j) -= c2 * integral
           enddo
   enddo
   !$OMP END DO NOWAIT
   !$OMP CRITICAL
   ao_bi_elec_integral_alpha += ao_bi_elec_integral_alpha_tmp
   !$OMP END CRITICAL
   !$OMP CRITICAL
   ao_bi_elec_integral_beta  += ao_bi_elec_integral_beta_tmp
   !$OMP END CRITICAL
   deallocate(keys,values,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)
   !$OMP END PARALLEL
 else
   PROVIDE ao_bielec_integrals_in_map 
   PROVIDE ao_bielec_integrals_erf_in_map 
           
   integer(omp_lock_kind) :: lck(ao_num)
   integer*8                      :: i8
   integer                        :: ii(8), jj(8), kk(8), ll(8), k2
   integer(cache_map_size_kind)   :: n_elements_max, n_elements
   integer(key_kind), allocatable :: keys(:)
   double precision, allocatable  :: values(:)
   integer(cache_map_size_kind)   :: n_elements_max_erf, n_elements_erf
   integer(key_kind), allocatable :: keys_erf(:)
   double precision, allocatable  :: values_erf(:)

   !$OMP PARALLEL DEFAULT(NONE)                                      &
       !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
       !$OMP  n_elements,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)&
       !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
       !$OMP  ao_integrals_map, ao_bi_elec_integral_alpha, ao_bi_elec_integral_beta) 

   call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
   allocate(keys(n_elements_max), values(n_elements_max))
   allocate(ao_bi_elec_integral_alpha_tmp(ao_num,ao_num), &
            ao_bi_elec_integral_beta_tmp(ao_num,ao_num))
   ao_bi_elec_integral_alpha_tmp = 0.d0
   ao_bi_elec_integral_beta_tmp  = 0.d0

   !$OMP DO SCHEDULE(dynamic,64)
   !DIR$ NOVECTOR
   do i8=0_8,ao_integrals_map%map_size
     n_elements = n_elements_max
     call get_cache_map(ao_integrals_map,i8,keys,values,n_elements)
     do k1=1,n_elements
       call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))

       do k2=1,8
         if (kk(k2)==0) then
           cycle
         endif
         i = ii(k2)
         j = jj(k2)
         k = kk(k2)
         l = ll(k2)
         integral = (SCF_density_matrix_ao_alpha(k,l)+SCF_density_matrix_ao_beta(k,l)) * values(k1)
         ao_bi_elec_integral_alpha_tmp(i,j) += integral
         ao_bi_elec_integral_beta_tmp (i,j) += integral
       enddo
     enddo
   enddo
   !$OMP END DO NOWAIT
   !$OMP CRITICAL
   ao_bi_elec_integral_alpha += ao_bi_elec_integral_alpha_tmp
   !$OMP END CRITICAL
   !$OMP CRITICAL
   ao_bi_elec_integral_beta  += ao_bi_elec_integral_beta_tmp
   !$OMP END CRITICAL
   deallocate(keys,values,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)
   !$OMP END PARALLEL

   !$OMP PARALLEL DEFAULT(NONE)                                      &
       !$OMP PRIVATE(i,j,l,k1,k,integral_erf,ii,jj,kk,ll,i8,keys_erf,values_erf,n_elements_max_erf, &
       !$OMP  n_elements_erf,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)&
       !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
       !$OMP  ao_integrals_erf_map, ao_bi_elec_integral_alpha, ao_bi_elec_integral_beta) 


   call get_cache_map_n_elements_max(ao_integrals_erf_map,n_elements_max_erf)
   allocate(ao_bi_elec_integral_alpha_tmp(ao_num,ao_num), &
            ao_bi_elec_integral_beta_tmp(ao_num,ao_num))
   allocate(keys_Erf(n_elements_max_erf), values_erf(n_elements_max_erf))

   ao_bi_elec_integral_alpha_tmp = 0.d0
   ao_bi_elec_integral_beta_tmp  = 0.d0
   !$OMP DO SCHEDULE(dynamic,64)
   !DIR$ NOVECTOR
   do i8=0_8,ao_integrals_erf_map%map_size
     n_elements_erf = n_elements_max_erf
     call get_cache_map(ao_integrals_erf_map,i8,keys_erf,values_erf,n_elements_erf)
     do k1=1,n_elements_erf
       call bielec_integrals_index_reverse(kk,ii,ll,jj,keys_erf(k1))

       do k2=1,8
         if (kk(k2)==0) then
           cycle
         endif
         i = ii(k2)
         j = jj(k2)
         k = kk(k2)
         l = ll(k2)
         double precision :: integral_erf
         integral_erf = values_erf(k1)
         ao_bi_elec_integral_alpha_tmp(l,j) -= (SCF_density_matrix_ao_alpha(k,i) * integral_erf)
         ao_bi_elec_integral_beta_tmp (l,j) -= (SCF_density_matrix_ao_beta (k,i) * integral_erf)
       enddo
     enddo
   enddo

   !$OMP END DO NOWAIT
   !$OMP CRITICAL
   ao_bi_elec_integral_alpha += ao_bi_elec_integral_alpha_tmp
   !$OMP END CRITICAL
   !$OMP CRITICAL
   ao_bi_elec_integral_beta  += ao_bi_elec_integral_beta_tmp
   !$OMP END CRITICAL
   deallocate(ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)
   deallocate(keys_erf,values_erf)
   !$OMP END PARALLEL

 endif

END_PROVIDER

 
 BEGIN_PROVIDER [ double precision, Fock_matrix_ao_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_ao_beta,  (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Alpha Fock matrix in AO basis set
 END_DOC
 
 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     Fock_matrix_ao_alpha(i,j) = Fock_matrix_alpha_no_xc_ao(i,j) + ao_potential_alpha_xc(i,j)
     Fock_matrix_ao_beta (i,j) = Fock_matrix_beta_no_xc_ao(i,j)  + ao_potential_beta_xc(i,j)
   enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [ double precision, Fock_matrix_alpha_no_xc_ao, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_beta_no_xc_ao,  (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Mono electronic an Coulomb matrix in AO basis set
 END_DOC
 
 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     Fock_matrix_alpha_no_xc_ao(i,j) = ao_mono_elec_integral(i,j) + ao_bi_elec_integral_alpha(i,j) 
     Fock_matrix_beta_no_xc_ao(i,j) = ao_mono_elec_integral(i,j) + ao_bi_elec_integral_beta (i,j) 
   enddo
 enddo

END_PROVIDER



 BEGIN_PROVIDER [ double precision, RS_KS_energy ]
!BEGIN_PROVIDER [ double precision, SCF_energy ]
&BEGIN_PROVIDER [ double precision, two_electron_energy]
&BEGIN_PROVIDER [ double precision, one_electron_energy]
&BEGIN_PROVIDER [ double precision, Fock_matrix_energy]
&BEGIN_PROVIDER [ double precision, trace_potential_xc ]
 implicit none
 BEGIN_DOC
 ! Range-separated Kohn-Sham energy
 END_DOC
 RS_KS_energy = nuclear_repulsion
 
 integer                        :: i,j
 double precision :: accu_mono,accu_fock
 one_electron_energy = 0.d0
 two_electron_energy = 0.d0
 Fock_matrix_energy = 0.d0
 trace_potential_xc = 0.d0
 do j=1,ao_num
   do i=1,ao_num
    Fock_matrix_energy +=   Fock_matrix_ao_alpha(i,j) * SCF_density_matrix_ao_alpha(i,j) + & 
                            Fock_matrix_ao_beta(i,j) * SCF_density_matrix_ao_beta(i,j) 
    two_electron_energy += 0.5d0 * ( ao_bi_elec_integral_alpha(i,j) * SCF_density_matrix_ao_alpha(i,j) & 
                +ao_bi_elec_integral_beta(i,j) * SCF_density_matrix_ao_beta(i,j) ) 
    one_electron_energy += ao_mono_elec_integral(i,j) * (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
! possible bug fix for open-shell
!    trace_potential_xc += (ao_potential_alpha_xc(i,j) + ao_potential_beta_xc(i,j) ) *  (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
    trace_potential_xc += ao_potential_alpha_xc(i,j) * SCF_density_matrix_ao_alpha(i,j) + ao_potential_beta_xc(i,j) *  SCF_density_matrix_ao_beta (i,j)
   enddo
 enddo
 RS_KS_energy +=  e_exchange_dft + e_correlation_dft + one_electron_energy + two_electron_energy
!SCF_energy = RS_KS_energy 
END_PROVIDER 

BEGIN_PROVIDER [double precision, extra_energy_contrib_from_density]
 implicit none
! possible bug fix for open-shell:
! extra_energy_contrib_from_density = e_exchange_dft + e_correlation_dft - 0.25d0 * trace_potential_xc
 extra_energy_contrib_from_density = e_exchange_dft + e_correlation_dft - 0.5d0 * trace_potential_xc
END_PROVIDER 

!BEGIN_PROVIDER [ double precision, SCF_energy ]
! implicit none
! SCF_energy = RS_KS_energy 
!END_PROVIDER 
