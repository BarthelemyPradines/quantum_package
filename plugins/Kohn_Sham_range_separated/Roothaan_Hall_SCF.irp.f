subroutine Roothaan_Hall_SCF

BEGIN_DOC
! Roothaan-Hall algorithm for SCF Hartree-Fock calculation
END_DOC

  implicit none

  double precision               :: energy_SCF,energy_SCF_previous,Delta_energy_SCF
  double precision               :: max_error_DIIS,max_error_DIIS_alpha,max_error_DIIS_beta
  double precision, allocatable  :: Fock_matrix_DIIS(:,:,:),error_matrix_DIIS(:,:,:)

  integer                        :: iteration_SCF,dim_DIIS,index_dim_DIIS
 
  integer                        :: i,j
  double precision, allocatable :: mo_coef_save(:,:)
 
  allocate(mo_coef_save(ao_num,mo_tot_num),                          &
      Fock_matrix_DIIS (ao_num,ao_num,max_dim_DIIS),                 &
      error_matrix_DIIS(ao_num,ao_num,max_dim_DIIS)                  &
      )

  call write_time(output_kohn_sham)

  write(output_kohn_sham,'(A4, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================'
  write(output_kohn_sham,'(A4, 1X, A16, 1X, A16, 1X, A16)')  &
    '  N ', 'Energy  ', 'Energy diff  ', 'DIIS error  '
  write(output_kohn_sham,'(A4, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================'

! Initialize energies and density matrices

  energy_SCF_previous = RS_KS_energy
  Delta_energy_SCF    = 1.d0
  iteration_SCF       = 0
  dim_DIIS            = 0
  max_error_DIIS      = 1.d0


!
! Start of main SCF loop
!
  do while(( (max_error_DIIS > threshold_DIIS_nonzero).or.(dabs(Delta_energy_SCF) > thresh_SCF) ) .and. (iteration_SCF < n_it_SCF_max))

! Increment cycle number

    iteration_SCF += 1

! Current size of the DIIS space

    dim_DIIS = min(dim_DIIS+1,max_dim_DIIS)
 
    if (scf_algorithm == 'DIIS') then
      
      ! Store Fock and error matrices at each iteration
      do j=1,ao_num
        do i=1,ao_num
          index_dim_DIIS = mod(dim_DIIS-1,max_dim_DIIS)+1
          Fock_matrix_DIIS (i,j,index_dim_DIIS) = Fock_matrix_ao(i,j)
          error_matrix_DIIS(i,j,index_dim_DIIS) = FPS_SPF_matrix_ao(i,j)
        enddo
      enddo
      
      ! Compute the extrapolated Fock matrix

      call extrapolate_Fock_matrix(                                    &
          error_matrix_DIIS,Fock_matrix_DIIS,                          &
          Fock_matrix_ao,size(Fock_matrix_ao,1),                       &
          iteration_SCF,dim_DIIS                                       &
          )

      Fock_matrix_ao_alpha = Fock_matrix_ao*0.5d0
      Fock_matrix_ao_beta  = Fock_matrix_ao*0.5d0
      TOUCH Fock_matrix_ao_alpha Fock_matrix_ao_beta

    endif

    mo_coef = eigenvectors_Fock_matrix_mo

    TOUCH mo_coef

   !Calculate error vectors

    max_error_DIIS = maxval(Abs(FPS_SPF_Matrix_mo))

   !SCF energy

    energy_SCF = RS_KS_energy
    print*,'RS_KS_energy main =',RS_KS_energy
    Delta_Energy_SCF = energy_SCF - energy_SCF_previous
    if ( (SCF_algorithm == 'DIIS').and.(Delta_Energy_SCF > 0.d0) ) then
      Fock_matrix_ao(1:ao_num,1:ao_num) = Fock_matrix_DIIS (1:ao_num,1:ao_num,index_dim_DIIS)
      Fock_matrix_ao_alpha = Fock_matrix_ao*0.5d0
      Fock_matrix_ao_beta  = Fock_matrix_ao*0.5d0
      TOUCH Fock_matrix_ao_alpha Fock_matrix_ao_beta
    endif

    double precision :: level_shift_save
    level_shift_save = level_shift
    mo_coef_save(1:ao_num,1:mo_tot_num) = mo_coef(1:ao_num,1:mo_tot_num)
    do while (Delta_Energy_SCF > 0.d0)
      mo_coef(1:ao_num,1:mo_tot_num) = mo_coef_save
      TOUCH mo_coef
      level_shift = level_shift + 0.1d0
      mo_coef(1:ao_num,1:mo_tot_num) = eigenvectors_Fock_matrix_mo(1:ao_num,1:mo_tot_num)
      TOUCH mo_coef level_shift
      Delta_Energy_SCF = RS_KS_energy - energy_SCF_previous
      print*,'RS_KS_energy = ',RS_KS_energy,Delta_Energy_SCF
      energy_SCF = RS_KS_energy
      if (level_shift-level_shift_save > 1.d0) exit
      dim_DIIS=0
    enddo
    level_shift = level_shift_save
    SOFT_TOUCH level_shift
    energy_SCF_previous = energy_SCF

   !Print results at the end of each iteration

    write(output_kohn_sham,'(I4, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, I3)')  &
      iteration_SCF, energy_SCF, Delta_energy_SCF, max_error_DIIS, dim_DIIS

    if (Delta_energy_SCF < 0.d0) then
      call save_mos
    endif

  enddo

 
 !End of Main SCF loop
 

  write(output_kohn_sham,'(A4, 1X, A16, 1X, A16, 1X, A16)') &
    '====','================','================','================'
  write(output_kohn_sham,*)
  
  if(.not.no_oa_or_av_opt)then
   call mo_as_eigvectors_of_mo_matrix(Fock_matrix_mo,size(Fock_matrix_mo,1),size(Fock_matrix_mo,2),mo_label,1,.true.)
  endif

  call write_double(output_kohn_sham, Energy_SCF, 'Hartree-Fock energy')
  call ezfio_set_hartree_fock_energy(Energy_SCF)

  call write_time(output_kohn_sham)

end


subroutine extrapolate_Fock_matrix(      &
  error_matrix_DIIS,Fock_matrix_DIIS,    &
  Fock_matrix_ao_,size_Fock_matrix_ao,   &
  iteration_SCF,dim_DIIS                 &
)

BEGIN_DOC
! Compute the extrapolated Fock matrix using the DIIS procedure
END_DOC

  implicit none

  double precision,intent(in)   :: Fock_matrix_DIIS(ao_num,ao_num,*),error_matrix_DIIS(ao_num,ao_num,*)
  integer,intent(in)            :: iteration_SCF, size_Fock_matrix_ao
  double precision,intent(inout):: Fock_matrix_ao_(size_Fock_matrix_ao,ao_num)
  integer,intent(inout)         :: dim_DIIS

  double precision,allocatable  :: B_matrix_DIIS(:,:),X_vector_DIIS(:)
  double precision,allocatable  :: C_vector_DIIS(:)

  double precision,allocatable  :: scratch(:,:)
  integer                       :: i,j,k,i_DIIS,j_DIIS

  allocate(                               &
    B_matrix_DIIS(dim_DIIS+1,dim_DIIS+1), &
    X_vector_DIIS(dim_DIIS+1),            &
    C_vector_DIIS(dim_DIIS+1),            &
    scratch(ao_num,ao_num)                &
  )

! Compute the matrices B and X
  do j=1,dim_DIIS
    do i=1,dim_DIIS

      j_DIIS = mod(iteration_SCF-j,max_dim_DIIS)+1
      i_DIIS = mod(iteration_SCF-i,max_dim_DIIS)+1

! Compute product of two errors vectors

      call dgemm('N','N',ao_num,ao_num,ao_num,                      &
           1.d0,                                                    &
           error_matrix_DIIS(1,1,i_DIIS),size(error_matrix_DIIS,1), &
           error_matrix_DIIS(1,1,j_DIIS),size(error_matrix_DIIS,1), &
           0.d0,                                                    &
           scratch,size(scratch,1))

! Compute Trace

      B_matrix_DIIS(i,j) = 0.d0
      do k=1,ao_num
        B_matrix_DIIS(i,j) = B_matrix_DIIS(i,j) + scratch(k,k)
      enddo
    enddo
  enddo

! Pad B matrix and build the X matrix

  do i=1,dim_DIIS
    B_matrix_DIIS(i,dim_DIIS+1) = -1.d0
    B_matrix_DIIS(dim_DIIS+1,i) = -1.d0
    C_vector_DIIS(i) = 0.d0
  enddo
  B_matrix_DIIS(dim_DIIS+1,dim_DIIS+1) = 0.d0
  C_vector_DIIS(dim_DIIS+1) = -1.d0

! Solve the linear system C = B.X

  integer              :: info
  integer,allocatable  :: ipiv(:)

  allocate(          &
    ipiv(dim_DIIS+1) &
  )

  double precision, allocatable :: AF(:,:)
  allocate (AF(dim_DIIS+1,dim_DIIS+1))
  double precision :: rcond, ferr, berr
  integer :: iwork(dim_DIIS+1), lwork

  call dsysvx('N','U',dim_DIIS+1,1,      &
    B_matrix_DIIS,size(B_matrix_DIIS,1), &
    AF, size(AF,1),                      &
    ipiv,                                &
    C_vector_DIIS,size(C_vector_DIIS,1), &
    X_vector_DIIS,size(X_vector_DIIS,1), &
    rcond,                               &
    ferr,                                &
    berr,                                &
    scratch,-1,                          &
    iwork,                               &
    info                                 &
  )
  lwork = int(scratch(1,1))
  deallocate(scratch)
  allocate(scratch(lwork,1))

  call dsysvx('N','U',dim_DIIS+1,1,      &
    B_matrix_DIIS,size(B_matrix_DIIS,1), &
    AF, size(AF,1),                      &
    ipiv,                                &
    C_vector_DIIS,size(C_vector_DIIS,1), &
    X_vector_DIIS,size(X_vector_DIIS,1), &
    rcond,                               &
    ferr,                                &
    berr,                                &
    scratch,size(scratch),               &
    iwork,                               &
    info                                 &
  )

 if(info < 0) then
   stop 'bug in DIIS'
 endif

 if (rcond > 1.d-12) then

  ! Compute extrapolated Fock matrix


      !$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(SHARED)
      do j=1,ao_num
        do i=1,ao_num
          Fock_matrix_ao_(i,j) = 0.d0
        enddo
        do k=1,dim_DIIS
          do i=1,ao_num
            Fock_matrix_ao_(i,j) = Fock_matrix_ao_(i,j) +            &
                X_vector_DIIS(k)*Fock_matrix_DIIS(i,j,dim_DIIS-k+1)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

  else
    dim_DIIS = 0
  endif

end

