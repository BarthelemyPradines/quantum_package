 character*20 function str(k)
!   "Convert an integer to string."
 implicit none
 integer, intent(in) :: k
 write (str, *) k
 str = adjustl(str)
 end 


double precision function binom_func(i,j)
  implicit none
  BEGIN_DOC
  !.. math                       :: 
  !
  !  \frac{i!}{j!(i-j)!}
  !
  END_DOC
  integer,intent(in)             :: i,j
  double precision               :: logfact
  integer, save                  :: ifirst
  double precision, save         :: memo(0:15,0:15)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: memo
  integer                        :: k,l
  if (ifirst == 0) then
    ifirst = 1
    do k=0,15
      do l=0,15
        memo(k,l) = dexp( logfact(k)-logfact(l)-logfact(k-l) )
      enddo
    enddo
  endif
  if ( (i<=15).and.(j<=15) ) then
    binom_func = memo(i,j)
  else
    binom_func = dexp( logfact(i)-logfact(j)-logfact(i-j) )
  endif
end

double precision function inverse_erf(y,thr)
 implicit none
 BEGIN_DOC
! y  = erf(inverser_erf) +- thr
 END_DOC
 double precision, intent(in) :: thr, y
 double precision :: xmin,xmax,xminbefore
 xmin = 0.d0
 xmax = 6.d0
!do while (dabs((xmin)-(xmax)).gt.thr)
 do while (dabs(erf(xmin) - y).gt.thr)
  if(erf(xmin).le.y)then
   xminbefore = xmin
   xmin += (xmax-xmin) * 0.5d0
  else if(erf(xmin).gt.y)then
   xmax = xmin
   xmin = xminbefore + (xmin-xminbefore) * 0.5d0
  endif
 enddo
 inverse_erf = xmin

end


 BEGIN_PROVIDER [ double precision, binom, (0:40,0:40) ]
&BEGIN_PROVIDER [ double precision, binom_transp, (0:40,0:40) ]
  implicit none
  BEGIN_DOC
  ! Binomial coefficients
  END_DOC
  integer                        :: k,l
  double precision               :: logfact
  do k=0,40
    do l=0,40
      binom(k,l) = dexp( logfact(k)-logfact(l)-logfact(k-l) )
      binom_transp(l,k) = binom(k,l)
    enddo
  enddo
END_PROVIDER



double precision function fact(n)
  implicit none
  BEGIN_DOC
  ! n!
  END_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  
  if (n<=memomax) then
    if (n<2) then
      fact = 1.d0
    else
      fact = memo(n)
    endif
    return
  endif
  
  integer                        :: i
  memo(1) = 1.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)*dble(i)
  enddo
  memomax = min(n,100)
  double precision :: logfact
  fact = dexp(logfact(n))
end function

double precision function logfact(n)
  implicit none
  BEGIN_DOC
  ! n!
  END_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  
  if (n<=memomax) then
    if (n<2) then
      logfact = 0.d0
    else
      logfact = memo(n)
    endif
    return
  endif
  
  integer                        :: i
  memo(1) = 0.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)+dlog(dble(i))
  enddo
  memomax = min(n,100)
  logfact = memo(memomax)
  do i=101,n
    logfact += dlog(dble(i))
  enddo
end function



BEGIN_PROVIDER [ double precision, fact_inv, (128) ]
  implicit none
  BEGIN_DOC
  ! 1/n!
  END_DOC
  integer                        :: i
  double precision               :: fact
  do i=1,size(fact_inv)
    fact_inv(i) = 1.d0/fact(i)
  enddo
END_PROVIDER


double precision function dble_fact(n)
  implicit none
  integer :: n
  double precision :: dble_fact_even, dble_fact_odd

  dble_fact = 1.d0

  if(n.lt.0) return

  if(iand(n,1).eq.0)then
    dble_fact = dble_fact_even(n)
  else
    dble_fact= dble_fact_odd(n)
  endif

end function

double precision function dble_fact_even(n) result(fact2)
  implicit none
  BEGIN_DOC
  ! n!!
  END_DOC
  integer                        :: n,k
  double precision, save         :: memo(0:100)
  integer, save                  :: memomax = 0
  double precision               :: prod

  ASSERT (iand(n,1) /= 1)

!  prod=1.d0
!  do k=2,n,2
!   prod=prod*dfloat(k)
!  enddo
!  fact2=prod
!  return
!
  if (n <= memomax) then
    if (n < 2) then
      fact2 = 1.d0
    else
      fact2 = memo(n)
    endif
    return
  endif

  integer                        :: i
  memo(0)=1.d0
  memo(1)=1.d0
  do i=memomax+2,min(n,100),2
    memo(i) = memo(i-2)* dble(i)
  enddo
  memomax = min(n,100)
  fact2 = memo(memomax)
  
  if (n > 100) then
    double precision :: dble_logfact
    fact2 = dexp(dble_logfact(n))
  endif

end function

double precision function dble_fact_odd(n) result(fact2)
  implicit none
  BEGIN_DOC
  ! n!!
  END_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  
  ASSERT (iand(n,1) /= 0)
  if (n<=memomax) then
    if (n<3) then
      fact2 = 1.d0
    else
      fact2 = memo(n)
    endif
    return
  endif
  
  integer                        :: i
  memo(1) = 1.d0
  do i=memomax+2,min(n,99),2
    memo(i) = memo(i-2)* dble(i)
  enddo
  memomax = min(n,99)
  fact2 = memo(memomax)
  
  if (n > 99) then
    double precision :: dble_logfact
    fact2 = dexp(dble_logfact(n))
  endif

end function

double precision function dble_logfact(n) result(logfact2)
  implicit none
  BEGIN_DOC
  ! n!!
  END_DOC
  integer                        :: n
  integer :: k
  double precision :: prod
  prod=0.d0
  do k=2,n,2
   prod=prod+dlog(dfloat(k))
  enddo
  logfact2=prod
  return
  
end function

subroutine write_git_log(iunit)
  implicit none
  BEGIN_DOC
  ! Write the last git commit in file iunit.
  END_DOC
  integer, intent(in)            :: iunit
  write(iunit,*) '----------------'
  write(iunit,*) 'Last git commit:'
  BEGIN_SHELL [ /bin/bash ]
  git log -1 2>/dev/null | sed "s/'//g"| sed "s/^/    write(iunit,*) '/g" | sed "s/$/'/g" || echo "Unknown"
  END_SHELL
  write(iunit,*) '----------------'
end

BEGIN_PROVIDER [ double precision, inv_int, (128) ]
  implicit none
  BEGIN_DOC
  ! 1/i
  END_DOC
  integer                        :: i
  do i=1,128
    inv_int(i) = 1.d0/dble(i)
  enddo
END_PROVIDER

subroutine wall_time(t)
  implicit none
  BEGIN_DOC
  ! The equivalent of cpu_time, but for the wall time.
  END_DOC
  double precision, intent(out)  :: t
  integer                        :: c
  integer, save                  :: rate = 0
  if (rate == 0) then
    CALL SYSTEM_CLOCK(count_rate=rate)
  endif
  CALL SYSTEM_CLOCK(count=c)
  t = dble(c)/dble(rate)
end

BEGIN_PROVIDER [ integer, nproc ]
  implicit none
  BEGIN_DOC
  ! Number of current OpenMP threads
  END_DOC
  
  integer                        :: omp_get_num_threads
  nproc = 1
  !$OMP PARALLEL
  !$OMP MASTER
  !$ nproc = omp_get_num_threads()
  !$OMP END MASTER
  !$OMP END PARALLEL
END_PROVIDER


double precision function u_dot_v(u,v,sze)
  implicit none
  BEGIN_DOC
  ! Compute <u|v>
  END_DOC
  integer, intent(in)            :: sze
  double precision, intent(in)   :: u(sze),v(sze)
  double precision, external     :: ddot
  
  !DIR$ FORCEINLINE
  u_dot_v = ddot(sze,u,1,v,1)
  
end

double precision function u_dot_u(u,sze)
  implicit none
  BEGIN_DOC
  ! Compute <u|u>
  END_DOC
  integer, intent(in)            :: sze
  double precision, intent(in)   :: u(sze)
  double precision, external     :: ddot
  
  !DIR$ FORCEINLINE
  u_dot_u = ddot(sze,u,1,u,1)
  
end

subroutine normalize(u,sze)
  implicit none
  BEGIN_DOC
  ! Normalizes vector u
  END_DOC
  integer, intent(in)            :: sze
  double precision, intent(inout):: u(sze)
  double precision               :: d
  double precision, external     :: dnrm2
  integer                        :: i
  
  !DIR$ FORCEINLINE
  d = dnrm2(sze,u,1)
  if (d /= 0.d0) then
    d = 1.d0/d
  endif
  if (d /= 1.d0) then
    !DIR$ FORCEINLINE
    call dscal(sze,d,u,1)
  endif
end

double precision function approx_dble(a,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: a
  double precision :: f
  integer :: i

  if (a == 0.d0) then
    approx_dble = 0.d0
    return
  endif
  f = 1.d0
  do i=1,-int(dlog10(dabs(a)))+n
    f = f*.1d0
  enddo
  do i=1,int(dlog10(dabs(a)))-n
    f = f*10.d0
  enddo
  approx_dble = dnint(a/f)*f

end



subroutine lowercase(txt,n)
  implicit none
  BEGIN_DOC
! Transform to lower case
  END_DOC
  character*(*), intent(inout)   :: txt
  integer, intent(in)            :: n
  character( * ), PARAMETER      :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character( * ), PARAMETER      :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer                        :: i, ic
  do i=1,n
    ic = index( UPPER_CASE, txt(i:i) )
    if (ic /= 0) then
      txt(i:i) = LOWER_CASE(ic:ic)
    endif
  enddo
end


 double precision function power(n,x)
 implicit none
 integer :: i,n
 double precision :: x
 power = 1 
 do i = 1, n
  power *= x
 enddo
 return
 end

