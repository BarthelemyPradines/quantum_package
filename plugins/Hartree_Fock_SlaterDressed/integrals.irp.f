!*****************************************************************************
subroutine GauSlaOverlap(expGau,cGau,aGau,expSla,cSla,result)
  implicit none

  BEGIN_DOC
  ! Compute the overlap integral between a Gaussian function
  ! with arbitrary angular momemtum and a s-type Slater function
  END_DOC

! Input variables 
  double precision,intent(in)   :: expGau,expSla
  double precision,intent(in)   :: cGau(3),cSla(3)
  integer,intent(in)            :: aGau(3)
  double precision,intent(out)  :: result

! Final value of the integrals
  double precision              :: ss,ps,ds
  double precision              :: pxs,pys,pzs
  double precision              :: dxxs,dyys,dzzs,dxys,dxzs,dyzs

  double precision              :: pi,E,AB,AxBx,AyBy,AzBz,t,u,k

  pi = 4d0*atan(1d0)
  
! calculate the length AB between the two centers and other usful quantities

  AB = (cGau(1)-cSla(1))**2d0 + (cGau(2)-cSla(2))**2d0 + (cGau(3)-cSla(3))**2d0
  AB = sqrt(AB)

  AxBx = (cGau(1)-cSla(1))/2d0
  AyBy = (cGau(2)-cSla(2))/2d0
  AzBz = (cGau(3)-cSla(3))/2d0

! intermediate variables

  t = expSla*sqrt(0.25d0/expGau)
  u = sqrt(expGau)*AB

  if(AB > 0d0) then 

!   (s|s) 
    ss = (t+u)*erfc(t+u)*exp(2d0*t*(t+u)) - (t-u)*erfc(t-u)*exp(2d0*t*(t-u))

!   (p|s)
    ps = (exp(t**2d0-u**2d0)*(-4d0*t+sqrt(pi)*(exp((t-u)**2d0)*(1d0+2d0*t*(t-u))*erfc(t-u) & 
       + exp((t+u)**2d0)*(1d0+2d0*t*(t+u))*erfc(t+u))))/sqrt(pi)

!   (d|s)
    ds = 4d0*exp(2d0*t*(t-u))*t*(-((1d0+t**2d0-t*u)*erfc(t-u))+exp(4d0*t*u)*(1d0+t*(t+u))*erfc(t+u))

!   backward scaling
    ds = 3d0*ss/u**5d0 - 3d0*ps/u**4d0 + ds/u**3d0
    ps = ps/u**2d0-ss/u**3d0
    ss = ss/u

  else

!   concentric case
    ss = 2d0*exp(t**2d0)*((-2d0*t)/sqrt(pi)+exp(t**2d0)*(1d0+2d0*t**2d0)*erfc(t))
    ps = (8d0*exp(t**2d0)*t*(-2d0*(1d0+t**2d0)+exp(t**2d0)*sqrt(pi)*t*(3d0+2d0*t**2d0)*erfc(t)))/(3d0*sqrt(pi))

  endif

  k = t**3d0*exp(-t**2d0)*4d0*pi/expSla**(3d0/2d0)

! (s|s) 
  ss = k*ss

! (p|s) 
  ps = k*ps

  pxs = AxBx*ps
  pys = AyBy*ps
  pzs = AzBz*ps

! (d|s) 
  ds = k*ds

  dxxs = (2d0*ss+ps)/(4d0*expGau) + AxBx**2d0*ds
  dyys = (2d0*ss+ps)/(4d0*expGau) + AyBy**2d0*ds
  dzzs = (2d0*ss+ps)/(4d0*expGau) + AzBz**2d0*ds

  dxys = AxBx*AyBy*ds
  dxzs = AxBx*AzBz*ds
  dyzs = AyBy*AzBz*ds

  select case (sum(aGau))
    case (0)
      result = ss

    case (1)
      if (aGau(1) == 1) then
        result = pxs
      else if (aGau(2) == 1) then
        result = pys
      else if (aGau(3) == 1) then
        result = pzs
      endif

    case (2)
      if (aGau(1) == 2) then
        result = dxxs
      else if (aGau(2) == 2) then
        result = dyys
      else if (aGau(3) == 2) then
        result = dzzs
      else if (aGau(1)+aGau(2) == 2) then
        result = dxys
      else if (aGau(1)+aGau(3) == 2) then
        result = dxzs
      else if (aGau(2)+aGau(3) == 2) then
        result = dyzs
      endif

    case default
      stop 'GauSlaOverlap not implemented'

  end select

end
!*****************************************************************************


!*****************************************************************************
subroutine GauSlaKinetic(expGau,cGau,aGau,expSla,cSla)

  implicit none

  BEGIN_DOC
  ! Compute the kinetic energy integral between a Gaussian function
  ! with arbitrary angular momemtum and a s-type Slater function
  END_DOC

! Input variables 
  double precision,intent(in)   :: expGau,expSla
  double precision,intent(in)   :: cGau(3),cSla(3)
  integer,intent(in)            :: aGau(3)

! Final value of the integrals
  double precision              :: ss,ps,ds
  double precision              :: pxs,pys,pzs
  double precision              :: dxxs,dyys,dzzs,dxys,dxzs,dyzs

  double precision              :: pi,E,AB,AxBx,AyBy,AzBz,t,u,k

  pi = 4d0*atan(1d0)
  
! calculate the length AB between the two centers

  AB = (cGau(1)-cSla(1))**2d0 + (cGau(2)-cSla(2))**2d0 + (cGau(3)-cSla(3))**2d0
  AB = sqrt(AB)

  AxBx = (cGau(1)-cSla(1))/2d0
  AyBy = (cGau(2)-cSla(2))/2d0
  AzBz = (cGau(3)-cSla(3))/2d0

! intermediate variables

  t = expSla*sqrt(0.25d0/expGau)
  u = sqrt(expGau)*AB

  if(AB > 0d0) then 

!   (s|s) 
    ss = (1d0+t*(t-u))*erfc(t-u)*exp(2d0*t*(t-u)) - (1d0+t*(t+u))*erfc(t+u)*exp(2d0*t*(t+u))

!   (p|s)
    ps = (exp(t**2d0-2d0*t*u-u**2d0)*(4d0*exp(2d0*t*u)*(1d0+t**2d0)    & 
       + sqrt(pi)*t*(-(exp(t**2d0+u**2d0)*(3d0+2d0*t*(t-u))*erfc(t-u)) & 
       - exp(2d0*t*u+(t+u)**2d0)*(3d0+2d0*t*(t+u))*erfc(t+u))))/sqrt(pi)

!   (d|s)
    ds = (-8d0*exp(t**2d0-u**2d0)*u+4d0*exp(2d0*t*(t-u))*sqrt(pi)*t**2d0*((2d0+t**2d0-t*u)*erfc(t-u) &
       - exp(4d0*t*u)*(2d0+t*(t+u))*erfc(t+u)))/sqrt(pi)

!   backward scaling
    ds = 3d0*ss/u**5d0 - 3d0*ps/u**4d0 + ds/u**3d0
    ps = ps/u**2d0-ss/u**3d0
    ss = ss/u

  else

!   concentric case
    ss = (4d0*exp(t**2d0)*(1d0+t**2d0))/sqrt(pi)-2d0*exp(2d0*t**2d0)*t*(3d0+2d0*t**2d0)*erfc(t)
    ps = (8d0*exp(t**2d0)*(-1d0+4d0*t**2d0+2d0*t**4d0-exp(t**2)*sqrt(pi)*t**3d0*(5d0+2d0*t**2d0)*erfc(t)))/(3d0*sqrt(pi))

  endif

  k = expSla*sqrt(expGau)*t**3d0*exp(-t**2d0)*4d0*pi/expSla**(3d0/2d0)

! (s|s) 
  ss = k*ss

! (p|s) 
  ps = k*ps

  pxs = AxBx*ps
  pys = AyBy*ps
  pzs = AzBz*ps

! (d|s) 
  ds = k*ds

  dxxs = (2d0*ss+ps)/(4d0*expGau) + AxBx**2d0*ds
  dyys = (2d0*ss+ps)/(4d0*expGau) + AyBy**2d0*ds
  dzzs = (2d0*ss+ps)/(4d0*expGau) + AzBz**2d0*ds

  dxys = AxBx*AyBy*ds
  dxzs = AxBx*AzBz*ds
  dyzs = AyBy*AzBz*ds

! Print result
  write(*,'(A12,F16.10)') & 
    '(s|T|s) = ',ss
  write(*,'(A12,F16.10,3X,A12,F16.10,3X,A12,F16.10)') & 
    '(px|T|s) = ',pxs,'(py|T|s) = ',pys,'(pz|T|s) = ',pzs
  write(*,'(A12,F16.10,3X,A12,F16.10,3X,A12,F16.10,3X,A12,F16.10,3X,A12,F16.10,3X,A12,F16.10)') & 
    '(dx2|T|s) = ',dxxs,'(dy2|T|s) = ',dyys,'(dz2|T|s) = ',dzzs,'(dxy|T|s) = ',dxys,'(dxz|T|s) = ',dxzs,'(dyz|T|s) = ',dyzs

end
!*****************************************************************************


!*****************************************************************************
subroutine GauSlaNuclear(expGau,cGau,aGau,expSla,cSla,ZNuc,cNuc)

  implicit none

  BEGIN_DOC
  ! Compute the nuclear attraction integral between a Gaussian function
  ! with arbitrary angular momemtum and a s-type Slater function
  END_DOC

! Input variables 
  double precision,intent(in)   :: expGau,expSla
  double precision,intent(in)   :: cGau(3),cSla(3)
  integer,intent(in)            :: aGau(3)
  double precision,intent(in)   :: cNuc(3)
  double precision,intent(in)   :: ZNuc

! Final value of the overlap integral
  double precision              :: ss,ps,ds,fs
  double precision              :: pxs,pys,pzs

  double precision              :: pi,E,AB,x,y,k

  pi = 4d0*atan(1d0)
  E = exp(1d0)
  
! calculate the length AB between the two centers

  AB = (cGau(1)-cSla(1))**2d0 + (cGau(2)-cSla(2))**2d0 + (cGau(3)-cSla(3))**2d0
  AB = sqrt(AB)

! intermediate variables

  x = sqrt(expSla**2d0/(4d0*expGau))
  y = sqrt(expGau)*AB

  if(AB > 0d0) then 
    ss = (1d0+x*(x+y))*erfc(x+y)*exp(2d0*x*(x+y)) - (1d0+x*(x-y))*erfc(x-y)*exp(2d0*x*(x-y))
    ss = ss/y
  else
    ss = (4d0*E**x**2d0*(1d0+x**2d0))/sqrt(Pi)-2d0*E**(2d0*x**2d0)*x*(3d0+2d0*x**2d0)*Erfc(x)
  endif

  k = expSla*sqrt(expGau)*x**3d0*exp(-x**2)*4d0*pi/expSla**(3d0/2d0)
  ss = k*ss

! Print result
  write(*,*) ss

end
!*****************************************************************************
double precision function BoysF0(t)
  implicit none
  double precision, intent(in)  :: t
  double precision              :: pi

  pi = 4d0*atan(1d0)

  if(t > 0d0) then
    BoysF0 = 0.5d0*sqrt(pi/t)*erf(sqrt(t))
  else
    BoysF0 = 1d0
  endif

end
!*****************************************************************************


BEGIN_PROVIDER [ double precision, GauSlaOverlap_matrix, (ao_num, nucl_num) ]
  implicit none
  BEGIN_DOC
  ! <Gaussian | Slater> overlap matrix
  END_DOC
  integer                        :: i,j,k
  double precision               :: cGau(3)
  double precision               :: cSla(3)
  double precision               :: expSla, res, expGau
  integer                        :: aGau(3)

  do k=1,nucl_num
    cSla(1:3) = nucl_coord_transp(1:3,k)
    expSla    = slater_expo(k)

    do i=1,ao_num
      cGau(1:3) = nucl_coord_transp(1:3, ao_nucl(i))
      aGau(1:3) = ao_power(i,1:3)
      GauSlaOverlap_matrix(i,k) = 0.d0

      do j=1,ao_prim_num(i)
        expGau = ao_expo_ordered_transp(j,i)
        call GauSlaOverlap(expGau,cGau,aGau,expSla,cSla,res)
        GauSlaOverlap_matrix(i,k) += ao_coef_normalized_ordered_transp(j,i) * res
      enddo

    enddo

  enddo
  
END_PROVIDER

BEGIN_PROVIDER [ double precision, MOSlaOverlap_matrix, (mo_tot_num, nucl_num) ]
  implicit none
  BEGIN_DOC
! <MO | Slater>
  END_DOC
  call dgemm('N','N',mo_tot_num,nucl_num,ao_num,1.d0,                &
      mo_coef_transp, size(mo_coef_transp,1),                        &
      GauSlaOverlap_matrix, size(GauSlaOverlap_matrix,1),                    &
      0.d0, MOSlaOverlap_matrix, size(MOSlaOverlap_matrix,1))
END_PROVIDER
