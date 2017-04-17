
subroutine create_selection_buffer(N, siz, res)
  use selection_types
  implicit none

  integer, intent(in) :: N, siz
  type(selection_buffer), intent(out) :: res

  allocate(res%det(N_int, 2, siz), res%val(siz))

  res%val = 0d0
  res%det = 0_8
  res%N = N
  res%mini = 0d0
  res%cur = 0
end subroutine

  
subroutine add_to_selection_buffer(b, det, val)
  use selection_types
  implicit none

  type(selection_buffer), intent(inout) :: b
  integer(bit_kind), intent(in) :: det(N_int, 2)
  double precision, intent(in) :: val
  integer :: i

  if(dabs(val) >= b%mini) then
    b%cur += 1
    b%det(1:N_int,1:2,b%cur) = det(1:N_int,1:2)
    b%val(b%cur) = val
    if(b%cur == size(b%val)) then
      call sort_selection_buffer(b)
    end if
  end if
end subroutine


subroutine sort_selection_buffer(b)
  use selection_types
  implicit none

  type(selection_buffer), intent(inout) :: b
  double precision, allocatable::  absval(:)
  integer, allocatable :: iorder(:)
  double precision, pointer :: vals(:)
  integer(bit_kind), pointer :: detmp(:,:,:)
  integer :: i, nmwen
  logical, external :: detEq
  nmwen = min(b%N, b%cur)


  allocate(iorder(b%cur), detmp(N_int, 2, size(b%det,3)), absval(b%cur), vals(size(b%val)))
  absval = -dabs(b%val(:b%cur))
  do i=1,b%cur
    iorder(i) = i
  end do
  ! Optimal for almost sorted data
!  call sorted_dnumber(absval, b%cur, i)
!  if (b%cur/i > 
!  call insertion_dsort(absval, iorder, b%cur)
  call dsort(absval, iorder, b%cur)
  do i=1, nmwen
    detmp(1:N_int,1,i) = b%det(1:N_int,1,iorder(i))
    detmp(1:N_int,2,i) = b%det(1:N_int,2,iorder(i))
    vals(i) = b%val(iorder(i))
  end do
  do i=nmwen+1, size(vals)
    vals(i) = 0.d0
  enddo
  deallocate(b%det, b%val)
  b%det => detmp
  b%val => vals
  b%mini = max(b%mini,dabs(b%val(b%N)))
  b%cur = nmwen
end subroutine

