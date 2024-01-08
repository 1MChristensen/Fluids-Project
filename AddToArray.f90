module DynamicalArrays
  use iso_fortran_env

  contains

  subroutine AddToList(list, element)

    IMPLICIT NONE

    integer :: i, isize
    real(real64), intent(in) :: element
    real(real64), dimension(:), allocatable, intent(inout) :: list
    real(real64), dimension(:), allocatable :: clist


    if(allocated(list)) then
        isize = size(list)
        allocate(clist(isize+1))
        do i=1,isize          
        clist(i) = list(i)
        end do
        clist(isize+1) = element

        deallocate(list)
        call move_alloc(clist, list)

    else
        allocate(list(1))
        list(1) = element
    end if


  end subroutine AddToList


end module DynamicalArrays
