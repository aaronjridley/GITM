!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModSort

  implicit none

  private ! except

  public :: sort_quick ! the quick sort algoritm (Numerical Recipes)
  public :: sort_sum   ! add up numbers but first sort them by magnitude
  public :: sort_test  ! unit test

  public :: sort_quick_func

  ! Local variables for the unit test

  integer, parameter :: n=5
  real    :: a_I(n)

contains
  !===========================================================================
  subroutine sort_quick(n, arr, indx)

    ! Quick sort algorithm: sorts indx array according to 'arr'
    ! so that arr(indx(i)) <= arr(indx(i+1)) for i=1..n-1

    ! Based on the F77 code indexx from Numerical Recipes

    integer, intent(in)  :: n
    real, intent(in)     :: arr(n)
    integer, intent(out) :: indx(n)

    integer, parameter   :: M=7, NSTACK=1000 ! should be > 2*log_2(n)

    integer :: i,indxt,ir,itemp,j,jstack,k,l, iStack(NSTACK)
    real    :: a

    !-------------------------------------------------------------------
    do j=1,n
       indx(j)=j
    end do
    jstack=0
    l=1
    ir=n
1   if(ir-l < M)then
       do j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do i=j-1,1,-1
             if(arr(indx(i)).le.a)goto 2
             indx(i+1)=indx(i)
          end do
          i=0
2         indx(i+1)=indxt
       end do
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       itemp=indx(k)
       indx(k)=indx(l+1)
       indx(l+1)=itemp
       if(arr(indx(l+1)) > arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
       endif
       if(arr(indx(l)) > arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
       endif
       if(arr(indx(l+1)) > arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
       endif
       i=l+1
       j=ir
       indxt=indx(l)
       a=arr(indxt)
3      continue
       i=i+1
       if(arr(indx(i)) < a)goto 3
4      continue
       j=j-1
       if(arr(indx(j)) > a)goto 4
       if(j < i)goto 5
       itemp=indx(i)
       indx(i)=indx(j)
       indx(j)=itemp
       goto 3
5      indx(l)=indx(j)
       indx(j)=indxt
       jstack=jstack+2
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif
    goto 1

  end subroutine sort_quick

  !============================================================================

  subroutine sort_quick_func(n, is_larger, iSort_I)

    ! Quick sort algorithm: sorts iSort_I array according to is_larger function
    ! so that is_larger(iSort_I(i+1), iSort_I(i)) is true for i=1..n-1

    ! Based on the F77 code indexx from Numerical Recipes

    integer, intent(in)  :: n
    interface
       logical function is_larger(i,j)
         integer, intent(in):: i,j
       end function is_larger
    end interface
    integer, intent(out) :: iSort_I(n)

    integer, parameter   :: M=7, NSTACK=1000 ! should be > 2*log_2(n)

    integer :: i,indxt,ir,itemp,j,jstack,k,l, iStack(NSTACK)

    !-------------------------------------------------------------------
    do j=1,n
       iSort_I(j)=j
    end do
    jstack=0
    l=1
    ir=n
1   if(ir-l < M)then
       do j=l+1,ir
          indxt=iSort_I(j)
          do i=j-1,1,-1
             if(.not.is_larger(iSort_I(i), indxt)) goto 2
             iSort_I(i+1)=iSort_I(i)
          end do
          i=0
2         iSort_I(i+1)=indxt
       end do
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       itemp=iSort_I(k)
       iSort_I(k)=iSort_I(l+1)
       iSort_I(l+1)=itemp
       if(is_larger(iSort_I(l+1), iSort_I(ir)))then
          itemp=iSort_I(l+1)
          iSort_I(l+1)=iSort_I(ir)
          iSort_I(ir)=itemp
       endif
       if(is_larger(iSort_I(l), iSort_I(ir)))then
          itemp=iSort_I(l)
          iSort_I(l)=iSort_I(ir)
          iSort_I(ir)=itemp
       endif
       if(is_larger(iSort_I(l+1), iSort_I(l)))then
          itemp=iSort_I(l+1)
          iSort_I(l+1)=iSort_I(l)
          iSort_I(l)=itemp
       endif
       i=l+1
       j=ir
       indxt=iSort_I(l)
3      continue
       i=i+1
       if(is_larger(indxt, iSort_I(i))) goto 3
4      continue
       j=j-1
       if(is_larger(iSort_I(j), indxt)) goto 4
       if(j < i)goto 5
       itemp=iSort_I(i)
       iSort_I(i)=iSort_I(j)
       iSort_I(j)=itemp
       goto 3
5      iSort_I(l)=iSort_I(j)
       iSort_I(j)=indxt
       jstack=jstack+2
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif
    goto 1

  end subroutine sort_quick_func

  !============================================================================

  real function sort_sum(a_I)

    ! add up a_I but sort it first based on the magnitude of elements

    real, intent(in):: a_I(:)

    real :: SortSum
    integer :: i, n
    integer, allocatable:: i_I(:)
    !-----------------------------------------------------------------------
    n = size(a_I)
    allocate(i_I(n))
    call sort_quick(n, abs(a_I), i_I)
    SortSum = 0.0
    do i = n, 1, -1
       if(i_I(i) < 0)write(*,*)'This avoids ifort 12 optimization error!'
       SortSum = SortSum + a_I(i_I(i))
    end do
    sort_sum = SortSum
    deallocate(i_I)

  end function sort_sum

  !======================================================================
  logical function is_larger_test(i, j)
    integer, intent(in):: i, j

    is_larger_test = a_I(i) > a_I(j)

  end function is_larger_test
  !======================================================================
  subroutine sort_test

    real    :: b_I(n), SortSum
    integer :: i_I(n), i, iTest
    logical :: IsError
    !------------------------------------------------------------------------
    a_I = (/0.0, 2.0, 1.0, 4.0, 0.0/)

    do iTest = 1, 2
       if(iTest == 1)then
          write(*,'(a)')'Testing sort_quick'
          call sort_quick(n,a_I,i_I)
       else
          write(*,'(a)')'Testing sort_quick_func'
          call sort_quick_func(n,is_larger_test,i_I)
       end if
       b_I = a_I(i_I)

       IsError = .false.
       do i=2,n
          if(b_I(i-1) > b_I(i))then
             write(*,*)'Error at index i=',i,': sorted b_I(i-1)=',b_I(i-1),&
                  ' should not exceed b_I(i)=',b_I(i)
             IsError = .true.
          end if
       end do

       if(IsError)then
          write(*,'(a,5f5.0)')'original array =',a_I
          write(*,'(a,5i5)'  )'sorted index   =',i_I
          write(*,'(a,5f5.0)')'sorted array   =',b_I
       end if
    end do

    write(*,'(a)')'Testing sort_sum'
    a_I = (/1.e-6, 1.e-7, 1.e10, -3.e9, -7.e9/)

    SortSum = sort_sum(a_I)
    if(abs(SortSum - 1.1e-6) > 1.1e-12) &
         write(*,*)'Error: SortSum should be 1.1e-6, but it is ',SortSum

  end subroutine sort_test

end module ModSort
