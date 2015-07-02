!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine Insert_into_Indices_Array(array, label_)

  use ModIndices
  use ModTimeConvert, ONLY: time_int_to_real
  implicit none

  real, dimension(6,MaxIndicesEntries), intent(in) :: array
  integer, intent(in) :: label_

  real*8 :: time_now
  integer, dimension(7) :: itime

  integer :: iStart, iEnd, i, j, nPts, nPtsSub
  real    :: dt, f107a

  npts = Number_of_Points(array)

  if (npts == -1) then

     Indices_TV(:,label_) = -1.0e32
     nIndices_V(label_) = 0

  elseif (npts == 0) then

     itime(1) = array(1,1)
     itime(2) = array(2,1)
     itime(3) = array(3,1)
     itime(4) = array(4,1)
     itime(5) = array(5,1)
     itime(6) = 0
     itime(7) = 0
     call time_int_to_real(itime, IndexTimes_TV(1,label_))
     Indices_TV(1,label_) = array(6,1)
     nIndices_V(label_) = 1

  else

     do i = 1, npts

        itime(1) = array(1,i)
        itime(2) = array(2,i)
        itime(3) = array(3,i)
        itime(4) = array(4,i)
        itime(5) = array(5,i)
        itime(6) = 0
        itime(7) = 0
        call time_int_to_real(itime, IndexTimes_TV(i,label_))
        Indices_TV(i,label_) = array(6,i)

     enddo

     nIndices_V(label_) = npts

  endif

  if (label_ == f107_) then

     if (npts == -1) then

        Indices_TV(:,f107a_) = -1.0e32
        nIndices_V(f107a_) = 0

     elseif (npts == 0) then

        IndexTimes_TV(1,f107a_) = IndexTimes_TV(1,f107_)
        Indices_TV(1,f107a_) = Indices_TV(1,f107_)
        nIndices_V(f107a_) = 1

     else

        do i = 1, npts

           dt = 40.0 * 24.0 * 3600.0  ! want 40 days on either side

           iStart = i
           do while (iStart > 1 .and. &
                IndexTimes_TV(i,f107_)-dt < IndexTimes_TV(iStart,f107_))
              iStart = iStart - 1
           enddo

           iEnd = i
           do while (iEnd < npts .and. &
                IndexTimes_TV(i,f107_)+dt > IndexTimes_TV(iEnd,f107_))
              iEnd = iEnd + 1
           enddo

           nPtsSub = 0
           f107a = 0
           time_now = 0
           do j = iStart, iEnd
              f107a = f107a + Indices_TV(j,f107_)
              time_now = time_now + IndexTimes_TV(j,f107_)
              nPtsSub = nPtsSub + 1
           enddo

           IndexTimes_TV(i,f107a_) = time_now/nPtsSub
           Indices_TV(i,f107a_) = f107a/nPtsSub

        enddo

        nIndices_V(f107a_) = npts

     endif

  endif

contains

  integer function Number_of_Points(array)

    real, dimension(6,MaxIndicesEntries), intent(in) :: array
    integer :: n

    n = 0

    if (maxval(abs(array)) == 0) then

       Number_of_Points = -1

    else

       if (maxval(abs(array(:,2:MaxIndicesEntries))) == 0 .and. &
            array(1,1) /= 0.0 .and. array(4,1) == 0.0) then

          Number_of_Points = 0

       else

          Number_of_Points = 1

          do while (maxval(abs(array(:,Number_of_Points))) > 0)
             Number_of_Points = Number_of_Points + 1
          enddo

          Number_of_Points = Number_of_Points - 1

       endif

    endif

  end function Number_of_Points

end subroutine Insert_into_Indices_Array
