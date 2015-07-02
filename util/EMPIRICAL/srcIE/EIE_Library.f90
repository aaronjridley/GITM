!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!\
! This routine finds a point on in the Spherical system, given
! a Theta, Phi: 
! LocIn(1) = Phi
! LocIn(2) = Theta
! It returns a 4 element array:
! LocOut(1) = Index of Block
! LocOut(2) = Index of Longitude
! LocOut(3) = Index of Latitude
! LocOut(4) = Multiplication factor for Longitude
! LocOut(5) = Multiplication factor for Latitude
!/

subroutine EIE_FindPoint(LocIn, LocOut, iError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  real, dimension(2), intent(in)  :: LocIn
  real, dimension(5), intent(out) :: LocOut
  integer, intent(out) :: iError
  real :: MLTIn, LatIn
  integer :: j,i, iBLK

  logical :: IsFound

  LocOut = -1.0

  iError = 0

  !\
  ! Check to see if the point is even on the grid.
  !/

  MLTIn = mod(LocIn(1),24.0)

  LatIn = LocIn(2)
  if (LatIn > 90.0) then
     LatIn = 180.0 - LatIn
     MLTIn = mod(MLTIn+12.0,24.0)
  endif
  if (LatIn < -90.0) then
     LatIn = -180.0 - LatIn
     MLTIn = mod(MLTIn+12.0,24.0)
  endif

  if (MLTIn > 24.0 .or. MLTIn < 0 .or. LatIn > 90.0 .or. LatIn < -90.0) then
     iError = ecPointOutofRange_
     return
  endif

  iBLK = 1
  do while (iBLK <= EIEi_HavenBLKs)
     j = 1
     do while (j < EIEi_HavenMLTs)
        i = 1
        do while (i < EIEi_HavenLats)

           !\
           ! Check to see if the point is within the current cell
           !/

           if (LatIn <  EIEr3_HaveLats(j,i+1,iBLK) .and. &
               LatIn >= EIEr3_HaveLats(j,i,iBLK) .and. &
               MLTIn <  EIEr3_HaveMLTs(j+1,i,iBLK) .and. &
               MLTIn >= EIEr3_HaveMLTs(j,i,iBLK)) then

              !\
              ! If it is, then store the cell number and calculate
              ! the interpolation coefficients.
              !/

              LocOut(1) = iBLK
              LocOut(2) = j
              LocOut(3) = i

              LocOut(4) = (MLTIn                    -EIEr3_HaveMLTs(j,i,iBLK))/&
                          (EIEr3_HaveMLTs(j+1,i,iBLK)-EIEr3_HaveMLTs(j,i,iBLK))
              LocOut(5) = (LatIn                    -EIEr3_HaveLats(j,i,iBLK))/&
                          (EIEr3_HaveLats(j,i+1,iBLK)-EIEr3_HaveLats(j,i,iBLK))

              iBLK = EIEi_HavenBLKs
              j = EIEi_HavenMLTs
              i = EIEi_HavenLats

           endif

           i = i + 1

        enddo

        j = j + 1

     enddo

     iBLK = iBLK + 1

  enddo

end subroutine EIE_FindPoint
