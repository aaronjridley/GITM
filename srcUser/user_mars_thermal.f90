!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! ----------------------------------------------------------------
! If you want to output some specific variables, then do that here.
! In ModUserGITM, there are two variables defined, UserData2D and UserData3D.
! To output variables:
! 1. Figure out which variable you want to output.
! 2. Go into the code where the variable is set and copy it into
!    UserData3D or UserData2D.
! 3. Do this for each variable you want to output.
! 4. Edit output_header_user below, making a list of all the variables
!    that you have added.  Make sure you leave longitude, latitude, and
!    altitude in the list of variables.
! 5. Count the number of variables that you output (including 
!    Lon, Lat, and Alt). Change nVarsUser3d or nVarsUser2d in the 
!    subroutines towards the top of this file.
! 6. If you add more than 40 variables, you probably should check 
!    nUserOutputs in ModUserGITM.f90 and make sure that this number is
!    larger than the number of variables that you added.
! 7. Recompile and run. Debug. Repeat 7.
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser3d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser3d = 4

  if (nVarsUser3d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser3d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser1d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser1d = 12

  if (nVarsUser1d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser1d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser2d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser2d = 6

  if (nVarsUser2d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser2d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine output_header_user(cType, iOutputUnit_)

  use ModUserGITM

  implicit none

  character (len=5), intent(in) :: cType
  integer, intent(in)           :: iOutputUnit_

  ! ------------------------------------------
  ! 3D Output Header
  ! ------------------------------------------
 if (cType(1:2) == '3D') then 

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser3d, " nvars"
     write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") nLats+4, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") nLons+4, " nLongitudes"

     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Joule Heating"

  endif
  
  ! ------------------------------------------
  ! 1D Output Header
  ! ------------------------------------------
  if (cType(1:2) == '1D') then 

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser1d, " nvars"
     write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") 1, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") 1, " nLongitudes"

     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Temperature"
     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "VelGradTemp"
     write(iOutputUnit_,"(I7,A1,a)")  6, " ", "TempDivVel"
     write(iOutputUnit_,"(I7,A1,a)")  7, " ", "NumDiffTemp"
     write(iOutputUnit_,"(I7,A1,a)")  8, " ", "StressHeating"
     write(iOutputUnit_,"(I7,A1,a)")  9, " ", "LowAtmosRad"
     write(iOutputUnit_,"(I7,A1,a)")  10, " ", "EuvHeating"
     write(iOutputUnit_,"(I7,A1,a)")  11, " ", "RadCooling"     
     write(iOutputUnit_,"(I7,A1,a)")  12, " ", "Conduction"     

     
  endif

  ! ------------------------------------------
  ! 2D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '2D') then 

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser2d, " nvars"
     write(iOutputUnit_,"(I7,7A)")     1, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") nLats, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") nLons, " nLongitudes"

     write(iOutputUnit_,*) ""
     write(iOutputUnit_,*) "NO GHOSTCELLS"
     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Potential (kV)"
     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Average Energy (keV)"
     write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Total Energy (ergs)"

  endif

  write(iOutputUnit_,*) ""

end subroutine output_header_user

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon

  do iAlt=-1,nAlts+2
     do iLat=-1,nLats+2
        do iLon=-1,nLons+2
           write(iOutputUnit_)       &
                Longitude(iLon,iBlock), &
                Latitude(iLat,iBlock), &
                Altitude_GB(iLon, iLat, iAlt, iBlock),&
                UserData3D(iLon,iLat,iAlt,1:nVarsUser3d-3,iBlock)
        enddo
     enddo
  enddo

end subroutine output_3dUser

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_1dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon, iiAlt

  write(*,*) "temp: ",UserData1D(1,1,42,1),Temperature(1,1,42,1),TempUnit(1,1,42)
  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)

     write(iOutputUnit_)       &
          Longitude(1,iBlock), &
          Latitude(1,iBlock), &
          Altitude_GB(1, 1, iAlt, iBlock),&
          UserData1D(1,1,iiAlt,1:nVarsUSer1d-3)
  enddo

end subroutine output_1dUser
!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_2dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon

  iAlt = 1
  do iLat=1,nLats
     do iLon=1,nLons
        write(iOutputUnit_)       &
             Longitude(iLon,iBlock), &
             Latitude(iLat,iBlock), &
             Altitude_GB(iLon, iLat, iAlt, iBlock),&
             UserData2D(iLon,iLat,iAlt,1:nVarsUser2d-3,iBlock)
     enddo
  enddo

end subroutine output_2dUser

