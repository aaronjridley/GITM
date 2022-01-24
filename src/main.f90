! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!-----------------------------------------------------------------------------
! $Id: main.f90,v 1.10 2013/10/24 18:36:35 agburr Exp $
!
! Author: Aaron Ridley
!
! Comments: GITM main
!
! AGB 3/31/13: Added call for RCMR data assimilation
! AGB 10/23/13: Adapted RCMR call to new format
!-----------------------------------------------------------------------------

program GITM
  use ModRCPE
  use ModInputs
  use ModTime
  use ModGITM
  use ModMpi
  use ModSatellites, only: SatCurrentDat, SatAltDat, nRCMRSat
  use ModEUV, only: sza

  implicit none

  integer :: iBlock, iTimeToCheck
  real :: currentTimeTemp, tTemp
  

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------
  
  call init_mpi
  call start_timing("GITM")
  call delete_stop

  call init_planet
  call set_defaults
   
  call read_inputs(cInputFile)

  call set_inputs

  call initialize_gitm(CurrentTime)

  call write_output

  call report("Starting Main Time Loop",0)

  ! ------------------------------------------------------------------------
  ! Run for a few iterations
  ! ------------------------------------------------------------------------

  do while (CurrentTime < EndTime)

     call calc_pressure

     !!! We may have to split cMax and Dt calculation!!!
     if(RCMRFlag) then
        Dt = 1 !was 2
     else
        Dt = FixedDt
     end if

     !tTemp = tSimulation + dt
     !currentTimeTemp = StartTime + tTemp
     !call time_real_to_int2(currentTimeTemp, iTimeToCheck)

     !if (iTimeToCheck .eq. 1) then
     !  Dt = 1
     !endif
    
     call calc_timestep_vertical
     if (.not. Is1D) call calc_timestep_horizontal

     if(RCMRFlag .and. iStep .le. rcmrStepToTurnOff) then
        call mainRCPE
        if (mod((istep-1)*Dts,Measure_Dts) == 0.0 .and. iProc .eq. 5) then
          if (coefficientToEstimate == "AO2") then
            write(*,*) "ThermalConduction_" // coefficientToEstimate // ":", &
                       ThermalConduction_AO2
          else if (coefficientToEstimate == "AO") then
            write(*,*) "ThermalConduction_" // coefficientToEstimate // ":", &
                       ThermalConduction_AO
          else if (coefficientToEstimate == "s") then
            write(*,*) "ThermalConduction_" // coefficientToEstimate // ":", &
                       ThermalConduction_s
          endif
        endif
     endif  

     call advance

     if (.not.IsFramework) then
        call check_stop
     endif

     iStep = iStep + 1

     call write_output
  end do

  ! ------------------------------------------------------------------------
  ! Finish run
  ! ------------------------------------------------------------------------

  call finalize_gitm

end program GITM

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  call stop_gitm(StringError)

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  implicit none
  integer, intent(in) :: iUnit

  return

end subroutine CON_io_unit_new

subroutine time_real_to_int2(timereal, itime)

  implicit none

  integer, dimension(1:12) :: dayofmon
  integer :: itime
  double precision :: timereal
  integer :: nyear, nleap, nmonth, nday, nhour, nmin, nsec
  double precision :: speryear = 31536000.0
  double precision :: sperday = 86400.0
  double precision :: sperhour = 3600.0
  double precision :: spermin = 60.0
  double precision :: timeleft

  dayofmon(1) = 31
  dayofmon(2) = 28
  dayofmon(3) = 31
  dayofmon(4) = 30
  dayofmon(5) = 31
  dayofmon(6) = 30
  dayofmon(7) = 31
  dayofmon(8) = 31
  dayofmon(9) = 30
  dayofmon(10) = 31
  dayofmon(11) = 30
  dayofmon(12) = 31

  nyear = int(timereal/speryear)
  nleap = nyear/4
  nday = int((timereal - (dble(nyear)*speryear))/sperday)

  if (nday.le.nleap) then
     nyear = int((timereal - (dble(nleap)*sperday))/speryear)
     nleap = nyear/4
     nday = int((timereal - (dble(nyear)*speryear))/sperday)
     if (nday.le.nleap) then
        nyear = int((timereal - (dble(nleap)*sperday))/speryear)
        nleap = nyear/4
        nday = int((timereal - (dble(nyear)*speryear))/sperday)
     endif
  endif

  if (mod((nyear+65),4).eq.0) dayofmon(2) = dayofmon(2) + 1

  nday = nday - nleap

  timeleft = timereal - dble(nyear)*speryear
  timeleft = timeleft - dble(nday+nleap)*sperday

  nhour = int(timeleft/sperhour)
  timeleft = timeleft - dble(nhour)*sperhour

  nmin = int(timeleft/spermin)
  timeleft = timeleft - dble(nmin)*spermin

  nsec = int(timeleft)

  nmonth = 1;

  do while (nday.ge.dayofmon(nmonth))
     nday = nday - dayofmon(nmonth)
     nmonth = nmonth + 1
  end do

  itime = nsec


end subroutine time_real_to_int2

!---------------------------------------------------------------------------

