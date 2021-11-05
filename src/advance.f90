! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!------------------------------------------------------------------------------
! $Id: advance.f90,v 1.18 2017/10/30 14:05:36 ridley Exp $
! Author: Aaron Ridley, UMichigan
!
! Modified:
!           AGB, Oct 2013 - Adapted to allow new RCMR format
!           Asad, Feb 2013 - Adapted to update F10.7 using estimated value if
!                            RCMR data assimilation is being used.
!
! Comments: A routine to initialize the model for the next timestep
!------------------------------------------------------------------------------

subroutine advance

  use ModRCMR, only: RCMRFlag, RCMROutType
  use ModConstants
  use ModGITM
  use ModTime
  use ModEUV
  use ModInputs
  use ModIndicesInterfaces
  use ModMpi

  implicit none

  integer, external :: jday

  integer :: iBlock, iAlt, iLat, iLon,ispecies, iError
  real*8 :: DTime

  if(RCMRFlag) then
     if(RCMROutType == 'F107') then
        ! Asad: When running RCAC both the F10.7 and F10.7A should have the same
        !       number at this point because they are being estimated together

        call IO_set_f107_single(f107_est)
        call IO_set_f107a_single(f107_est)
     else if(RCMROutType == "PHOTOELECTRON") then
        PhotoElectronHeatingEfficiency = PhotoElectronHeatingEfficiency_est 
     else if(RCMROutType == "EDC") then
        EddyDiffusionCoef = EDC_est(1,1)  !Ankit23May16: Added EDC_est out
     end if
  end if

  call report("advance",1)
  call start_timing("advance")

  if (UseGSWMTides)  call update_tides
  if (UseWACCMTides) call update_waccm_tides
  if (UseBcPerturbation) call get_mean_bcs
  if (UsePerturbation) call user_perturbation

  if (.not. UseStatisticalModelsOnly) then

     call advance_vertical_all
     call add_sources 
     if (.not. Is1D) call advance_horizontal_all

  else

     Dt = DtStatisticalModels

  endif

  if (iDebugLevel > 0) &
       write(*,*) "=> MaxTemp : ",maxval(temperature)*TempUnit(1,1,nalts)

  tSimulation = tSimulation + dt
  CurrentTime = StartTime + tSimulation

  call time_real_to_int(CurrentTime, iTimeArray)

  DTime = CurrentTime - VernalTime
  do while (DTime > DaysPerYearInput*RotationPeriodInput)
     VernalTime = VernalTime+int(DaysPerYearInput)*RotationPeriodInput
     DTime = CurrentTime - VernalTime
  enddo
  iDay  = DTime / RotationPeriodInput
  uTime = (DTime / RotationPeriodInput - iDay) * RotationPeriodInput

  iJulianDay = jday(iTimeArray(1), iTimeArray(2), iTimeArray(3)) 

  if (UseStatisticalModelsOnly) then

     iError = 0
     call get_f107(CurrentTime, f107, iError)
     if (iError /= 0) then
        write(*,*) "Error in getting F107 value.  Is this set?"
        write(*,*) "Code : ",iError
        call stop_gitm("Stopping in advance")
     endif

     call get_f107a(CurrentTime, f107a, iError)
     if (iError /= 0) then
        write(*,*) "Error in getting F107a value.  Is this set?"
        write(*,*) "Code : ",iError
        call stop_gitm("Stopping in advance")
     endif

     write(*,*) "F10.7 = ", f107, f107a
     call init_msis
     call init_iri
     call init_b0

  endif

  call end_timing("advance")

end subroutine advance
