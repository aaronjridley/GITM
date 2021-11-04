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

  use ModRCMR, only: RCMRFlag
  use ModGITM, only: dt
  use ModInputs

  implicit none

  call report("advance",1)
  call start_timing("advance")

  if (RCMRFlag) call set_RCMR_estimations
  if (UseGSWMTides) call update_tides
  if (UseWACCMTides) call update_waccm_tides
  if (UseBcPerturbation) call get_mean_bcs
  if (UsePerturbation) call user_perturbation

  if (.not. UseStatisticalModelsOnly) then

     call advance_vertical_all
     call add_sources 
     if (.not. Is1D) call advance_horizontal_all

  else
     dt = DtStatisticalModels
  endif

  ! Increment time and all the time associated variables
  call update_time
  
  if (UseStatisticalModelsOnly) then
     call init_msis
     call init_iri
     call init_b0
  endif

  call end_timing("advance")

end subroutine advance
