! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine aurora(iBlock)

  use ModGITM
  use ModSources
  use ModTime, only : tSimulation
  use ModInputs
  use ModConstants
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock


  if (floor((tSimulation - dT)/dTAurora) == &
       floor(tSimulation/dTAurora)) return

  AuroralIonRateS                  = 0.0
  AuroralHeatingRate(:,:,:,iBlock) = 0.0

  return

end subroutine aurora




