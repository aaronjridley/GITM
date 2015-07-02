!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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




