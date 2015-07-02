!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine aurora(iBlock)

  use ModGITM
  use ModSources

  implicit none

  integer, intent(in) :: iBlock

  call report("Aurora",1)
  call start_timing("Aurora")

  ! No aurora by default

  AuroralHeatingRate(:,:,:,iBlock)  = 0.0
  AuroralIonRateS(:,:,:,:,iBlock) = 0.0

  call end_timing("Aurora")

end subroutine aurora




