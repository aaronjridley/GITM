! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

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




