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




