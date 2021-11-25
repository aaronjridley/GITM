! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_electron_temperature(iBlock)

  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  call report("Calc_electron_temperature",1)

  eTemperature(:,:,:,iBlock)= Temperature(:,:,:,iBlock)*TempUnit * 2.0
  ITemperature(:,:,:,iBlock)= Temperature(:,:,:,iBlock)*TempUnit * 1.5

  return

end subroutine calc_electron_temperature
