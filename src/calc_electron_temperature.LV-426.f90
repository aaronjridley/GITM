!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_electron_temperature(iBlock)

  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  call report("Calc_electron_temperature",1)

  eTemperature(:,:,:,iBlock)= Temperature(:,:,:,iBlock)*TempUnit * 2.0
  ITemperature(:,:,:,iBlock)= Temperature(:,:,:,iBlock)*TempUnit * 1.5

  return

end subroutine calc_electron_temperature
