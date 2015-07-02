!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

Module ModTimingGITM

  integer, parameter :: nTimingMax = 100
  character (LEN = 20), dimension(nTimingMax) :: cTimingNames
  real*8, dimension(nTimingMax, 2) :: Timings
  logical, dimension(nTimingMax) :: IsTiming
  integer, dimension(nTimingMax) :: iTimingLevel
  integer :: nTiming = 0

end Module ModTimingGITM

