! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

Module ModTimingGITM

  integer, parameter :: nTimingMax = 100
  character (LEN = 20), dimension(nTimingMax) :: cTimingNames
  real*8, dimension(nTimingMax, 2) :: Timings
  logical, dimension(nTimingMax) :: IsTiming
  integer, dimension(nTimingMax) :: iTimingLevel
  integer :: nTiming = 0

end Module ModTimingGITM

