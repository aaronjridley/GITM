! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

module ModTime

  use ModKind, ONLY: Real8_
  implicit none

  ! Time variables

  real                  :: tSimulation = 0.0
  integer, dimension(7) :: iTimeArray
  real(Real8_)          :: CurrentTime, EndTime, StartTime, VernalTime
  real(Real8_)          :: RestartTime
  real(Real8_)          :: PauseTime
  real                  :: utime
  integer               :: iJulianDay, iDay
  integer               :: iStep = 1

  integer, parameter :: iYear_   = 1
  integer, parameter :: iMonth_  = 2
  integer, parameter :: iDay_    = 3
  integer, parameter :: iHour_   = 4
  integer, parameter :: iMinute_ = 5
  integer, parameter :: iSecond_ = 6

end module ModTime
