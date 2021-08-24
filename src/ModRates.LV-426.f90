! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

module ModRates

  use ModGITM, only: nLons, nLats, nAlts

  implicit none

!    Temperature Independent (i.e. Static) Reaction Rates


  integer, parameter  :: nRrTempInd = 1

  real :: RrTempInd(nRrTempInd)

!    Temperature Depenedent

  integer, parameter  :: nRrTempDep = 1

  real :: RrTempDep(nLons, nLats, nAlts, nRrTempInd)

!    Collision Frequencies

  integer, parameter  :: iCfIN_ = 1
  integer, parameter  :: iCfEN_ = 2
  integer, parameter  :: iCfEI_ = 3

  integer, parameter  :: nCf = 3

  real :: Cf(nLons, nLats, nAlts, nCf)

contains

  subroutine set_RrTempInd

    ! Do nothing yet!

    return


  end subroutine set_RrTempInd


end module ModRates
