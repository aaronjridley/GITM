
module ModRates

  use ModGITM, only: nLons, nLats, nAlts

  implicit none

!    Temperature Independent (i.e. Static) Reaction Rates

  ! From Haider, 1997
  integer, parameter  :: iRrR1_    =  1
  integer, parameter  :: iRrR2_    =  2
  integer, parameter  :: iRrR3_    =  3
  integer, parameter  :: iRrR4_    =  4
  integer, parameter  :: iRrR5_    =  5
  integer, parameter  :: iRrR7_    =  6
  integer, parameter  :: iRrR9_    =  7
  integer, parameter  :: iRrR10_   =  8
  integer, parameter  :: iRrR11a_  =  9
  integer, parameter  :: iRrR11b_  =  10
  integer, parameter  :: iRrR12_   =  11
  integer, parameter  :: iRrR13_   =  12
  integer, parameter  :: iRrR14_   =  13
  integer, parameter  :: iRrR15_   =  14
  integer, parameter  :: iRrR16_   =  15
  integer, parameter  :: iRrR17_   =  16

  integer, parameter  :: nRrTempInd = 25

  real :: RrTempInd(nRrTempInd)

!    Temperature Depenedent

  integer, parameter  :: iRrK1_     =  1
  integer, parameter  :: nRrTempDep = 2

  real :: RrTempDep(nLons, nLats, nAlts, nRrTempDep)

!    Collision Frequencies

  integer, parameter  :: iCfIN_ = 1
  integer, parameter  :: iCfEN_ = 2
  integer, parameter  :: iCfEI_ = 3

  integer, parameter  :: nCf = 3

  real :: Cf(nLons, nLats, nAlts, nCf)

contains

  subroutine set_RrTempInd

    ! These are from Haider, 1997

    !O2+ + N -> NO+ + O
    RrTempInd(iRrR1_)    = 1.8e-10 / 1.0e6 ! cm-3 to m-3

    !O2+ + NO -> NO+ + O2
    RrTempInd(iRrR2_)    = 4.8e-10 / 1.0e6 ! cm-3 to m-3

    !N2+ + O -> NO+ + N
    RrTempInd(iRrR3_)    = 1.4e-10 / 1.0e6 ! cm-3 to m-3

    !CO2+ + N -> NO+ + CO
    RrTempInd(iRrR4_)    = 1.0e-11 / 1.0e6 ! cm-3 to m-3

    !CO2+ + NO -> NO+ + CO2
    RrTempInd(iRrR5_)    = 1.2e-10 / 1.0e6 ! cm-3 to m-3

    !CO+ + NO -> NO+ + CO
    RrTempInd(iRrR7_)    = 3.3e-10 / 1.0e6 ! cm-3 to m-3

    !CO+ + CO2 -> CO2+ + CO
    RrTempInd(iRrR9_)    = 1.1e-9 / 1.0e6 ! cm-3 to m-3

    !N2+ + CO2 -> CO2+ + N2
    RrTempInd(iRrR10_)    = 8.0e-10 / 1.0e6 ! cm-3 to m-3

  end subroutine set_RrTempInd


end module ModRates
