!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModRates

  use ModGITM, only: nLons, nLats, nAlts

  implicit none

!    Temperature Independent (i.e. Static) Reaction Rates

  ! From Bougher
  integer, parameter  :: iRrR4_    =  1
  integer, parameter  :: iRrR5_    =  2
  integer, parameter  :: iRrR6_    =  3
  integer, parameter  :: iRrR7_    =  4
  integer, parameter  :: iRrR8_    =  5
  integer, parameter  :: iRrR9_    =  6
  integer, parameter  :: iRrR10_   =  7
  integer, parameter  :: iRrR13_   =  8
  integer, parameter  :: nRrTempInd = 8

  real :: RrTempInd(nRrTempInd)

!    Temperature Depenedent

  integer, parameter  :: iRrK1_     =  1
  integer, parameter  :: iRrK2_     =  2
  integer, parameter  :: iRrK3_     =  3
  integer, parameter  :: iRrK11_    =  4
  integer, parameter  :: iRrK12_    =  5
  integer, parameter  :: nRrTempDep =  5

  real :: RrTempDep(nLons, nLats, nAlts, nRrTempDep)

!    Collision Frequencies

  integer, parameter  :: iCfIN_ = 1
  integer, parameter  :: iCfEN_ = 2
  integer, parameter  :: iCfEI_ = 3

  integer, parameter  :: nCf = 3

  real :: Cf(nLons, nLats, nAlts, nCf)

contains

  subroutine set_RrTempInd

    ! These are from Bougher

    !H2+ + H -> H+ + H2
    RrTempInd(iRrR4_)    = 6.4e-10 / 1.0e6 ! cm-3 to m-3

    !H2+ + H2 -> H3+ + H
    RrTempInd(iRrR5_)    = 2.0e-9 / 1.0e6 ! cm-3 to m-3

    !H2+ + He -> HeH+ + H
    RrTempInd(iRrR6_)    = 1.4e-10 / 1.0e6 ! cm-3 to m-3

    !H2+ + CH4 -> CH5+ + H
    RrTempInd(iRrR7_)    = 1.1e-11 / 1.0e6 ! cm-3 to m-3

    !H2+ + CH4 -> CH4+ + H2
    RrTempInd(iRrR8_)    = 1.4e-10 / 1.0e6 ! cm-3 to m-3

    !H2+ + CH4 -> CH3+H + H2
    RrTempInd(iRrR9_)    = 2.3e-10 / 1.0e6 ! cm-3 to m-3

    !H+ + 2H2 -> H3+ + H2
    RrTempInd(iRrR10_)    = 3.2e-9 / 1.0e6 ! cm-3 to m-3

    !H3+ + CH4 -> CH5+ + H2
    RrTempInd(iRrR13_)    = 2.4e-9 / 1.0e6 ! cm-3 to m-3

  end subroutine set_RrTempInd

end module ModRates
