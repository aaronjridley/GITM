!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModRates

  use ModGITM, only: nLons, nLats, nAlts

  implicit none

!    Temperature Independent (i.e. Static) Reaction Rates

  integer, parameter  :: iRrK4_    =  1
  integer, parameter  :: iRrK5_    =  2
  integer, parameter  :: iRrK6_    =  3
  integer, parameter  :: iRrK7_    =  4
  integer, parameter  :: iRrK8_    =  5
  integer, parameter  :: iRrK9_    =  6
  integer, parameter  :: iRrK10_   =  7
  integer, parameter  :: iRrK16_   =  8
  integer, parameter  :: iRrK17_   =  9
  integer, parameter  :: iRrK18_   = 10
  integer, parameter  :: iRrK21_   = 11
  integer, parameter  :: iRrK22_   = 12
  integer, parameter  :: iRrK23_   = 13
  integer, parameter  :: iRrK24_   = 14
  integer, parameter  :: iRrK26_   = 15
  integer, parameter  :: iRrK27_   = 16
  integer, parameter  :: iRrK32_   = 17
  integer, parameter  :: iRrBeta2_ = 18
  integer, parameter  :: iRrBeta4_ = 19
  integer, parameter  :: iRrBeta6_ = 20
  integer, parameter  :: iRrBeta7_ = 21
  integer, parameter  :: iRrJ1_    = 22
  integer, parameter  :: iRrG9_    = 23
  integer, parameter  :: iRrG12_   = 24
  integer, parameter  :: iRrG13_   = 25
  integer, parameter  :: iRrG21_   = 26
  integer, parameter  :: iRrG22_   = 27
  integer, parameter  :: iRrG24_   = 28
  integer, parameter  :: iRrG25_   = 29
  integer, parameter  :: iRrG26_   = 30
  integer, parameter  :: iRrG27_   = 31
  integer, parameter  :: iRrG30_   = 32
  integer, parameter  :: iRrEA3_   = 33
  integer, parameter  :: iRrEA4_   = 34
  integer, parameter  :: iRrEA5_   = 35
  integer, parameter  :: iRrAlpha3_ = 36

  integer, parameter  :: nRrTempInd = 36

  real :: RrTempInd(nRrTempInd)

!    Temperature Depenedent

  integer, parameter  :: iRrK1_     =  1
  integer, parameter  :: iRrK2_     =  2
  integer, parameter  :: iRrK3_     =  3
  integer, parameter  :: iRrK12_    =  4
  integer, parameter  :: iRrK13_    =  5
  integer, parameter  :: iRrK19_    =  6
  integer, parameter  :: iRrK20_    =  7
  integer, parameter  :: iRrK25_    =  8
  integer, parameter  :: iRrKM12_   =  9
  integer, parameter  :: iRrA1_     = 10
  integer, parameter  :: iRrA2_     = 11
  integer, parameter  :: iRrA3_     = 12
  integer, parameter  :: iRrBeta1_  = 13
  integer, parameter  :: iRrBeta3_  = 14
  integer, parameter  :: iRrBeta5_  = 15
  integer, parameter  :: iRrBeta8_  = 16
  integer, parameter  :: iRrBeta9_  = 17
  integer, parameter  :: iRrBeta9N_ = 18
  integer, parameter  :: iRrBeta17_ = 19
  integer, parameter  :: iRrAlpha5_ = 20
  integer, parameter  :: iRrAlpha7_ = 21
  integer, parameter  :: iRrAlpha8_ = 22
  integer, parameter  :: iRrG19_    = 23

  integer, parameter  :: nRrTempDep = 23

  real :: RrTempDep(nLons, nLats, nAlts, nRrTempInd)

!    Collision Frequencies

  integer, parameter  :: iCfIN_ = 1
  integer, parameter  :: iCfEN_ = 2
  integer, parameter  :: iCfEI_ = 3

  integer, parameter  :: nCf = 3

  real :: Cf(nLons, nLats, nAlts, nRrTempInd)

contains

  subroutine set_RrTempInd

    RrTempInd(iRrK4_)    = 1.0e-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK5_)    = 4.4E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK6_)    = 4.0E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK7_)    = 2.0E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK8_)    = 1.0E-12 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK9_)    = 6.0E-11 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK10_)   = 1.3E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK16_)   = 4.8E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK17_)   = 1.0E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK18_)   = 4.0E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK21_)   = 0.047 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK22_)   = 0.171 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK23_)   = 8.E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK24_)   = 5.0E-12/ 1.0e6 ! cm-3 to m-3 
    RrTempInd(iRrK26_)   = 7.E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK27_)   = 7.7E-5 / 1.0e6 ! cm-3 to m-3

    ! Torr et al 79 ( O2+ sink) :
    RrTempInd(iRrK32_)   = 1.8e-10 / 1.0e6 ! cm-3 to m-3

    RrTempInd(iRrBETA2_) = 5.0E-12 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrBETA4_) = 5.0E-13 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrBETA6_) = 7.0E-11 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrBETA7_) = 1.06E-5 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrJ1_)    = 5.0e-9 !/ 1.0e6 ! cm-3 to m-3
    !  RrTempInd(J1)    = 0.0

    ! These are all from Reese:
    RrTempInd(iRrG9_)  = 4.8e-10 / 1.0e6
    RrTempInd(iRrG12_) = 6.0e-10 / 1.0e6
    RrTempInd(iRrG13_) = 1.0e-11 / 1.0e6
    RrTempInd(iRrG21_) = 8.0e-13 / 1.0e6
    RrTempInd(iRrG22_) = 7.0e-10 / 1.0e6
    RrTempInd(iRrG24_) = 8.0e-10 / 1.0e6
    RrTempInd(iRrG25_) = 5.2e-11 / 1.0e6
    RrTempInd(iRrG26_) = 1.3e-10 / 1.0e6
    RrTempInd(iRrG27_) = 3.0e-11 / 1.0e6
    RrTempInd(iRrG30_) = 1.0e-10 / 1.0e6
    RrTempInd(iRrEA3_) = 0.047
    RrTempInd(iRrEA4_) = 0.171
    RrTempInd(iRrEA5_) = 7.7e-5
    RrTempInd(iRrAlpha3_) = 1.0e-10 / 1.0e6

  end subroutine set_RrTempInd


end module ModRates
