!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModPlanet

  use ModConstants
  use ModSizeGITM, only: nAlts

  implicit none

! Modified (01/18/07) : SWB :   Aij, s-exponents for mutual diffusion
! Majors (4):  COntrol the Pressures Gradients and winds
  integer, parameter :: nSpecies = 3
  integer, parameter :: iH2_     = 1
  integer, parameter :: iHe_     = 2
  integer, parameter :: iH_      = 3

! Minors (6) : Ride on background of mean winds derived from majors
  integer, parameter :: iCH3_  =  4
  integer, parameter :: iCH4_  =  5
  integer, parameter :: iC2H2_ =  6
  integer, parameter :: iC2H3_ =  7
  integer, parameter :: iC2H4_ =  8
  integer, parameter :: iC2H6_ =  9
  integer, parameter :: nSpeciesTotal = iC2H6_

! Major Ions (4):  Most Important to MWACM code
  integer, parameter  :: iH3P_    = 1
  integer, parameter  :: iHP_     = 2
  integer, parameter  :: iH2P_    = 3
  integer, parameter  :: iCH5P_   = 4
  integer, parameter  :: iHC2H5P_ = 5
  integer, parameter  :: iHeHP_   = 6
  integer, parameter  :: ie_      = 7
  integer, parameter  :: nIons    = ie_
  integer, parameter  :: nIonsAdvect = 2

  character (len=20) :: cSpecies(nSpeciesTotal)
  character (len=20) :: cIons(nIons)

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

  ! When you want to program in emissions, you can use this...
  integer, parameter :: nEmissions = 10

  !! CHANGE
  real, parameter :: GC_Mars                = 24.5                    ! m/s^2
  real, parameter :: RP_Mars                = 35730.0                 ! seconds
  real, parameter :: R_Mars                 = 71400.0*1000.0          ! meters
  !! CHANGE
  real, parameter :: DP_Mars                = 311000.0e-9             ! nT

  real, parameter :: Gravitational_Constant = GC_Mars
  real, parameter :: Rotation_Period        = RP_Mars
  real, parameter :: RBody                  = R_Mars
  real, parameter :: DipoleStrength         = DP_Mars

  real, parameter :: OMEGABody              = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period / 3600.0
  !! CHANGE
  real, parameter :: Tilt = 11.0

  ! This is the Vernal Equinox at Midnight (Ls = 0!!!)
  ! Earth-Mars clocks are set from this epoch at vernal equinox
  integer, parameter :: iVernalYear   = 1998
  integer, parameter :: iVernalMonth  =    7
  integer, parameter :: iVernalDay    =   14
  integer, parameter :: iVernalHour   =   16
  integer, parameter :: iVernalMinute =    0
  integer, parameter :: iVernalSecond =    0

  real, parameter :: SunOrbit_A = 1.52
  real, parameter :: SunOrbit_B = 0.04
  real, parameter :: SunOrbit_C = 0.15
  real, parameter :: SunOrbit_D = 0.00
  real, parameter :: SunOrbit_E = 0.00

  real, parameter :: DaysPerYear = 670.0
  real, parameter :: SecondsPerYear = DaysPerYear * Rotation_Period

  !Used as a damping term in Vertical solver.
  real, dimension(nAlts) :: VertTau = 1.0e9 

  logical :: IsEarth = .false.
  character (len=10) :: cPlanet = "Mars"

!  real :: KappaTemp0 = 2.22e-4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are Modified for Mars by SWB: 1/18/07
  ! -- Most source state Dij = Dji (check)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! These are Aij coefficients from B&K (1973) formulation: Aij*1.0E+17
  ! Use Jared Bell's Titan GITM formulation for Mars GITM

  real, parameter, dimension(4, 4) :: Diff0 = 1.0e17 * reshape( (/ &
! integer, parameter :: iCO2_    = 1
! integer, parameter :: iCO_     = 2
! integer, parameter :: iO_      = 3
! integer, parameter :: iN2_     = 4
     !---------------------------------+
     ! i=C02      CO      O     N2  
     !---------------------------------+
!       0.0000, 0.7762, 0.2119,  0.6580, &            ! CO2
!       0.4940, 0.0000, 0.9466,  0.9280, &            ! CO
!       0.7703, 0.9481, 0.0000,  0.9690, &            ! O
!       0.6580, 0.9280, 0.9690,  0.0000 /), (/4,4/) ) ! N2
!
       0.0000, 0.7762, 0.2119,  0.6580, &            ! CO2
       0.7762, 0.0000, 0.9466,  0.9280, &            ! CO
       0.2219, 0.9466, 0.0000,  0.9690, &            ! O
       0.6580, 0.9280, 0.9690,  0.0000 /), (/4,4/) ) ! N2

  ! These are s-exponents from B&K (1973) formulation: T**s
  real, parameter, dimension(4, 4) :: DiffExp = reshape( (/ &
       !---------------------------------+
       !i= CO2    CO      O      N2
       !---------------------------------+
       0.000,  0.750,  0.750, 0.752, &             ! CO2
       0.750,  0.000,  0.750, 0.710, &             ! CO
       0.750,  0.750,  0.000, 0.774, &             ! O
       0.752,  0.710,  0.774, 0.000 /), (/4,4/) )  ! N2

!     Arrays filled in init_radcool in data statements (np = 68)
  integer, parameter :: np=68
  real,dimension(np) :: pnbr,ef1,ef2,co2vmr,o3pvmr,n2covmr

  !! Stuff for initial conditions

  real , Dimension(-1:nAlts + 2) :: newalt
  real , Dimension(-1:nAlts + 2) :: InTemp
  real , Dimension(-1:nAlts + 2) :: IneTemp
  real , Dimension(-1:nAlts + 2) :: InITemp
  real , Dimension(-1:nAlts + 2,nSpeciesTotal) :: InNDensityS 
  real , Dimension(-1:nAlts + 2,nIons) :: InIDensityS

contains

  subroutine init_planet

    use ModTime
    use ModIoUnit, only : UnitTmp_

    integer :: iTime(7), iiAlt

!   Mass = AMU * mean molecular weight  (unlike TGCM codes)

    Mass(iH_)    = AMU
    Mass(iH2_)   = 2.0 * AMU
    Mass(iHe_)   = 4.0 * AMU
    Mass(iCH3_)  = 12.0 * AMU + 3*Mass(iH_)
    Mass(iCH4_)  = 12.0 * AMU + 4*Mass(iH_)
    Mass(iC2H2_) = 2 * 12.0 * AMU + 2*Mass(iH_)
    Mass(iC2H3_) = 2 * 12.0 * AMU + 3*Mass(iH_)
    Mass(iC2H4_) = 2 * 12.0 * AMU + 4*Mass(iH_)
    Mass(iC2H6_) = 2 * 12.0 * AMU + 6*Mass(iH_)

    cSpecies(iH_)    = "H"
    cSpecies(iH2_)   = "H!D2!N"
    cSpecies(iHe_)   = "He"
    cSpecies(iCH3_)  = "CH!D3!N"
    cSpecies(iCH4_)  = "CH!D4!N"
    cSpecies(iC2H3_) = "C!D2!NH!D2!N"
    cSpecies(iC2H4_) = "C!D2!NH!D4!N"
    cSpecies(iC2H6_) = "C!D2!NH!D6!N"

    cIons(iH3P_)   = "CH!D3!U+!N"
    cIons(iHP_)    = "H!U+!N"
    cIons(iH2P_)   = "H!D2!U+!N"
    cIons(iCH5P_)  = "CH!D5!U+!N"
    cIons(iC2H5P_) = "C!D2!NH!D5!U+!N"
    cIons(iHeHP_)  = "HeH!U+!N"
    cIons(ie_)     = "e-"

    !! CHANGE
    Vibration(iH_)    = 5.0
    Vibration(iH2_)   = 7.0
    Vibration(iHe_)   = 5.0
    Vibration(iCH4_)  = 7.0
    Vibration(iC2H4_) = 7.0

    MassI(iH3P_)   = 3*Mass(iH_)
    MassI(iHP_)    = Mass(iH_)
    MassI(iH2P_)   = 2*Mass(iH_)
    MassI(iCH5P_)  = Mass(iCH4_) + Mass(iH_)
    MassI(iC2H5P_) = Mass(iC2H4_) + Mass(iH_)
    MassI(iHeHP_)  = Mass(iH_) + Mass(iHe_)
    MassI(ie_)     = Mass_Electron

    itime = 0
    itime(1) = iVernalYear
    itime(2) = iVernalMonth
    itime(3) = iVernalDay
    itime(4) = iVernalHour
    itime(5) = iVernalMinute
    itime(6) = iVernalSecond
    call time_int_to_real(itime, VernalTime)

  end subroutine init_planet

end module ModPlanet
