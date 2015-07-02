!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModPlanet

  use ModConstants
  use ModSizeGITM, only: nAlts

  implicit none

  integer, parameter :: nSpecies = 1
  integer, parameter :: iH2_  = 1

  integer, parameter :: nSpeciesTotal = 2
  integer, parameter :: iH_   =  2

  integer, parameter  :: iH2P_   = 1
  integer, parameter  :: iHP_    = 2
  integer, parameter  :: ie_     = 3
  integer, parameter  :: nIons   = ie_
  integer, parameter  :: nIonsAdvect = 1
  integer, parameter  :: nSpeciesAll = nSpeciesTotal + nIons - 1
  
  character (len=20) :: cSpecies(nSpeciesTotal)
  character (len=20) :: cIons(nIons)

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

  integer, parameter :: nEmissions = 1
  
  real, parameter :: Gravitational_Constant = 10.0
  real, parameter :: Rotation_Period        = 100000.0
  real, parameter :: RBody                  = 1000000.0
  real, parameter :: DipoleStrength         = -10000.0e-9

  real, parameter :: OMEGABody              = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period / 3600.0
  real, parameter :: Tilt = 10.0

  real, parameter :: DaysPerYear = 365.25
  real, parameter :: SecondsPerYear = DaysPerYear * Rotation_Period

  integer, parameter :: iVernalYear   = 1999
  integer, parameter :: iVernalMonth  =    3
  integer, parameter :: iVernalDay    =   21
  integer, parameter :: iVernalHour   =    0
  integer, parameter :: iVernalMinute =    0
  integer, parameter :: iVernalSecond =    0

  ! Old orbital parameters
 !real, parameter :: SunOrbit_A = 1.000110
 !real, parameter :: SunOrbit_B = 0.034221
 !real, parameter :: SunOrbit_C = 0.001280
 !real, parameter :: SunOrbit_D = 0.000719
 !real, parameter :: SunOrbit_E = 0.000077

  !New Orbital Parameters
  !A: semi-major axis in AU
  !B: eccentricity
  !C: Longitude of perihelion
  !D: Mean Longitude
  !E: For calulating actual Longitude
 real, parameter :: SunOrbit_A = 1.00
 real, parameter :: SunOrbit_B = 0.00
 real, parameter :: SunOrbit_C = 0.00
 real, parameter :: SunOrbit_D = 0.00
 real, parameter :: SunOrbit_E = 129597740.63

  !Used as a damping term in Vertical solver.
  real :: VertTau(nAlts)

  logical :: IsEarth = .false.
  logical :: IsMars = .false.
  logical :: IsTitan = .false.
  logical :: NonMagnetic = .false.
  real, parameter :: PlanetNum = 0.0426

  character (len=10) :: cPlanet = "LV-426"
  
  integer, parameter :: nEmissionWavelengths = 20
  integer, parameter :: nPhotoBins = 190

  ! These are for the neutral friction routine...

  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2

 real, parameter, dimension(2, 2) :: Diff0 = 1.0e4 * reshape( (/ &
       ! 0      02     N2      N     NO
       !---------------------------------+
       0.00,   0.260,&            ! O
       0.26,   0.000/), (/2,2/) )  ! N

  ! These are the exponents
  real, parameter, dimension(2, 2) :: DiffExp = reshape( (/ &
       ! 0      02     N2
       !---------------------------------+
       0.00,  0.75, &             ! H
       0.75,  0.00/), (/2,2/) )  ! H2

contains

  subroutine init_planet

    use ModTime

    integer :: itime(7)

    Mass(iH_)       = 1.0 * AMU
    Mass(iH2_)      = 2.0 * AMU

    cSpecies(iH_)   = "H"
    cSpecies(iH2_)  = "H2"

    cIons(iHP_)     = "H!U+!N"
    cIons(iH2P_)    = "H2!U+!N"
    cIons(ie_)      = "e-"

    Vibration(iH_)    = 5.0
    if (nSpecies > 1) Vibration(iH2_)   = 7.0

    MassI(iHP_) = Mass(iH_)
    MassI(iH2P_) = Mass(iH2_)
    MassI(ie_) = Mass_Electron

    VertTau = 1.0e9

    itime = 0
    itime(1) = iVernalYear
    itime(2) = iVernalMonth
    itime(3) = iVernalDay
    itime(4) = iVernalHour
    itime(5) = iVernalMinute
    itime(6) = iVernalSecond
    call time_int_to_real(itime, VernalTime)

  end subroutine init_planet

!! Placeholder subroutines (for Titan specific Phyisics)

  subroutine init_radcooling
  return
  end subroutine init_radcooling

  subroutine init_magheat
  return
  end subroutine init_magheat

  subroutine init_isochem
  return
  end subroutine init_isochem

  subroutine init_aerosol
  return
  end subroutine init_aerosol

end module ModPlanet
