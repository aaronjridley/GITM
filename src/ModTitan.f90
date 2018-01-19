!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModPlanet

  use ModConstants
  use ModOrbital
  use ModSizeGITM

  implicit none

  integer, parameter :: iN2_      = 1
  integer, parameter :: iCH4_     = 2
  integer, parameter :: iAr_      = 3
  integer, parameter :: iHCN_     = 4
  integer, parameter :: iC2H4_    = 5
  integer, parameter :: i15N2_    = 6
  integer, parameter :: i13CH4_   = 7
  integer, parameter :: iH_       = 8
  integer, parameter :: iH2_      = 9
  integer, parameter :: iH2CN_    = 10
  integer, parameter :: iN4S_     = 11
  integer, parameter :: nSpecies  = 11

! C-Compounds
  integer, parameter :: iCH3_   =  12
  integer, parameter :: i3CH2_  =  13
  integer, parameter :: i1CH2_  =  14
  integer, parameter :: iCH_    =  15
  integer, parameter :: iN2D_   =  16
  integer, parameter :: nSpeciesTotal = iN2D_

!  integer, parameter :: iC_     =  16
! C2-Compounds
!  integer, parameter :: iC2_    =  17
!  integer, parameter :: iC2H_   =  18
!  integer, parameter :: iC2H2_  =  19
!  integer, parameter :: iC2H3_  =  20
!  integer, parameter :: iC2H5_  =  21
!  integer, parameter :: iC2H6_  =  22
! C3-Compounds
!  integer, parameter :: iC3H8_  =  23
! C4-Compounds
!  integer, parameter :: iC4H2_  =  24
! Nitrogen Compounds
!  integer, parameter :: iCN_    =  25
!  integer, parameter :: iC2N_   =  26
!  integer, parameter :: iC3N_   =  27
!
!  integer, parameter :: iC2N2_  =  28
!  integer, parameter :: iC4N2_  =  29
!
!  integer, parameter :: iHC3N_  =  30
!
!  integer, parameter :: iNH_    =  31
!  integer, parameter :: iN2D_   =  32
!  integer, parameter :: iN_2P_  =  33

!  integer, parameter :: nSpeciesTotal = iN_2P_

! Major Ions (4):  Most Important to MWACM code

! Major Ions
   integer, parameter  :: iHCNHP_ = 1
   integer, parameter  :: iC2H5P_ = 2
   integer, parameter  :: iCH3P_  = 3
   integer, parameter  :: iCH4P_  = 4
   integer, parameter  :: iNP_    = 5
   integer, parameter  :: iN2P_   = 6
   integer, parameter  :: iCH5P_  = 7
   integer, parameter  :: ie_     = 7

! Methane-type Ions
!   integer, parameter  :: iCH2P_  = 4
!   integer, parameter  :: iCH3P_  = 5
!   integer, parameter  :: iCH4P_  = 6
!   integer, parameter  :: iCHP_   = 7
!   integer, parameter  :: iCP_    = 8
! Nitrogen Ions
!
!   integer, parameter  :: iNP_    = 9
!   integer, parameter  :: iN2P_   = 10
!
!   integer, parameter  :: iNHP_   = 11
!   integer, parameter  :: iN2HP_  = 12
!   integer, parameter  :: iHCNP_  = 13
!
!   integer, parameter  :: iC2HP_  = 14
!   integer, parameter  :: iC2H2P_ = 15
!   integer, parameter  :: iC2H3P_ = 16
!   integer, parameter  :: iC2H4P_ = 17
!
!   integer, parameter  :: iH3P_   = 18
!   integer, parameter  :: iH2P_   = 19
!   integer, parameter  :: iHP_    = 20
!   integer, parameter  :: ie_     = 21

   integer, parameter  :: nIons   = ie_
   integer, parameter  :: nIonsAdvect = 2
   integer, parameter  :: nSpeciesAll = nSpeciesTotal + nIons - 1

  character (len=20) :: cSpecies(nSpeciesTotal)
  character (len=20) :: cIons(nIons)

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

  integer, parameter :: iE2470_ = 1
  integer, parameter :: iE7320_ = 2
  integer, parameter :: iE3726_ = 3
  integer, parameter :: iE5200_ = 4
  integer, parameter :: iE10400_ = 5
  integer, parameter :: iE6300_ = 6
  integer, parameter :: iE6364_ = 7

  integer, parameter :: nEmissions = 10

  integer, parameter :: i3371_ = 1
  integer, parameter :: i4278_ = 2
  integer, parameter :: i5200_ = 3
  integer, parameter :: i5577_ = 4
  integer, parameter :: i6300_ = 5
  integer, parameter :: i7320_ = 6
  integer, parameter :: i10400_ = 7
  integer, parameter :: i3466_ = 8
  integer, parameter :: i7774_ = 9
  integer, parameter :: i8446_ = 10
  integer, parameter :: i3726_ = 11

  real, parameter :: GC_Titan               = 1.354                  ! m/s^2
  real, parameter :: RP_Titan               = 1377684.0              ! seconds
  real, parameter :: R_Titan                = 2575.0*1000.0          ! meters
  real, parameter :: DP_Titan               = 1.0e-32                ! nT

  real, parameter :: Gravitational_Constant = GC_Titan
  real, parameter :: Rotation_Period        = RP_Titan
  real, parameter :: RBody                  = R_Titan
  real, parameter :: DipoleStrength         = DP_Titan

  real, parameter :: OMEGABody              = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period / 3600.0
  real, parameter :: Tilt = 26.70

  real, parameter :: DaysPerYear = 10759.53  ! Saturn has 29.5 Earth years
  real, parameter :: SecondsPerYear = DaysPerYear * Rotation_Period

  integer, parameter :: iVernalYear   = 1980
  integer, parameter :: iVernalMonth  =    2
  integer, parameter :: iVernalDay    =   14
  integer, parameter :: iVernalHour   =   12
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
 real, parameter :: SunOrbit_A = 9.5400001124
 real, parameter :: SunOrbit_B = 0.04
 real, parameter :: SunOrbit_C = 0.15
 real, parameter :: SunOrbit_D = 0.00
 real, parameter :: SunOrbit_E = 0.00

 real :: semimajoraxis_0 = semimajor_Titan
 real :: eccentricity_0 = eccentricity_Titan
 real :: inclination_0 = inclination_Titan
 real :: longitudePerihelion_0 = longitudePerihelion_Titan
 real :: longitudeNode_0 = longitudeNode_Titan
 real :: meanLongitude_0 = meanLongitude_Titan
 real :: semimajordot = semimajordot_Titan
 real :: eccentricitydot = eccentricitydot_Titan
 real :: inclinationdot = inclinationdot_Titan
 real :: longitudePeriheliondot = longitudePeriheliondot_Titan
 real :: longitudeNodedot = longitudeNodedot_Titan
 real :: meanLongitudedot = meanLongitudedot_Titan

  !Used as a damping term in Vertical solver.
  real :: VertTau(nAlts)

  logical :: IsEarth = .false.
  logical :: IsMars = .false.
  logical :: IsTitan = .true.
  logical :: NonMagnetic = .true.
  real, parameter :: PlanetNum = 0.061

  character (len=10) :: cPlanet = "Titan"

  integer, parameter :: nEmissionWavelengths = 20
  integer, parameter :: nPhotoBins = 190


  ! These are for the neutral friction routine...

  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
  ! JMB: Updated the Diff0 and DiffExp to be dynamic use (nspecies, nspecies)
  ! in the shape of the array

!  real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e16 * reshape( (/ &
!       !-----------------------
!       ! i=N2  CH4     Ar
!       !-----------------------
!       0.0000, 7.3000, 6.6200,   & ! N2
!       7.3000, 0.0000, 5.7400,   & ! CH4
!       6.6200, 5.7400, 0.0000/), & ! Ar
!       (/nSpecies,nSpecies/) )
!
!  real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
!       !----------------------
!       ! N2      CH4     Ar
!       !----------------------
!       0.0000, 0.7500, 0.7520,   & ! N2
!       0.7500, 0.0000, 0.7850,   & ! CH4
!       0.7520, 0.7850, 0.0000/), & ! Ar
!       (/nSpecies,nSpecies/) )


  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
  ! JMB: Updated the Diff0 and DiffExp to be dynamic use (nspecies, nspecies)
  ! in the shape of the array
!  real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e16 * reshape( (/ &
!  !----------------------------------------------------
!  ! i=N2  CH4     Ar      HCN     C2H4   15N2   13CH4
!  !----------------------------------------------------
!  0.0000, 7.3000, 6.6200, 5.1834, 5.090, 5.000, 7.300,& ! N2
!  7.3000, 0.0000, 5.7400, 5.1132, 5.080, 7.300, 5.645,&! CH4
!  6.6200, 5.7400, 0.0000, 6.5184, 5.090, 4.605, 4.623,& ! Ar
!  5.1834, 5.1132, 6.5184, 0.0000, 5.090, 5.139, 4.944,& ! HCN
!  5.0900, 5.0800, 5.0900, 5.0900, 0.000, 5.000, 7.300,& ! C2H4
!  5.0000, 7.3000, 4.6050, 5.1390, 5.000, 0.000, 7.300,& ! 15N2
!  7.3000, 5.6450, 4.6230, 4.9440, 7.300, 7.300, 0.000/),& ! 13CH4
!  (/nSpecies,nSpecies/) )
!!
!  real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
!  !---------------------------------------------------
!  ! i=N2  CH4     Ar      HCN      C2H4   15N2   13CH4
!  !---------------------------------------------------
!  0.0000, 0.7500, 0.7520, 0.8100, 0.810, 0.750, 0.750,&  ! N2
!  0.7500, 0.0000, 0.7850, 0.7650, 0.765, 0.750, 0.750,&  ! CH4
!  0.7520, 0.7850, 0.0000, 0.7520, 0.750, 0.750, 0.750,&  ! Ar
!  0.8100, 0.7650, 0.7520, 0.0000, 0.750, 0.750, 0.750,&  ! HCN
!  0.8100, 0.7650, 0.7500, 0.7500, 0.000, 0.750, 0.750,&  ! C2H4
!  0.8100, 0.7650, 0.7500, 0.7500, 0.750, 0.000, 0.750,&  ! 15N2
!  0.8100, 0.7650, 0.7500, 0.7500, 0.750, 0.750, 0.000/),&  ! 13CH4
! (/nSpecies,nSpecies/) )


  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
  ! JMB: Updated the Diff0 and DiffExp to be dynamic use (nspecies, nspecies)
  ! in the shape of the array
  real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e16 * reshape( (/ &
  !----------------------------------------------------------------------------------
  ! i=N2  CH4     Ar      HCN     C2H4   15N2   13CH4   H      H2    H2CN   N(4S)
  !----------------------------------------------------------------------------------
  0.0000, 7.3000, 6.6200, 5.1834, 5.090, 5.000, 7.300,48.700,33.700, 5.09, 9.690,& ! N2
  7.3000, 0.0000, 5.7400, 5.1132, 5.080, 7.300, 5.645,22.900,22.920, 5.09, 6.125,&! CH4
  6.6200, 5.7400, 0.0000, 6.5184, 5.090, 4.605, 4.623,48.700,15.500, 5.09, 9.690,& ! Ar
  5.1834, 5.1132, 6.5184, 0.0000, 5.090, 5.139, 4.944,48.700,17.451, 5.09, 9.690,& ! HCN
  5.0900, 5.0800, 5.0900, 5.0900, 0.000, 5.000, 7.300,48.700,33.700, 5.09, 9.690,& ! C2H4
  5.0000, 7.3000, 4.6050, 5.1390, 5.000, 0.000, 7.300,48.700,33.700, 5.09, 9.690,& ! 15N2
  7.3000, 5.6450, 4.6230, 4.9440, 7.300, 7.300, 0.000,22.920,22.900, 5.09, 6.125,& ! 13CH4
  48.700,22.9000,48.7000,48.7000,48.700,48.700,22.920, 0.000,48.700,48.70, 48.70,& ! H
  33.700,22.9200,15.5000,17.4510,33.700,33.700,22.900,48.700, 0.000,33.70, 33.70,& ! H2
  5.0900, 5.0800, 5.0900, 5.0900, 5.090, 5.090, 5.090,48.700,33.700, 0.00,  9.69,& ! H2CN
  9.6900, 6.1250, 9.6900, 9.6900, 9.690, 9.690, 6.125,48.700,33.700, 9.69, 0.000/), & ! N(4S)
  (/nSpecies,nSpecies/) )
!!
  real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
  !---------------------------------------------------------------------------------
  ! i=N2  CH4     Ar      HCN      C2H4   15N2   13CH4  H      H2     H2CN   N(4S)
  !---------------------------------------------------------------------------------
  0.0000, 0.7500, 0.7520, 0.8100, 0.810, 0.750, 0.750, 0.698, 0.710, 0.810, 0.764,&  ! N2
  0.7500, 0.0000, 0.7850, 0.7650, 0.765, 0.750, 0.750, 0.750, 0.765, 0.765, 0.750,&  ! CH4
  0.7520, 0.7850, 0.0000, 0.7520, 0.750, 0.750, 0.750, 0.750, 0.850, 0.750, 0.750,&  ! Ar
  0.8100, 0.7650, 0.7520, 0.0000, 0.750, 0.750, 0.750, 0.750, 0.750, 0.750, 0.750,&  ! HCN
  0.8100, 0.7650, 0.7500, 0.7500, 0.000, 0.750, 0.750, 0.750, 0.750, 0.750, 0.750,&  ! C2H4
  0.8100, 0.7650, 0.7500, 0.7500, 0.750, 0.000, 0.750, 0.750, 0.750, 0.750, 0.750,&  ! 15N2
  0.8100, 0.7650, 0.7500, 0.7500, 0.750, 0.750, 0.000, 0.750, 0.750, 0.750, 0.764,&  ! 13CH4
  0.6980, 0.7500, 0.7500, 0.7500, 0.750, 0.750, 0.750, 0.000, 0.750, 0.750, 0.750,&  ! H
  0.7100, 0.7650, 0.8500, 0.7500, 0.750, 0.750, 0.750, 0.750, 0.000, 0.750, 0.750,&  ! H2
  0.8100, 0.7650, 0.7500, 0.7500, 0.750, 0.750, 0.750, 0.750, 0.750, 0.000, 0.750,&   ! C2H4
  0.7640, 0.7500, 0.7500, 0.7500, 0.750, 0.764, 0.750, 0.750, 0.750, 0.750, 0.000 /), &! N4S
 (/nSpecies,nSpecies/) )


! Titan-Specific Parameters
   integer, parameter :: rotlines = 60       ! # of HCN Rotational Lines
   integer, parameter :: newrotfreqs = 9     ! # of Frequency Gaussian Points
   integer, parameter :: rotfreqs = 8        ! # of Frequency Gaussian Points
   integer, parameter :: newrotpts   = 17    ! # of Frequency Gaussian Points
   integer, parameter :: rotpts   = 16       ! # of Altitude Gassuian Points

   real, dimension(rotlines) :: freqw        ! Rotational Lines in Wave Numbers (cm^-1)
   real, dimension(rotlines) :: freqhz       ! Rotational Lines (in Hz)
   real, dimension(rotlines) :: ALPHAL0, ALPHAS0     !! Lorentz Halfwidth (in cm^-1)
   real, dimension(rotlines) :: ALPHALExp            !! Lorentz Halfwidth (in cm^-1)
   real, dimension(rotlines) :: STref, Eref          !! Lorentz Halfwidth (in cm^-1)
   real, dimension(rotlines) :: Einstein             !! Lorentz Halfwidth (in cm^-1)
   real, dimension(rotlines,9) :: pmat               !! Intensity Interpolation Matrix
   real, parameter :: atmtopas = 1.0/(1.01295e+05)
   integer, parameter :: nQuadPointsAlts    = 17 !! # of Altitude Gassuian Points
   integer, parameter :: nQuadPointsFreqs   = 5 !! # of Altitude Gassuian Points
! Gaussian Quadrature Points
   real, dimension(8,4) :: Qd, Qf1
! Guassian Quadrature Points
   real, dimension( nQuadPointsAlts) :: AltQuadAbscissas, AltQuadWeights
   real, dimension( nQuadPointsFreqs) :: FreqQuadAbscissas, FreqQuadWeights
! 86400 24 hours
! 43200 12 hours
! 28800 8  hours
! 21600 6  hours
   real :: dTCooling = 28800.00
   real :: dTRatio   = 3600.00

!! Read in Variables for Titan Startup
   real , Dimension(-1:nAlts + 2) :: newalt
   real , Dimension(-1:nAlts + 2) :: junk
   real , Dimension(-1:nAlts + 2) :: InTemp
   real , Dimension(-1:nAlts + 2) :: IneTemp
   real , Dimension(-1:nAlts + 2) :: InITemp
   real , Dimension(-1:nAlts + 2,nSpeciesTotal) :: InNDensityS
   real , Dimension(-1:nAlts + 2,nIons) :: InIDensityS
   real , Dimension(-1:nAlts + 2) :: InOP_Heating
   real , Dimension(-1:nAlts + 2) :: InHP_Heating
   real , Dimension(-1:nAlts + 2) :: InPickup_Heating
   real , Dimension(1:nLons,1:nLats,-1:nAlts+2,nBlocksMax) :: OP_Heating3D
   real , Dimension(1:nLons,1:nLats,-1:nAlts+2,nBlocksMax) :: HP_Heating3D
   real , Dimension(1:nLons,1:nLats,-1:nAlts+2,nBlocksMax) :: Pickup_Heating3D

!! Aerosol Calculation Variables
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: AerosolClassA
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: AerosolClassB
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: AerosolClassC

!! Aerosol Loss Rates
   real , Dimension(1:nLons,1:nLats,1:nAlts,nSpeciesTotal,nBlocksMax) :: AerosolTrappingLoss
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: H2AerosolProduction


!! Isotope Chemistry Variables
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: IsotopeScaling

!!! JMB eTemp Vars

   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: PhotoEHeat
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Collision_ee
!
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Collision_ch5pe
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Collision_c2h5pe
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Collision_hcnhpe
!
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Collision_n2e
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Collision_ch4e
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Collision_h2e

   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Collision_Totale
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: eLambda_Temp
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Qei
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Qen_elas
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Qen_N2rot
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Qen_N2vib
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Qen_inelas

!!! Time independent part and Time-dep part -- due to eccentricity
   real , Dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:3,nBlocksMax) :: TideAIndp
   real , Dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:3,nBlocksMax) :: TideADep

   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: HCNCoolRatio
   real , Dimension(1:nLons,1:nLats,1:nAlts,1:nSpecies,nBlocksMax) :: MeshRatio
   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: MeanMeshRatio

   real , Dimension(1:nLons,1:nLats,1:nAlts,1:nSpecies,nBlocksMax) :: YChar, XpS
   real , Dimension(1:nLons,1:nLats,1:nAlts,1:nSpecies,nBlocksMax) :: SmithErf, NewErf
   real , Dimension(1:nLons,1:nLats,1:nAlts,1:nSpecies,nBlocksMax) :: OldChap, NewChap
   real , Dimension(1:nLons,1:nLats,1:nAlts,1:nSpecies,nBlocksMax) :: OldHs, NewHs
   real , Dimension(1:nLons,1:nLats,1:nAlts,1:nSpecies,nBlocksMax) :: ChapIntegralS
contains

  subroutine init_planet

    use ModTime

    integer :: itime(7)

    Mass(iN2_)    = 14.00674 * AMU * 2.0
    Mass(iCH4_)   =  1.00790 * AMU * 4 + 12.011 * AMU
    Mass(iAr_)    = 40.00000 * AMU
    Mass(iHCN_)   =  1.00790 * AMU + 14.00674 * AMU + 12.011 * AMU
    Mass(iC2H4_)  =  1.00790 * AMU * 4.0 + 12.011 * AMU * 2
    Mass(i15N2_)  = 14.00674 * AMU + 15.00674 * AMU
    Mass(i13CH4_) = 13.01100 * AMU + 1.0079 * AMU * 4
    Mass(iH_)     =  1.00790 * AMU
    Mass(iH2_)    =  1.00790 * AMU * 2.0
    Mass(iH2CN_)  =  1.00790 * AMU * 2.0 + 14.00674 * AMU + 12.011 * AMU
    Mass(iN4S_)   = 14.00674 * AMU

    Mass(iN2D_)  =  Mass(iN4S_)

    Mass(i3CH2_)    =  1.0079 * AMU * 2 + 12.011 * AMU
    Mass(i1CH2_)    =  1.0079 * AMU * 2 + 12.011 * AMU
    Mass(iCH3_)     =  1.0079 * AMU * 3 + 12.011 * AMU
    Mass(iCH_)      =  1.0079 * AMU + 12.011 * AMU

!    Mass(iN_2P_)  =  Mass(iN4S_)
!    Mass(iC_)       = 12.011 * AMU
!
!! C2 Compounds
!    Mass(iC2_)      =                     12.011 * AMU * 2
!    Mass(iC2H_)     =  1.0079 * AMU * 1 + 12.011 * AMU * 2
!    Mass(iC2H2_)    =  1.0079 * AMU * 2 + 12.011 * AMU * 2
!    Mass(iC2H3_)    =  1.0079 * AMU * 3 + 12.011 * AMU * 2
!    Mass(iC2H5_)    =  1.0079 * AMU * 5 + 12.011 * AMU * 2
!    Mass(iC2H6_)    =  1.0079 * AMU * 6 + 12.011 * AMU * 2
!! C3
!    Mass(iC3H8_)    =  1.0079 * AMU * 8 + 12.011 * AMU * 3
!! C4
!    Mass(iC4H2_)    =  1.0079 * AMU * 2 + 12.011 * AMU * 4
!! C-N
!    Mass(iCN_)      =  12.011 * AMU  + 14.00674 * AMU
!    Mass(iC2N_)     =  12.011*AMU*2  + 14.00674 * AMU
!    Mass(iC3N_)     =  12.011*AMU*3  + 14.00674 * AMU
! C2-N
!    Mass(iC2N2_)    =  12.011*AMU*2  + 14.00674 * AMU*2
! Cx-N
!    Mass(iC4N2_)    =  12.011*AMU*4  + 14.00674 * AMU*2
!    Mass(iHC3N_)    =  1.0079*AMU*1 + 12.011*AMU*3  + 14.00674 * AMU*1

! Nitrogen Compounds
!    Mass(iN4S_)     =  14.00674 * AMU
!    Mass(iNH_)      =  1.0079 * AMU + 14.00674 * AMU



!! Majors
    cSpecies(iN2_)    = "N2"
    cSpecies(iCH4_)   = "CH4"
    cSpecies(iAr_)    = "Ar"
    cSpecies(iHCN_)   = "HCN"
    cSpecies(iC2H4_)  = "C2H4"
    cSpecies(i15N2_)  = "15N2"
    cSpecies(i13CH4_) = "13CH4"
    cSpecies(iH_)     = "H"
    cSpecies(iH2_)    = "H2"
    cSpecies(iH2CN_)  = "H2CN"
    cSpecies(iN4S_ )  = "N(4S)"

!!!!! MINORS
!! C-Compounds
    cSpecies(iCH3_)   = "CH3"
    cSpecies(i3CH2_)  = "3CH2"
    cSpecies(i1CH2_)  = "1CH2"
    cSpecies(iCH_)    = "CH"
    cSpecies(iN2D_)  = "N(2D)"
!    cSpecies(iC_ )    = "C"
!! C2-Compounds
!    cSpecies(iC2_)   = "C2"
!    cSpecies(iC2H_)  = "C2H"
!    cSpecies(iC2H2_) = "C2H2"
!    cSpecies(iC2H3_) = "C2H3"
!    cSpecies(iC2H5_) = "C2H5"
!    cSpecies(iC2H6_) = "C2H6"
!! C3-Compounds
!    cSpecies(iC3H8_) = "C3H8"
!! C4-Compounds
!    cSpecies(iC4H2_) = "C4H2"
!! C-N Compounds
!    cSpecies(iCN_)   = "CN"
!    cSpecies(iC2N_)  = "C2N"
!    cSpecies(iC3N_)  = "C3N"
!    cSpecies(iC2N2_) = "C2N2"
!    cSpecies(iC4N2_) = "C4N2"
! HCN-compounds
!    cSpecies(iHC3N_) = "HC3N"
!    cSpecies(iNH_)   = "NH"
!    cSpecies(iN_2P_) = "N(2P)"

!! Major Ions
    cIons(iHCNHP_)  = "HCNH+"
    cIons(iC2H5P_)  = "C2H5+"
    cIons(iCH5P_ )  = " CH5+"
! Minor Ions
    cIons(iCH3P_)   = "CH3+"
    cIons(iCH4P_)   = "CH4+"
    cIons(iNP_)     = "N+"
    cIons(iN2P_)    = "N2+"

!    cIons(iCH2P_)   = "CH2+"
!    cIons(iCHP_ )   = "CH+"
!    cIons(iCP_  )   = "C+"
!
!    cIons(iNHP_)    = "NH+"
!    cIons(iN2HP_)   = "N2H+"
!    cIons(iHCNP_)   = "HCN+"
!
!    cIons(iC2HP_)   = "C2H+"
!    cIons(iC2H2P_)  = "C2H2+"
!    cIons(iC2H3P_)  = "C2H3+"
!    cIons(iC2H4P_)  = "C2H4+"
!
!    cIons(iH3P_)    = "H3+"
!    cIons(iH2P_)    = "H2+"
!    cIons(iHP_)     = "H+"
    cIons(ie_)      = "e-"

!! Used for CP Calculations.
    Vibration(iN2_)    = 5.0 + 2.0
    Vibration(iH2_)    = 5.0 + 2.0
    Vibration(iCH4_)   = 6.5374 + 2.0
    Vibration(iHCN_)   = 5.0 + 2.0
    Vibration(iAr_)    = 5.0
    Vibration(i13CH4_) = 6.5374 + 2.0
    Vibration(i15N2_)  = 5.0 + 2.0
    Vibration(iH_)     = 5.0

!! Ions
    MassI(iHCNHP_)  =  1.0079 * AMU * 2 + 14.00674 * AMU + 12.011 * AMU
    MassI(iC2H5P_)  =  1.0079 * AMU * 5 + 12.011 * AMU * 2
    MassI(iCH5P_)   =  1.0079 * AMU * 5 + 12.011 * AMU
!
!    MassI(iCH2P_ )  =  1.0079 * AMU * 2 + 12.011 * AMU
    MassI(iCH3P_)   =  1.0079 * AMU * 3 + 12.011 * AMU
    MassI(iCH4P_)   =  1.0079 * AMU * 4 + 12.011 * AMU
    MassI(iN2P_)    =  14.00674 * AMU * 2
    MassI(iNP_)     =  14.00674 * AMU
!
!    MassI(iC2H2P_)   =  1.0079 * AMU * 2 + 12.011 * AMU*2
!    MassI(iC2H3P_)   =  1.0079 * AMU * 3 + 12.011 * AMU*2
!    MassI(iC2H4P_)   =  1.0079 * AMU * 4 + 12.011 * AMU*2
!
!    MassI(iNHP_)     =  14.00674 * AMU + 1.0079 * AMU
!    MassI(iN2HP_)    =  14.00674 * AMU*2 + 1.0079 * AMU
!    MassI(iHCNP_)    =  14.00674 * AMU + 1.0079 * AMU + 12.011 * AMU
!
!
!    MassI(iH3P_)    =  1.0079 * AMU * 3
!    MassI(iH2P_)    =  1.0079 * AMU * 2
!    MassI(iHP_ )    =  1.0079 * AMU
!    MassI(ie_) = Mass_Electron
!
!    MassI(iCHP_)    =  1.0079 * AMU + 12.011 * AMU
!    MassI(iCP_ )    =  12.011 * AMU

    MassI(ie_) = Mass_Electron

!! Next, are the Text Designations for the Species

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
!\
! This is for titan GITM radiation code
! This routine simply reads in the .txt file
! Containing the appropriate line strenghs,
! Line halfwidths, and so forth.
!/
      Use ModIoUnit, only : UnitTmp_

      implicit none

      integer :: iline
      real, dimension(1:rotlines) :: junk1
      real, dimension(1:rotlines) :: junk2
      real, dimension(1:rotlines) :: junk3

  open(UNIT = UnitTmp_, FILE = 'DataIn/NewLogIntensity_Minimal.txt', STATUS='OLD', &
       ACTION = 'READ')
  110 FORMAT(9(ES10.3,1X))

  do iline = 1,rotlines
   read(UnitTmp_,110) &
       pmat(iline,1), &
       pmat(iline,2), &
       pmat(iline,3), &
       pmat(iline,4), &
       pmat(iline,5), &
       pmat(iline,6), &
       pmat(iline,7), &
       pmat(iline,8), &
       pmat(iline,9)
  enddo
  close(Unit = UnitTmp_)

  open(UNIT = UnitTmp_, FILE = 'DataIn/New116HCNLines_Minimal.txt', STATUS='OLD', &
       ACTION = 'READ')

  112 FORMAT(F10.6, 1X, E9.3, 1X, E9.3, 1X, F6.4, 1X, F5.3, 1X, F9.4, 1X, F4.2)
  do iline = 1,rotlines
   read(UnitTmp_,112) &
       freqw(iline), &
       STref(iline), &
    Einstein(iline), &
     ALPHAL0(iline), &
     ALPHAS0(iline), &
        Eref(iline), &
    ALPHALExp(iline)
       freqhz(iline) = Speed_Light*100.0*freqw(iline)
  end do
  close(Unit = UnitTmp_)

      Qd(1,1) = 1.0e0
      Qd(1:2,2) = (/.6521451548e0, .3478548451e0/)
      Qd(1:4,3) = (/.3626837833e0, .3137066458e0, .2223810344e0,  &
                   .1012285362e0/)
      Qd(1:8,4) = (/.1894506104e0, .1826034150e0, .1691565193e0, .1495959888e0,  &
                    .1246289712e0, .0951585116e0, .0622535239e0, .0271524594e0/)

      Qf1(1,1) = .5773502691e0
      Qf1(1:2,2) = (/.3399810435e0, .8611363115e0/)
      Qf1(1:4,3) = (/.1834346424e0, .5255324099e0, .7966664774e0,  &
                   .9602898564e0/)
      Qf1(1:8,4) = (/.0950125098e0, .2816035507e0, .4580167776e0, .6178762444e0,  &
                    .7554044084e0, .8656312023e0, .9445750230e0, .9894009350e0/)


! 9-point gaussian quadrature
     FreqQuadAbscissas = &
      (/ 0.0000000000000000,  &
         0.8360311073266358, &
         0.9681602395076261, &
         0.3242534234038089, &
         0.6133714327005904/)

     FreqQuadWeights = &
      (/ 0.3302393550012598, &
         0.1806481606948574, &
         0.0812743883615744, &
         0.3123470770400029, &
         0.2606106964029354/)

! Altitude Gaussian Quadrature Points
     AltQuadAbscissas = &
      (/ 0.0000000000,  0.0936310658547334, &
                        0.1864392988279916, &
                        0.2776090971524970, &
                        0.3663392577480734, &
                        0.4518500172724507, &
                        0.5333899047863476, &
                        0.6102423458363790, &
                        0.6817319599697428, &
                        0.7472304964495622, &
                        0.8061623562741665, &
                        0.8580096526765041, &
                        0.9023167677434336, &
                        0.9386943726111684, &
                        0.9668229096899927, &
                        0.9864557262306425, &
                        0.9974246942464552 /)

     AltQuadWeights = &
      (/ 0.0937684461602100, 0.0933564260655961, &
                             0.0921239866433168, &
                             0.0900819586606386, &
                             0.0872482876188443, &
                             0.0836478760670387, &
                             0.0793123647948867, &
                             0.0742798548439541, &
                             0.0685945728186567, &
                             0.0623064825303175, &
                             0.0554708466316636, &
                             0.0481477428187117, &
                             0.0404015413316696, &
                             0.0323003586323290, &
                             0.0239155481017495, &
                             0.0153217015129347, &
                             0.0066062278475874 /)

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

  subroutine init_topography
    return
  end subroutine init_topography

end module ModPlanet
