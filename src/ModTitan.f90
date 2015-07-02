!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPlanet

  use ModConstants
  use ModSizeGITM

  implicit none

!! Major Species Indices
!! Have Individual Species Velocities (Radial)

  integer, parameter :: nSpecies = 5
  integer, parameter :: iN2_    = 1
  integer, parameter :: iCH4_   = 2
  integer, parameter :: iAr_    = 3
  integer, parameter :: iHCN_   = 4
  integer, parameter :: iH2_    = 5

! Minor Species

  integer, parameter :: i15N2_  = 6
  integer, parameter :: i13CH4_ = 7
  integer, parameter :: i3CH2_  =  8
  integer, parameter :: i1CH2_  =  9
  integer, parameter :: iCH3_   =  10
  integer, parameter :: iCH_    =  11
  integer, parameter :: iN4S_   =  12
  integer, parameter :: iH_     =  13
  integer, parameter :: iC2H4_  =  14
  integer, parameter :: iH2CN_  =  15
  integer, parameter :: nSpeciesTotal = iH2CN_

! Major Ions (4):  Most Important to MWACM code

  integer, parameter  :: iH2CNP_ = 1
  integer, parameter  :: iC2H5P_ = 2
  integer, parameter  :: iCH3P_  = 3
  integer, parameter  :: iN2P_   = 4
  integer, parameter  :: iNP_    = 5
  integer, parameter  :: ie_     = 6
  integer, parameter  :: nIons   = ie_
  integer, parameter  :: nIonsAdvect = 2

  character (len=20) :: cSpecies(nSpeciesTotal)
  character (len=20) :: cIons(nIons)

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

!! Seems to be used only for GLOW Code
  integer, parameter :: nEmissionWavelengths = 1
  integer, parameter :: nPhotoBins = 1

  ! When you want to program in emissions, you can use this...
  integer, parameter :: nEmissions = 10

  real, parameter :: GC_Titan                = 1.354                   ! m/s^2
  real, parameter :: RP_Titan                = 1377684.3               ! seconds
  real, parameter :: R_Titan                 = 2575.00*1000.0          ! meters
  real, parameter :: DP_Titan                = 0.0

  real, parameter :: Gravitational_Constant = GC_Titan
  real, parameter :: Rotation_Period        = RP_Titan
  real, parameter :: RBody                  = R_Titan
  real, parameter :: DipoleStrength         = DP_Titan

  real, parameter :: OMEGABody              = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period / 3600.0
  real, parameter :: Tilt = 26.70                                       ! Assume the Obliquity of Saturn Here!

  ! This is the Vernal Equinox at Midnight (Ls = 0!!!)
  ! Earth-Mars clocks are set from this epoch at vernal equinox
  integer, parameter :: iVernalYear   = 1980
  integer, parameter :: iVernalMonth  =    2
  integer, parameter :: iVernalDay    =   14
  integer, parameter :: iVernalHour   =   12
  integer, parameter :: iVernalMinute =    0
  integer, parameter :: iVernalSecond =    0

  real, parameter :: SunOrbit_A = 9.54
  real, parameter :: SunOrbit_B = 0.04
  real, parameter :: SunOrbit_C = 0.15
  real, parameter :: SunOrbit_D = 0.00
  real, parameter :: SunOrbit_E = 0.00

  real, parameter :: DaysPerYear = 10759.53
  real, parameter :: SecondsPerYear = DaysPerYear * Rotation_Period
  
    !Used as a damping term in Vertical solver.
  real, dimension(nAlts) :: VertTau = 1.0e9 

!! Ensure that Planetary-specific Physics is Enabled
  logical :: IsEarth = .false.
  logical :: IsMars  = .false.
  logical :: IsJupiter  = .false.
  logical :: IsTitan = .true.
  real, parameter :: PlanetNum = 0.06

  logical :: NonMagnetic  = .true.

  character (len=10) :: cPlanet = "Titan"

!! Not Known What this is
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


!! Radiative Transfer Parameters Listed Here for HCN
!! JMB 1-22-2009

   integer, parameter :: rotlines = 116     !! # of HCN Rotational Lines
   integer, parameter :: rotfreqs = 16      !! # of Frequency Gaussian Points
   integer, parameter :: rotpts   = 16      !! # of Altitude Gassuian Points
   real, dimension(rotlines) :: freqw       !! Rotational Lines in Wave Numbers (cm^-1)
   real, dimension(rotlines) :: freqhz      !! Rotational Lines (in Hz)
   real, dimension(rotlines) :: gammal      !! Lorentz Halfwidth (in cm^-1)
   real, dimension(rotlines,9) :: pmat      !! Intensity Interpolation Matrix
   real, parameter :: atmtopas = 1.0/(1.01295e+05)

!! Gaussian Quadrature Points
   real, dimension(8,4) :: Qd, Qf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are Modified for Titan by JMB: 1/22/09
  ! -- Most source state Dij = Dji (check)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! These are Aij coefficients from B&K (1973) formulation: Aij*1.0E+17
  ! Use Jared Bell's Titan GITM formulation for Mars GITM

!!! Dimension (nSpecies, nSpecies)
!  real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e16 * reshape( (/ &
!     !-----------------------------------------------------------+
!     ! i=N2      CH4     Ar       HCN     H2      15N2    13CH4  
!     !-----------------------------------------------------------+
!       0.0000, 7.3400, 6.6344,  5.1834, 18.8000,  5.0460, 5.4766, &            ! N2
!       7.3400, 0.0000, 5.7537,  5.1132, 23.0000,  6.6750, 5.6450, &            ! CH4
!       6.6344, 5.7537, 0.0000,  6.5184, 17.2550,  4.6050, 4.6230, &            ! Ar
!       5.1834, 5.1132, 6.5184,  0.0000, 17.4510,  5.1390, 4.9440, &            ! HCN
!      18.8000,23.0000,17.2550, 17.4510,  0.0000, 18.8800,15.9670, &            ! H2
!       5.0460, 6.6750, 4.6050,  5.1390, 18.8800,  0.0000, 6.4758, &            ! 15N2
!       5.4766, 5.6450, 4.6230,  4.9440, 15.9670,  6.4768, 0.0000 /), &         ! 13CH4
!       (/nSpecies,nSpecies/) ) 
!
!  real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
!     !-----------------------------------------------------------+
!     ! i=N2      CH4     Ar       HCN     H2      15N2    13CH4  
!     !-----------------------------------------------------------+
!       0.0000, 0.7500, 0.7520,  0.8100,  0.8200,  0.7500, 0.7500, &            ! N2
!       0.7500, 0.0000, 0.7850,  0.7650,  0.7650,  0.7500, 0.7500, &            ! CH4
!       0.7520, 0.7850, 0.0000,  0.7520,  0.7500,  0.7500, 0.7500, &            ! Ar
!       0.8100, 0.7650, 0.7520,  0.0000,  0.7500,  0.7500, 0.7500, &            ! HCN
!       0.8200, 0.7650, 0.7500,  0.7500,  0.0000,  0.7500, 0.7500, &            ! H2
!       0.7500, 0.7500, 0.7500,  0.7500,  0.7500,  0.0000, 0.7500, &            ! 15N2
!       0.7500, 0.7500, 0.7500,  0.7500,  0.7500,  0.7500, 0.0000 /), &         ! 13CH4
!       (/nSpecies,nSpecies/) ) 

!!! REDUCED TITAN ATMOSPHERE (nSpecies = 4)
!
  !! Stuff for initial conditions

  real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e16 * reshape( (/ &
     !-----------------------------------------------------------+
     ! i=N2      CH4     Ar       HCN     H2    
     !-----------------------------------------------------------+
       0.0000, 7.3400, 6.6344,  5.1834, 18.8000, &            ! N2
       7.3400, 0.0000, 5.7537,  5.1132, 23.0000, &            ! CH4
       6.6344, 5.7537, 0.0000,  6.5184, 17.2550, &            ! Ar
       5.1834, 5.1132, 6.5184,  0.0000, 17.4510, &            ! HCN
      18.8000,23.0000,17.2550, 17.4510,  0.0000 /), &         ! H2
       (/nSpecies,nSpecies/) ) 

  real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
     !-----------------------------------------------------------+
     ! i=N2      CH4     Ar       HCN     H2      
     !-----------------------------------------------------------+
       0.0000, 0.7500, 0.7520,  0.8100,  0.8200, &            ! N2
       0.7500, 0.0000, 0.7850,  0.7650,  0.7650, &            ! CH4
       0.7520, 0.7850, 0.0000,  0.7520,  0.7500, &            ! Ar
       0.8100, 0.7650, 0.7520,  0.0000,  0.7500, &            ! HCN
       0.8200, 0.7650, 0.7500,  0.7500,  0.0000 /), &         ! H2
       (/nSpecies,nSpecies/) ) 


  real , Dimension(-1:nAlts + 2) :: newalt
  real , Dimension(-1:nAlts + 2) :: junk
  real , Dimension(-1:nAlts + 2) :: InTemp
  real , Dimension(-1:nAlts + 2) :: IneTemp
  real , Dimension(-1:nAlts + 2) :: InITemp
  real , Dimension(-1:nAlts + 2,nSpeciesTotal) :: InNDensityS 
  real , Dimension(-1:nAlts + 2,nIons) :: InIDensityS

!! Used in calc_electrodynamics
  real, parameter:: AltMinIono=500.0 ! in km

!! HCN Rotational Cooling Frquency
  real :: dTCooling = 86400.00

!! Magnetospheric Heating Sources (From Cravens et al [2005])
  real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: EHeatingRate
  real , Dimension(1:nLons,1:nLats,1:nAlts,nSpeciesTotal,nBlocksMax) :: MagDissRateS
  real , Dimension(1:nLons,1:nLats,1:nAlts,nIons,nBlocksMax) :: MagIonRateS

!! Aerosol Calculation Variables
  real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: AerosolClassA
  real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: AerosolClassB
  real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: AerosolClassC

!! Aerosol Loss Rates
  real , Dimension(1:nLons,1:nLats,1:nAlts,nSpeciesTotal,nBlocksMax) :: AerosolTrappingLoss
  real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: H2AerosolProduction

!! Isotope Chemistry Variables
  real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: IsotopeScaling


contains

  subroutine init_planet

    use ModTime
    use ModIoUnit, only : UnitTmp_

    integer :: iTime(7), iiAlt

!! Majors
    Mass(iN2_)    =  14.00674 * AMU * 2
    Mass(iCH4_)   =   1.0079 * AMU * 4 + 12.011 * AMU
    Mass(iH2_)    =   1.0079 * AMU * 2
    Mass(iAr_)    =  40.0000 * AMU
    Mass(iHCN_)   =   1.0079 * AMU + 14.00674 * AMU + 12.011 * AMU
    Mass(i15N2_)  =  14.00674 * AMU + 15.00674 * AMU
    Mass(i13CH4_) =  13.01100 * AMU + 1.0079 * AMU * 4

!! Minors
    Mass(i3CH2_)    =  1.0079 * AMU * 2 + 12.011 * AMU
    Mass(i1CH2_)    =  1.0079 * AMU * 2 + 12.011 * AMU
    Mass(iCH3_)     =  1.0079 * AMU * 3 + 12.011 * AMU
    Mass(iCH_)      =  1.0079 * AMU + 12.011 * AMU
    Mass(iN4S_)     =  14.00674 * AMU 
    Mass(iH_)       =  1.0079 * AMU 
    Mass(iC2H4_)    =  1.0079 * AMU * 4 + 12.011 * AMU * 2
    Mass(iH2CN_)    =  1.0079 * AMU * 2 + 14.00674 * AMU + 12.011 * AMU

!! Ions
    MassI(iN2P_)    =  14.00674 * AMU * 2
    MassI(iNP_)     =  14.00674 * AMU
    MassI(iH2CNP_)  =  1.0079 * AMU * 2 + 14.00674 * AMU + 12.011 * AMU
    MassI(iCH3P_)   =  1.0079 * AMU * 3 + 12.011 * AMU
    MassI(iC2H5P_)  =  1.0079 * AMU * 5 + 12.011 * AMU * 2
    MassI(ie_) = Mass_Electron

!! Next, are the Text Designations for the Species
!! Majors
    cSpecies(iN2_)   = "N!D2!N"
    cSpecies(iCH4_)  = "CH!D4!N"
    cSpecies(iH2_)   = "H!D2!N"
    cSpecies(iAr_)   = "Ar"
    cSpecies(iHCN_)  = "HCN"
    cSpecies(i15N2_)   = "!U15!NN!D2!N"
    cSpecies(i13CH4_)  = "!U13!NCH!D4!N"

!! Minors
    cSpecies(i3CH2_)  = "!U3!NCH!D2!N"
    cSpecies(i1CH2_)  = "!U1!NCH!D2!N"
    cSpecies(iCH3_)   = "CH!D3!N"
    cSpecies(iCH_)    = "CH"
    cSpecies(iN4S_ )  = "N(!U4!NS)"
    cSpecies(iCH_)    = "H"
    cSpecies(iC2H4_)  = "C!D2!NH!D4!N"
    cSpecies(iH2CN_)  = "H!D2!NCN"

!! Ions

    cIons(iH2CNP_)  = "H!D2!NCN!U+!N"
    cIons(iC2H5P_)  = "C!D2!NH!D5!N!U+!N"
    cIons(iCH3P_)   = "CH!D3!U+!N"
    cIons(iN2P_)   = "N!D2!U+!N"
    cIons(iNP_)   = "N!U+!N"
    cIons(ie_)     = "e!U-!N"

!! Used for CP Calculations.
!! Clearly, this does not allow for Temperature Variations.

    Vibration(iN2_)   = 5.0 + 2.0
    Vibration(iH2_)   = 5.0 + 2.0
    Vibration(iCH4_)  = 6.5374 + 2.0
    Vibration(iHCN_)  = 5.0 + 2.0
    Vibration(iAr_)   = 3.0
    Vibration(i13CH4_)  = 6.5374 + 2.0
    Vibration(i15N2_)   = 5.0 + 2.0

    itime = 0
    itime(1) = iVernalYear
    itime(2) = iVernalMonth
    itime(3) = iVernalDay
    itime(4) = iVernalHour
    itime(5) = iVernalMinute
    itime(6) = iVernalSecond
    call time_int_to_real(itime, VernalTime)

    write(*,*) 'Reading in the Titan_input.txt'

    open(UNIT = UnitTmp_, FILE = 'UA/DataIn/Titan_TA_500km_10km.txt', STATUS='OLD', ACTION = 'READ')
111 FORMAT(F7.2,1X,F8.4,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
              1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
              1X,ES10.3)

    InNDensityS(:,:) = 1.0e+3
    InIDensityS(:,:) = 1.0e+3

    do iiAlt = -1,nAlts+2
       read(UnitTmp_,111) &
            newalt(iiAlt), &
            InTemp(iiAlt), &
            !
            InNDensityS(iiAlt,iH_), &
            InNDensityS(iiAlt,iH2_), &
            InNDensityS(iiAlt,iCH_), &
            InNDensityS(iiAlt,i1CH2_), &
            InNDensityS(iiAlt,i3CH2_), &
            InNDensityS(iiAlt,iCH3_), &
            InNDensityS(iiAlt,iCH4_), &
            InNDensityS(iiAlt,iC2H4_), &
            InNDensityS(iiAlt,iN4S_), &
            InNDensityS(iiAlt,iN2_), &
            InNDensityS(iiAlt,iHCN_), &
            InNDensityS(iiAlt,iH2CN_)
    end do
    close(Unit = UnitTmp_)

    open(UNIT = UnitTmp_, FILE = 'UA/DataIn/HCN_500km_10km.txt', STATUS='OLD', ACTION = 'READ')
131 FORMAT(F7.2,1X,ES10.3)
    do iiAlt = -1,nAlts+2
       read(UnitTmp_,131) &
            newalt(iiAlt), &
            InNDensityS(iiAlt,iHCN_)
     enddo
    close(Unit = UnitTmp_)

    open(UNIT = UnitTmp_, FILE = 'UA/DataIn/TeKe_500km_10km.txt', STATUS='OLD', ACTION = 'READ')
125 FORMAT(1X,F9.4,1X,ES10.3)
    do iiAlt = -1,nAlts+2
       read(UnitTmp_,125) &
            IneTemp(iiAlt), &
               junk(iiAlt)
     enddo
    close(Unit = UnitTmp_)


  end subroutine init_planet

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

      open(UNIT = UnitTmp_, FILE = 'DataIn/intensity116.txt', STATUS='OLD', &
           ACTION = 'READ')

110 FORMAT(9ES21.10)

!\
! P is the polynomial interpolation coefficients for the Intensities
! as a function of temperature (K)
!/

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
!\
! Intensity Matrix is Properly Read into Code! 3-03-2006
!/


      open(UNIT = UnitTmp_, FILE = 'DataIn/HCN_116.txt', STATUS='OLD', &
           ACTION = 'READ')

111 FORMAT(F10.6, F7.4)

  do iline = 1,rotlines
   read(UnitTmp_,111) &
       freqw(iline), &
       gammal(iline)
       freqhz(iline) = Speed_Light*100.0*freqw(iline)
  end do
  close(Unit = UnitTmp_)

! Sets the Matrices D, F for Frequency Gaussian Points
!
      Qd(1,1) = 1.0e0
      Qd(1:2,2) = (/.6521451548e0, .3478548451e0/) 
      Qd(1:4,3) = (/.3626837833e0, .3137066458e0, .2223810344e0,  &
                   .1012285362e0/)
      Qd(1:8,4) = (/.1894506104e0, .1826034150e0, .1691565193e0, .1495959888e0,  &
                    .1246289712e0, .0951585116e0, .0622535239e0, .0271524594e0/)

      Qf(1,1) = .5773502691e0 
      Qf(1:2,2) = (/.3399810435e0, .8611363115e0/)
      Qf(1:4,3) = (/.1834346424e0, .5255324099e0, .7966664774e0,  &
                   .9602898564e0/)
      Qf(1:8,4) = (/.0950125098e0, .2816035507e0, .4580167776e0, .6178762444e0,  &
                    .7554044084e0, .8656312023e0, .9445750230e0, .9894009350e0/)

      end subroutine init_radcooling 

  subroutine init_topography
    return
  end subroutine init_topography

end module ModPlanet
