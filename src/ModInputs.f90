! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

module ModInputs

  use ModConstants
  use ModPlanet
  use ModIoUnit, only : UnitTmp_
  use ModKind, only:    Real8_

  implicit none

  logical :: iRhoOutputList=.true.
  logical :: iNeutralDensityOutputList(nSpeciesTotal)=.true.
  logical :: iNeutralWindOutputList(3)=.true.
  logical :: iIonDensityOutputList(nIons)=.true.
  logical :: iIonWindOutputList(3)=.true.
  logical :: iTemperatureOutputList(3)=.true.

  integer                   :: useDART = 0 !alexey, default is to not use DART

  integer, parameter        :: iCharLen_     = 400

  character (len=iCharLen_) :: outputDir = "UA/data"
  character (len=iCharLen_) :: logDir = "UA/data"
  character (len=iCharLen_) :: restartOutDir = "UA/restartOUT"
  character (len=iCharLen_) :: restartInDir = "UA/restartIN"
  
  integer                   :: iOutputUnit_  = UnitTmp_
  integer                   :: iInputUnit_   = UnitTmp_
  integer                   :: iRestartUnit_ = UnitTmp_

  integer                   :: iLogFileUnit_ = 92
  integer                   :: iCodeInfoFileUnit_ = 93
  logical                   :: IsOpenLogFile = .false.
  real                      :: DtLogFile = 60.0

  integer, parameter        :: nInputMaxLines = 10000
  integer                   :: nInputLines
  character (len=iCharLen_) :: cInputText(nInputMaxLines)

  character (len=iCharLen_) :: cInputFile = "UAM.in"

  character (len=iCharLen_) :: cAMIEFileSouth = "none"
  character (len=iCharLen_) :: cAMIEFileNorth = "none"

  character (len=iCharLen_) :: PotentialModel
  character (len=iCharLen_) :: AuroralModel

  logical :: UseCCMCFileName = .false.

  logical :: UseSecondsInFilename = .true.    !xianjing

  logical :: UseIMF = .true.
  logical :: UseHpi = .true.

  logical :: UseIRHeating      = .false.
  
  !!! Xing Meng Nov 2018 to use ISR E field in a local region + Weimer elsewhere
  logical :: UseRegionalAMIE = .false.
  logical :: UseTwoAMIEPotentials = .false.
  real(Real8_) :: AMIETimeStart, AMIETimeEnd
  real    :: AMIELonStart = 208.0
  real    :: AMIELonEnd = 216.0
  real    :: AMIELatStart = 65.0
  real    :: AMIELatEnd = 70.0
  real    :: AMIEBoundaryWidth = 4.0  ! lat and lon width to transit to Weimer solution

  logical :: UseNewellAurora   = .false.
  logical :: UseNewellAveraged = .true.
  logical :: UseNewellMono     = .false.
  logical :: UseNewellWave     = .false.
  logical :: DoNewellRemoveSpikes = .true.
  logical :: DoNewellAverage      = .true.

  logical :: UseOvationSME     = .false.
  logical :: UseOvationSMEMono = .false.
  logical :: UseOvationSMEWave = .false.
  logical :: UseOvationSMEIon  = .false.
  
  logical :: UseAeModel        = .false.

  logical :: UseFangEnergyDeposition = .false.
  
  logical :: IsKappaAurora = .false.
  real :: AuroraKappa = 3
  real :: AveEFactor = 1.0
  ! This is true because FR&E is default.
  logical :: NormalizeAuroraToHP = .true.

  logical :: xxx

  
  logical :: UseCusp = .false.
  real :: CuspAveE = 0.1
  real :: CuspEFlux = 2.0
  real :: CuspMltHalfWidth = 1.5
  real :: CuspLatHalfWidth = 1.0

  logical :: DoOverwriteIonosphere = .false.
  logical :: DoOverwriteWithIRI    = .true.
  logical :: DoOverwriteWithSami   = .false.
  character (len=iCharLen_) :: SamiInFile

  character (len=iCharLen_) :: TypeLimiter = "minmod"

  integer, dimension(7) :: iStartTime

  real :: CPUTimeMax = 1.0e32

  logical :: UseVariableInputs = .false.
  logical :: IsFramework = .false.

  logical :: DoRestart = .false.

  integer :: iAltTest = -1
  integer :: iDebugLevel = 0
  logical :: UseBarriers = .false.
  logical :: DoCheckForNans = .false.
  integer :: nSteps = 10

  logical :: UsePerturbation = .false., DuringPerturb = .false.

  real :: CFL = 0.25
  real :: FixedDt = 1.0e32
  
  integer :: nOutputTypes = 0
  integer, parameter :: nMaxOutputTypes = 50
  character (len=iCharLen_), dimension(nMaxOutputTypes) :: OutputType

  real :: DtPlot(nMaxOutputTypes)
  real :: DtPlotSave(nMaxOutputTypes)
  real :: PlotTimeChangeDt(nMaxOutputTypes)
  real(Real8_) :: PlotTimeChangeStart, PlotTimeChangeEnd

  logical :: DoAppendFiles = .false.

  ! Xing Meng March 2020 -- Save HIME type output in a user-specified region
  real    :: HIMEPlotLonStart = 208.0
  real    :: HIMEPlotLonEnd = 216.0
  real    :: HIMEPlotLatStart = 65.0
  real    :: HIMEPlotLatEnd = 70.0

  real :: DtRestart   = 60.0*60.0
  real :: DtReport    =  1.0*60.0
  real :: DtAurora    = 60.0*1.0
  real :: DtPotential = 60.0*1.0
  real :: DtGlow      = 60.0
  real :: TimeDelayHighLat = 0.0
  real :: TimeDelayEUV     = 0.0

  real :: f107  = 150.0
  real :: f107a = 150.0
  integer :: iModelSolar = 0

  real :: AltMin =  95.0 * 1000.0
  real :: AltMax = 500.0 * 1000.0

  real :: BetaLimiter = 2.0

  real :: ConcentrationLatitude = 45.0
  real :: StretchingPercentage = 0.0
  real :: StretchingFactor = 1.0
  logical :: NewStretchedGrid = .false.
  real :: StretchWidth = 10.0

  logical :: UseTopography = .false.
  real :: AltMinUniform = 0.0
  real :: AltMinIono=80.0 ! in km

  real :: TempMax = 1000.0
  real :: TempMin =  200.0
  real :: TempWidth    =  25.0*1e3
  real :: TempHeight   = 150.0*1e3

  real :: LogRho0

  integer :: nBlocksLon = 2
  integer :: nBlocksLat = 2

  logical :: Is1D = .false.
  logical :: IsFullSphere = .false.
  logical :: UseGlow = .false.

  real    :: LonStart = 0.0
  real    :: LonEnd   = 0.0
  real    :: LatStart = -pi/4.0
  real    :: LatEnd   =  pi/4.0

  logical :: UseStretchedAltitude = .true.

  logical :: UseApex = .true.
  logical :: UseMSIS = .true.
  real, dimension(25) :: sw_msis = 1.0
  logical :: UseIRI  = .true.
  logical :: UseMSISTides  = .true.
  logical :: UseMSISOnly   = .false.
  logical :: UseGSWMTides  = .false.
  logical :: UseWACCMTides = .false.
  logical :: UseMSISDiurnal = .true.
  logical :: UseMSISSemidiurnal = .true.
  logical :: UseMSISTerdiurnal = .true.
  logical :: UseStatisticalModelsOnly = .false.
  real    :: DtStatisticalModels = 3600.0
  logical :: UseOBCExperiment = .false.
  real :: MsisOblateFactor = 0.0

  logical :: UseGswmComp(4) = .true.

  real :: MagneticPoleRotation = 0.0
  real :: MagneticPoleTilt = 0.0
  real :: xDipoleCenter = 0.0
  real :: yDipoleCenter = 0.0
  real :: zDipoleCenter = 0.0

  logical :: IsFixedTilt = .false.

  !\
  ! These are things for the ion precipitation
  !/

  logical :: UseIonPrecipitation = .false.

  logical :: UseDamping = .false.

!! EddyVelocity Terms
  real    :: MaximumVerticalVelocity = 1000.0

  !\
  ! These are logicals to turn off and on forcing terms:
  !/

  logical :: UsePressureGradient = .true.
  logical :: UseGravity          = .true.
  logical :: UseIonDrag          = .true.
  logical :: UseViscosity        = .true.
  real    :: TestViscosityFactor = 1.0
  logical :: UseCoriolis         = .true.
  logical :: UseGravityWave      = .false.

  logical :: UseHorAdvection     = .true.
  logical :: UseVerAdvection     = .true.
  logical :: UseNeutralFriction  = .true.

  logical :: UseAUSMSolver     = .false.

  logical :: UseIonPressureGradient = .true.
  logical :: UseIonGravity          = .true.
  logical :: UseNeutralDrag         = .true.
  logical :: UseExB                 = .true.
  logical :: UseImplicitChemistry = .false.
  logical :: IsAsymmetric          = .false.
  Real :: BetaPointImpl         = 1.0

  logical :: UseDynamo              = .false.
  real    :: DynamoHighLatBoundary  = 65.0
  integer :: nItersMax              = 500
  real    :: MaxResidual            = 1.0
  logical :: IncludeCowling         = .false.
  real    :: DynamoLonAverage       = 10.0

  logical :: UseImprovedIonAdvection = .false.
  logical :: UseNighttimeIonBCs = .false.
  real :: MinTEC = 2.0
  logical :: UseImplicitFieldAlignedMomentum = .false.
  real    :: DivIonVelCoef = 1.0

  logical :: UseGitmBCs = .false.
  character(len=iCharLen_) :: GitmBCsDir

  logical :: UseSolarHeating   = .true.
  logical :: UseJouleHeating   = .true.
  logical :: UseAuroralHeating = .true.
  logical :: UseNOCooling      = .true.
  logical :: UseOCooling       = .true.
  logical :: UseConduction     = .true.
  logical :: UseTurbulentCond = .true.

  logical :: UseDiffusion      = .false.
  logical :: UseVerAdvectionT  = .true.

  logical :: UseCO2Cooling = .true.
  real    :: CO2ppm = 225.0

  ! Allow the user to change the planet's characteristics:
  real :: RotationPeriodInput = Rotation_Period
  real :: OmegaBodyInput = 2.0*pi/Rotation_Period
  real :: HoursPerDayInput = Rotation_Period/3600.0
  real :: DaysPerYearInput = DaysPerYear
  real :: PlanetTiltInput = Tilt
  

  real :: PhotoElectronHeatingEfficiency = 0.0
  real :: NeutralHeatingEfficiency = 0.05

  real :: KappaTemp0 = 5.6e-4
  real :: ThermalConduction_AO2 = 3.6e-4
  real :: ThermalConduction_AO  = 5.6e-4
  real :: ThermalConduction_s   = 0.726
  !! Pawlowski says AO2 = 3.6e-4 - 5.6e-4
  !!                AO  = 5.6e-4 - 7.6e-4
  !!                s   = 0.69 - 0.75

  real :: EddyDiffusionCoef = 0.0
  real :: EddyDiffusionPressure0 = 0.0
  real :: EddyDiffusionPressure1 = 0.0
  real :: Kappa1DCorrectionFactor = 45.0
  real :: Kappa1DCorrectionPower  = 1.75
  logical :: UseKappa1DCorrection = .false.

  logical :: UseIonChemistry     = .true.
  logical :: UseNeutralChemistry = .true.
  logical :: UseIonAdvection     = .true.

  logical :: DoCheckStopFile = .true.

  ! AGB: Setting physical limits for ionospheric dynamics
  real :: MaxVParallel = 100.0         
  real :: MaxEField = 0.1

  !\
  ! Methods for completing chemistry
  !/

  integer, parameter :: cSteadyStateChemType_ = 1
  integer, parameter :: cImplicitChemType_    = 2
  integer, parameter :: cSubCycleChemType_    = 3
  integer, parameter :: nChemTypes_       = 3

  character (len=100), dimension(nChemTypes_) :: sChemType 

  character (len=100)                         :: sInputIonChemType
  character (len=100)                         :: sInputNeutralChemType
  integer :: iInputIonChemType, iInputNeutralChemType

  real :: LogNS0(nSpecies)

  logical                   :: UseEUVAC     = .true.
  logical                   :: UseTobiska   = .true.
  logical                   :: UseAboveHigh = .true.
  logical                   :: UseBelowLow  = .true.

  logical                   :: UseRidleyEUV = .false.
  
  logical                   :: UseEUVData =.false.
  character (len=iCharLen_) :: cEUVFile

  !\
  ! Eclipse Information
  !/
  logical :: IncludeEclipse = .false.
  real(Real8_) :: EclipseStartTime
  real(Real8_) :: EclipseEndTime
  real :: EclipseStartY, EclipseEndY
  real :: EclipseStartZ, EclipseEndZ
  real :: EclipsePeak, EclipseMaxDistance
  real :: EclipseExpAmp, EclipseExpWidth
  
  ! These are Mars Specific, but ignored by other codes:
  ! Some are modified in Planet.f90 (set_planet_defaults)
  real :: DtLTERadiation = 100.0*Rotation_Period

  !!!!!!!!!!!!!!! DUST
  Logical :: UseDustDistribution = .False.
  ! Setting the depth to which dust is mixed based on a reference dust opacity
  !This can now be read in as a horizontal distribution, or as a constant value
  !using UAM.in
  character (len=iCharLen_) :: cDustFile="NotSet"
  character (len=iCharLen_) :: cConrathFile="NotSet"
  real :: CONRNU_temp = 0.03    ! Standard value  ~25km half-height
  !     real :: CONRNU = 0.003   ! ~50 km half-height
  !     real  :: CONRNU = 0.5     ! ~10 km half-height
  ! Global mean dust opacity
  real :: TAUTOT_temp  = 0.3 !do not set to 0.0 or less

  real :: dtDust
  character (len=iCharLen_) :: DustFileType

  !!!!!!!!!!!!!!!!!!!!!!

  !   RPTAU   :  Reference Pressure optical depth;  6.1 mbar for now
  real, parameter :: RPTAU  = 6.1

  ! Top of the shortwave calculation for lower atmosphere radiation code(mbars)
  real, parameter :: PRAD  = 1.0E-6 !mb (~115 km)

  ! Top of the longwave calculation for lower atmosphere radiation code(PASCALS)
  ! and where radcool begins
  real, parameter :: PLONG = 4.0E-1 ! Pa  (also 4.0  ubar)

  ! Ls Variables
  real, dimension(7) :: Ls_a
  real, dimension(7) :: Ls_tau
  real, dimension(7) :: Ls_phi
	  
  DATA Ls_a / 0.007, 0.006, 0.004, 0.004, 0.002, 0.002, 0.002 /
  DATA Ls_tau / 2.2353, 2.7543, 1.1177, 15.7866, 2.1354, 2.4694, 32.8493 /
  DATA Ls_phi / 49.409, 168.173, 191.837, 21.736, 15.704, 95.528, 49.095 /

  !\
  ! Variables for Wave Perturbation (WP) model (in user.f90)
  ! to characterize surface motions during tsunamis or earthquakes
  !/ 
  logical :: UseBcPerturbation = .false.
  integer :: iTypeBcPerturb
  real    :: RefLon, RefLat, EpicenterLon, EpicenterLat
  real    :: PerturbTimeDelay, PerturbDuration, SeisWaveTimeDelay   ! in sec 
  real    :: PerturbWaveDirection       ! in deg counter-clockwise to east
  real    :: PerturbWaveSpeed           ! in m/sec
  real    :: PerturbWaveHeight          ! in m
  real    :: PerturbWavePeriod          ! in sec
  real    :: EpiDistance                ! in m
  character (len=iCharLen_) :: cSurfacePerturbFileName 
  integer, parameter :: nMaxPerturbFreq = 100
  real    :: FFTReal(nMaxPerturbFreq), FFTImag(nMaxPerturbFreq), &
       PerturbWaveFreq(nMaxPerturbFreq)
  integer :: nPerturbFreq 

contains

  ! -------------------------------------------------------------
  ! set_strings
  !     Here we initialize all of the strings which need to be
  !     initialized in the code.
  ! -------------------------------------------------------------

  subroutine set_strings

    sChemType(cSteadyStateChemType_) = "steady"
    sChemType(cImplicitChemType_)    = "implicit"
    sChemType(cSubCycleChemType_)    = "subcycle"

  end subroutine set_strings

  ! -------------------------------------------------------------
  ! set_defaults
  ! -------------------------------------------------------------

  subroutine set_defaults

    use ModTime
    use ModPlanet, only:IsEarth

    call set_strings

    sInputIonChemType     = sChemType(cSubCycleChemType_)
    sInputNeutralChemType = sChemType(cSubCycleChemType_)
    iInputIonChemType     = cSubCycleChemType_
    iInputNeutralChemType = cSubCycleChemType_

    PotentialModel = "Weimer05"
    AuroralModel = "ihp"

    dTAurora = 120.0
    
    if (IsEarth) then 
       PhotoElectronHeatingEfficiency = 0.06
    endif

    tSimulation = 0.0

    iStartTime(1) = 1999
    iStartTime(2) =    3
    iStartTime(3) =   21
    iStartTime(4) =    0
    iStartTime(5) =    0
    iStartTime(6) =    0
    iStartTime(7) =    0

    iTimeArray = iStartTime

    call time_int_to_real(iStartTime, CurrentTime)

    nOutputTypes = 0
    DtPlot       = -1.0
    DtRestart    = -1.0

    ! Initiate Neutral Variables From MSIS90

    f107a = 150.0
    f107 = 150.0

    iModelSolar = 0

    UseApex   = .false.
    DoRestart = .false.

    call set_planet_defaults

  end subroutine set_defaults

end module ModInputs

