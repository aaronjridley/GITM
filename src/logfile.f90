! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!----------------------------------------------------------------------------
! $ Id: $
!
! Author: Aaron Ridley, UMichigan
!
! Comments: Routines to produce the log file output.  The log file contains
!           useful physical outputs in ascii form, allowing the user to see
!           how the model is running without reading the GITM binaries.  It is
!           also used in the compilation tests to ensure that GITM is set
!           up successfuly given the desired configuration options.
!
! AGB 1/9/13: Add subsolar point and VTEC to log
!-----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! get_log_info: computes the physical values that will be output in the log.
!               This routine goes to each processor and performs computations
!               as appropriate.
!
! AGB 1/9/13: Added computation of VTEC at the subsolar point
!----------------------------------------------------------------------------

subroutine get_log_info(SSLon, SSLat, GlobalMinTemp, GlobalMaxTemp, &
     GlobalMinVertVel, GlobalMaxVertVel, AverageTemp, AverageVertVel, &
     TotalVolume, SSVTEC)

  use ModGITM

  real, intent(in)  :: SSLon, SSLat 
  real, intent(out) :: GlobalMinTemp, GlobalMaxTemp
  real, intent(out) :: GlobalMinVertVel, GlobalMaxVertVel
  real, intent(out) :: AverageTemp, AverageVertVel
  real, intent(out) :: TotalVolume, SSVTEC

  integer :: iBlock, iSpecies, iLon, iLat
  real :: rLon, rLat
  !--------------------------------------------------------------------------

  GlobalMaxTemp    = 0.0
  GlobalMinTemp    = 1.0e32
  GlobalMaxVertVel = 0.0
  GlobalMinVertVel = 1.0e32
  AverageTemp      = 0.0
  AverageVertVel   = 0.0
  TotalVolume      = 0.0
  SSVTEC           = -1.0e32

  do iBlock = 1, nBlocks

     call calc_rates(iBlock)

     GlobalMaxTemp = max(GlobalMaxTemp, &
          maxval(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     GlobalMinTemp = min(GlobalMinTemp, &
          minval(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     AverageTemp = AverageTemp + &
          sum(Temperature(1:nLons,1:nLats,1:nAlts,iBlock) &
          *   TempUnit(1:nLons,1:nLats,1:nAlts) &
          *   CellVolume(1:nLons,1:nLats,1:nAlts,iBlock))

     GlobalMaxVertVel = max(GlobalMaxVertVel, &
          maxval(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock)))
     GlobalMinVertVel = min(GlobalMinVertVel, &
          minval(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock)))
     do iSpecies = 1, nSpecies
        AverageVertVel = AverageVertVel + &
             sum(VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) &
             *   CellVolume(1:nLons,1:nLats,1:nAlts,iBlock))
     enddo

     TotalVolume = TotalVolume + &
          sum(CellVolume(1:nLons,1:nLats,1:nAlts,iBlock))
  enddo

  call calc_single_vtec_interp(SSLon, SSLat, SSVTEC)

  AverageTemp    = AverageTemp
  AverageVertVel = AverageVertVel / nSpecies

end subroutine get_log_info

!==============================================================================
! logfile: A routine to write a log file
!
! AGB 1/9/13: Added Subsolar location and VTEC at the subsolar point to the
!             output.  The subsolar location is computed in this subroutine
!             and the location passed into get_log_info to compute the VTEC.
!             MPI_REDUCE is then used to retrieve the appropriate VTEC.
!==============================================================================

subroutine logfile(dir)

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi
  use ModIndices
  use ModIndicesInterfaces
  use ModIoUnit, ONLY: io_unit_new
  use ModUtilities, ONLY: flush_unit

  implicit none

  character (len=*), intent(in) :: dir
  character (len=8) :: cIter

  real    :: minTemp, maxTemp, localVar, minVertVel, maxVertVel
  real    :: AverageTemp, AverageVertVel, TotalVolume, Bx, By, Bz, Vx, Hpi
  real    :: HPn, HPs, SSLon, SSLat, SSVTEC
  integer :: iError

  if (.not. IsOpenLogFile .and. iProc == 0) then

     call write_code_information(dir)

     IsOpenLogFile = .true.
     iLogFileUnit_ = io_unit_new()

     write(cIter,"(i8.8)") iStep

     open(unit=iLogFileUnit_, &
          file=trim(dir)//"/log"//cIter//".dat",status="replace")

     write(iLogFileUnit_,'(a)') "GITM2 log file"
     write(iLogFileUnit_,'(a,L2)') "## Inputs from UAM.in" 
     write(iLogFileUnit_,'(a,L2)') "# Resart=", dorestart
     write(iLogFileUnit_,'(4(a,f9.3))') "# Eddy coef: ", EddyDiffusionCoef, &
          " Eddy P0: ",EddyDiffusionPressure0,&
          " Eddy P1: ",EddyDiffusionPressure1
     write(iLogFileUnit_,'(2(a,L2))') "# Statistical Models Only: ", &
          usestatisticalmodelsonly, " Apex: ",useApex
     if (useEUVdata) then
        write(iLogFileUnit_,'(a,L2,a)') "# EUV Data: ", useEUVdata, "File: ", &
             cEUVFile
     else
        write(iLogFileUnit_,'(a,L2)') "# EUV Data: ", useEUVdata
     endif
      write(iLogFileUnit_,'(a,a15)') "# AMIE: ", cAmieFileNorth, cAmieFileSouth
      write(iLogFileUnit_,'(3(a,L2))') "# Solar Heating: ",useSolarHeating, &
          " Joule Heating: ",useJouleHeating, &
          " Auroral Heating: ", useAuroralHeating
       write(iLogFileUnit_,'(2(a,L2))') "# NO Cooling: ", useNOCooling, &
          " O Cooling: ", useOCooling
       write(iLogFileUnit_,'(3(a,L2))') "# Conduction: ",useConduction, &
          " Turbulent Conduction: ", useTurbulentCond
       write(iLogFileUnit_,'(3(a,L2))') "# Pressure Grad: ", &
            usePressureGradient, " Ion Drag: ", useIonDrag, &
          " Neutral Drag: ", useNeutralDrag
       write(iLogFileUnit_,'(3(a,L2))') "# Viscosity: ", useViscosity,&
          " Coriolis: ", useCoriolis, " Gravity: ",useGravity 
       write(iLogFileUnit_,'(3(a,L2))') "# Ion Chemistry: ", useIonChemistry, &
          " Ion Advection: ", useIonAdvection, " Neutral Chemistry: ", &
          useNeutralChemistry
          
       write(iLogFileUnit_,'(a)') " "
       write(iLogFileUnit_,'(a)') "#START"
       write(iLogFileUnit_,'(a)') &
            "   iStep yyyy mm dd hh mm ss  ms      dt "// &
            "min(T) max(T) mean(T) min(VV) max(VV) mean(VV) F107 F107A "// &
            "By Bz Vx HP HPn HPs SubsolarLon SubsolarLat SubsolarVTEC"
  endif

  call get_subsolar(CurrentTime, VernalTime, SSLon, SSLat)
  call get_log_info(SSLon, SSLat, MinTemp, MaxTemp, MinVertVel, MaxVertVel, &
       AverageTemp, AverageVertVel, TotalVolume, SSVTEC)

  localVar = TotalVolume
  call MPI_REDUCE(localVar, TotalVolume, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  localVar = MinTemp
  call MPI_REDUCE(localVar, minTemp, 1, MPI_REAL, MPI_MIN, &
       0, iCommGITM, iError)

  localVar = MaxTemp
  call MPI_REDUCE(localVar, maxTemp, 1, MPI_REAL, MPI_MAX, &
       0, iCommGITM, iError)

  localVar = MinVertVel
  call MPI_REDUCE(localVar, minVertVel, 1, MPI_REAL, MPI_MIN, &
       0, iCommGITM, iError)

  localVar = MaxVertVel
  call MPI_REDUCE(localVar, maxVertVel, 1, MPI_REAL, MPI_MAX, &
       0, iCommGITM, iError)

  LocalVar = AverageTemp
  call MPI_REDUCE(LocalVar, AverageTemp, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  LocalVar = AverageVertVel
  call MPI_REDUCE(LocalVar, AverageVertVel, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  LocalVar = HemisphericPowerNorth
  call MPI_REDUCE(LocalVar, HPn, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  LocalVar = HemisphericPowerSouth
  call MPI_REDUCE(LocalVar, HPs, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  LocalVar = SSVTEC
  call MPI_REDUCE(LocalVar, SSVTEC, 1, MPI_REAL, MPI_MAX, &
       0, iCommGITM, iError) 

  if (iProc == 0) then

     AverageTemp = AverageTemp / TotalVolume
     AverageVertVel = AverageVertVel / TotalVolume

     call get_f107(CurrentTime, f107, iError)
     call get_f107A(CurrentTime, f107A, iError)
     call get_IMF_By(CurrentTime, by, iError)
     call get_IMF_Bz(CurrentTime, bz, iError)
     call get_sw_v(CurrentTime, Vx, iError)
     call get_hpi(CurrentTime,Hpi,iError)

     if (Is1D) SSVTEC = -1.0
     
     write(iLogFileUnit_,"(i8,i5,5i3,i4,f8.4,6f13.5,8f9.1,10f10.5,10f10.5,10f8.3)") &
          iStep, iTimeArray, dt, minTemp, maxTemp, AverageTemp, &
          minVertVel, maxVertVel, AverageVertVel,&
          f107,f107A,By,Bz,Vx, Hpi, HPn/1.0e9, &
          HPs/1.0e9, SSLon, SSLat, SSVTEC

     call flush_unit(iLogFileUnit_)
  endif

end subroutine logfile



subroutine write_code_information(dir)

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi
  use ModIndices
  use ModIndicesInterfaces
  use ModIoUnit,    ONLY: io_unit_new
  use ModUtilities, ONLY: flush_unit
  use ModRCMR

  implicit none

  character (len=*), intent(in) :: dir

  integer, dimension(7) :: iTime
  integer :: i

  if (iProc == 0) then

     iCodeInfoFileUnit_ = io_unit_new()

     open(unit=iCodeInfoFileUnit_, &
          file=trim(dir)//"/run_information.txt",status="replace")

     write(iCodeInfoFileUnit_,*) "GITM2 Run Information"
     write(iCodeInfoFileUnit_,*) "---------------------"
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "nSpecies", nSpecies
     write(iCodeInfoFileUnit_,*) "nSpeciesTotal", nSpeciesTotal
     do i=1,nSpeciesTotal
        write(iCodeInfoFileUnit_,*) "Mass of Species ",cSpecies(i), Mass(i)/AMU
     enddo
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) "nIons", nIons
     write(iCodeInfoFileUnit_,*) "nIonsAdvect", nIonsAdvect
     do i=1,nIons
        write(iCodeInfoFileUnit_,*) "Ion Species ",cIons(i)
     enddo

     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) "Inputs from UAM.in:"
     write(iCodeInfoFileUnit_,*) ""

     call time_real_to_int(StartTime, iTime)
     write(iCodeInfoFileUnit_,*) "#TIMESTART"
     write(iCodeInfoFileUnit_,*) iTime(1)
     write(iCodeInfoFileUnit_,*) iTime(2)
     write(iCodeInfoFileUnit_,*) iTime(3)
     write(iCodeInfoFileUnit_,*) iTime(4)
     write(iCodeInfoFileUnit_,*) iTime(5)
     write(iCodeInfoFileUnit_,*) iTime(6)
     write(iCodeInfoFileUnit_,*) ""

     call time_real_to_int(EndTime, iTime)
     write(iCodeInfoFileUnit_,*) "#TIMEEND"
     write(iCodeInfoFileUnit_,*) iTime(1)
     write(iCodeInfoFileUnit_,*) iTime(2)
     write(iCodeInfoFileUnit_,*) iTime(3)
     write(iCodeInfoFileUnit_,*) iTime(4)
     write(iCodeInfoFileUnit_,*) iTime(5)
     write(iCodeInfoFileUnit_,*) iTime(6)
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#GRID"
     write(iCodeInfoFileUnit_,*) nBlocksLon
     write(iCodeInfoFileUnit_,*) nBlocksLat
     write(iCodeInfoFileUnit_,*) LatStart
     write(iCodeInfoFileUnit_,*) LatEnd
     write(iCodeInfoFileUnit_,*) LonStart
     write(iCodeInfoFileUnit_,*) LonEnd
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#DIFFUSION"
     write(iCodeInfoFileUnit_,*) UseDiffusion
     write(iCodeInfoFileUnit_,*) EddyDiffusionCoef
     write(iCodeInfoFileUnit_,*) EddyDiffusionPressure0
     write(iCodeInfoFileUnit_,*) EddyDiffusionPressure1
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#THERMALCONDUCTION"
     write(iCodeInfoFileUnit_,*) ThermalConduction_AO2
     write(iCodeInfoFileUnit_,*) ThermalConduction_AO
     write(iCodeInfoFileUnit_,*) ThermalConduction_s
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#PHOTOELECTRON"
     write(iCodeInfoFileUnit_,*) PhotoElectronHeatingEfficiency
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#NEUTRALHEATING"
     write(iCodeInfoFileUnit_,*) NeutralHeatingEfficiency
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#CFL"
     write(iCodeInfoFileUnit_,*) cfl
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#LIMITER"
     write(iCodeInfoFileUnit_,*) trim(TypeLimiter)
     write(iCodeInfoFileUnit_,*) BetaLimiter
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#F107"
     write(iCodeInfoFileUnit_,*) f107
     write(iCodeInfoFileUnit_,*) f107a
     write(iCodeInfoFileUnit_,*) ""

!     write(iCodeInfoFileUnit_,*) "#SOLARWIND"
!     write(iCodeInfoFileUnit_,*) bx
!     write(iCodeInfoFileUnit_,*) by
!     write(iCodeInfoFileUnit_,*) bz
!     write(iCodeInfoFileUnit_,*) vx
!     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#THERMO"
     write(iCodeInfoFileUnit_,*) UseSolarHeating
     write(iCodeInfoFileUnit_,*) UseJouleHeating
     write(iCodeInfoFileUnit_,*) UseAuroralHeating
     write(iCodeInfoFileUnit_,*) UseNOCooling
     write(iCodeInfoFileUnit_,*) UseOCooling
     write(iCodeInfoFileUnit_,*) UseConduction
     write(iCodeInfoFileUnit_,*) UseTurbulentCond
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#FORCING"
     write(iCodeInfoFileUnit_,*) UsePressureGradient
     write(iCodeInfoFileUnit_,*) UseIonDrag
     write(iCodeInfoFileUnit_,*) UseNeutralFriction
     write(iCodeInfoFileUnit_,*) UseViscosity
     write(iCodeInfoFileUnit_,*) UseCoriolis
     write(iCodeInfoFileUnit_,*) UseGravity
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#APEX"
     write(iCodeInfoFileUnit_,*) UseApex
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#AMIEFILES"
     write(iCodeInfoFileUnit_,*) trim(cAMIEFileNorth)
     write(iCodeInfoFileUnit_,*) trim(cAMIEFileSouth)
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#STATISTICALMODELSONLY"
     write(iCodeInfoFileUnit_,*) UseStatisticalModelsOnly
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#TIDES"
     write(iCodeInfoFileUnit_,*) UseMSISOnly
     write(iCodeInfoFileUnit_,*) UseMSISTides
     write(iCodeInfoFileUnit_,*) UseGSWMTides
     write(iCodeInfoFileUnit_,*) UseWACCMTides
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#DUSTDATA"
     write(iCodeInfoFileUnit_,*) UseDustDistribution
     write(iCodeInfoFileUnit_,*) trim(cDustFile)
     write(iCodeInfoFileUnit_,*) trim(cConrathFile)
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#DUST"
     write(iCodeInfoFileUnit_,*) tautot_temp
     write(iCodeInfoFileUnit_,*) Conrnu_temp
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#GSWMCOMP"
     do i=1,4
        write(iCodeInfoFileUnit_,*) UseGswmComp(i)
     enddo
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#USEPERTURBATION"
     write(iCodeInfoFileUnit_,*) UsePerturbation
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#DAMPING"
     write(iCodeInfoFileUnit_,*) UseDamping
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#GRAVITYWAVE"
     write(iCodeInfoFileUnit_,*) UseGravityWave
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#AURORAMODS"
     write(iCodeInfoFileUnit_,*) NormalizeAuroraToHP
     write(iCodeInfoFileUnit_,*) AveEFactor
     write(iCodeInfoFileUnit_,*) IsKappaAurora
     write(iCodeInfoFileUnit_,*) AuroraKappa
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#NEWELLAURORA"
     write(iCodeInfoFileUnit_,*) UseNewellAurora
     write(iCodeInfoFileUnit_,*) UseNewellAveraged
     write(iCodeInfoFileUnit_,*) UseNewellMono
     write(iCodeInfoFileUnit_,*) UseNewellWave
     write(iCodeInfoFileUnit_,*) DoNewellRemoveSpikes
     write(iCodeInfoFileUnit_,*) DoNewellAverage
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#OVATIONSME"
     write(iCodeInfoFileUnit_,*) UseOvationSME
     write(iCodeInfoFileUnit_,*) UseOvationSMEMono
     write(iCodeInfoFileUnit_,*) UseOvationSMEWave
     write(iCodeInfoFileUnit_,*) UseOvationSMEIon
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#IONLIMITS"
     write(iCodeInfoFileUnit_,*) MaxVParallel
     write(iCodeInfoFileUnit_,*) MaxEField
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#VERTICALSOURCES"
     write(iCodeInfoFileUnit_,*) MaximumVerticalVelocity
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#DYNAMO"
     write(iCodeInfoFileUnit_,*) UseDynamo
     write(iCodeInfoFileUnit_,*) DynamoHighLatBoundary
     write(iCodeInfoFileUnit_,*) nItersMax
     write(iCodeInfoFileUnit_,*) MaxResidual
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#IONFORCING"
     write(iCodeInfoFileUnit_,*) UseExB
     write(iCodeInfoFileUnit_,*) UseIonPressureGradient
     write(iCodeInfoFileUnit_,*) UseIonGravity
     write(iCodeInfoFileUnit_,*) UseNeutralDrag
     write(iCodeInfoFileUnit_,*) UseDynamo
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#CHEMISTRY"
     write(iCodeInfoFileUnit_,*) UseIonChemistry
     write(iCodeInfoFileUnit_,*) UseIonAdvection
     write(iCodeInfoFileUnit_,*) UseNeutralChemistry
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#FIXTILT"
     write(iCodeInfoFileUnit_,*) IsFixedTilt
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#DIPOLE"
     write(iCodeInfoFileUnit_,*) MagneticPoleRotation
     write(iCodeInfoFileUnit_,*) MagneticPoleTilt
     write(iCodeInfoFileUnit_,*) xDipoleCenter
     write(iCodeInfoFileUnit_,*) yDipoleCenter
     write(iCodeInfoFileUnit_,*) zDipoleCenter
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#APEX"
     write(iCodeInfoFileUnit_,*) UseApex
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#NEWSTRETCH"
     write(iCodeInfoFileUnit_,*) NewStretchedGrid
     write(iCodeInfoFileUnit_,*) ConcentrationLatitude
     write(iCodeInfoFileUnit_,*) StretchWidth
     write(iCodeInfoFileUnit_,*) StretchingPercentage
     write(iCodeInfoFileUnit_,*) ""
     
     write(iCodeInfoFileUnit_,*) "#STRETCH"
     write(iCodeInfoFileUnit_,*) ConcentrationLatitude
     write(iCodeInfoFileUnit_,*) StretchingPercentage
     write(iCodeInfoFileUnit_,*) StretchingFactor
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#TOPOGRAPHY"
     write(iCodeInfoFileUnit_,*) UseTopography
     write(iCodeInfoFileUnit_,*) AltMinUniform
     write(iCodeInfoFileUnit_,*) ""

     !write(iCodeInfoFileUnit_,*) "#RCMR"
     !write(iCodeInfoFileUnit_,*) RCMRInType
     !write(iCodeInfoFileUnit_,*) RCMROutType
     !write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#ELECTRODYNAMICS"
     write(iCodeInfoFileUnit_,*) dTPotential
     write(iCodeInfoFileUnit_,*) dTAurora
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#INPUTTIMEDELAY"
     write(iCodeInfoFileUnit_,*) TimeDelayHighLat
     write(iCodeInfoFileUnit_,*) TimeDelayEUV
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#LTERadiation"
     write(iCodeInfoFileUnit_,*) DtLTERadiation
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#EUV_DATA"
     write(iCodeInfoFileUnit_,*) UseEUVData
     write(iCodeInfoFileUnit_,*) trim(cEUVFile)
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#RESTART"
     write(iCodeInfoFileUnit_,*) DoRestart
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "UseCO2Cooling", UseCO2Cooling
     write(iCodeInfoFileUnit_,*) "CO2ppm", CO2ppm
     write(iCodeInfoFileUnit_,*) "iInputIonChemType", iInputIonChemType
     write(iCodeInfoFileUnit_,*) "iInputNeutralChemType", iInputNeutralChemType
     write(iCodeInfoFileUnit_,*) "RPTAU", RPTAU
     write(iCodeInfoFileUnit_,*) "PRAD", PRAD
     write(iCodeInfoFileUnit_,*) "PLONG", PLONG
     write(iCodeInfoFileUnit_,*) "Ls_a", Ls_a
     write(iCodeInfoFileUnit_,*) "Ls_tau", Ls_tau
     write(iCodeInfoFileUnit_,*) "Ls_phi", Ls_phi
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""

     close(iCodeInfoFileUnit_)

  endif

end subroutine write_code_information
