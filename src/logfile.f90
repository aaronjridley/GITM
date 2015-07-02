!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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
  use ModUtilities, ONLY: flush_unit

  implicit none

  character (len=*), intent(in) :: dir
  character (len=8) :: cIter

  real    :: minTemp, maxTemp, localVar, minVertVel, maxVertVel
  real    :: AverageTemp, AverageVertVel, TotalVolume, Bx, By, Bz, Vx, Hpi
  real    :: HPn, HPs, SSLon, SSLat, SSVTEC
  integer :: iError

  if (.not. IsOpenLogFile .and. iProc == 0) then

     IsOpenLogFile = .true.
     call CON_io_unit_new(iLogFileUnit_)

     write(cIter,"(i8.8)") iStep

     open(unit=iLogFileUnit_, &
          file=dir//"/log"//cIter//".dat",status="replace")

     write(iLogFileUnit_,'(a)') "GITM2 log file"
     write(iLogFileUnit_,'(a,L2)') "## Inputs from UAM.in" 
      write(iLogFileUnit_,'(a,L2)') "# Resart=", dorestart
     write(iLogFileUnit_,'(4(a,f9.3))') "# Eddy coef: ", EddyDiffusionCoef, &
          " Eddy P0: ",EddyDiffusionPressure0,&
          " Eddy P1: ",EddyDiffusionPressure1,&
          " Eddy Scaling: ",EddyScaling
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
          " Turbulent Conduction: ", useTurbulentCond, &
          " Updated Turbulent Conduction: ",useUpdatedTurbulentCond
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
  use ModUtilities, ONLY: flush_unit

  implicit none

  character (len=*), intent(in) :: dir

  integer, dimension(7) :: iTime

  if (iProc == 0) then

     call CON_io_unit_new(iCodeInfoFileUnit_)

     open(unit=iCodeInfoFileUnit_, &
          file=dir//"/run_information.txt",status="replace")

     write(iCodeInfoFileUnit_,*) "GITM2 Run Information"
     write(iCodeInfoFileUnit_,*) "---------------------"
     write(iCodeInfoFileUnit_,*) "Inputs from UAM.in:"
     write(iCodeInfoFileUnit_,*) ""

     call time_real_to_int(StartTime, iTime)
     write(iCodeInfoFileUnit_,*) "#TIMESTART"
     write(iCodeInfoFileUnit_,*) iTimeArray(1)
     write(iCodeInfoFileUnit_,*) iTimeArray(2)
     write(iCodeInfoFileUnit_,*) iTimeArray(3)
     write(iCodeInfoFileUnit_,*) iTimeArray(4)
     write(iCodeInfoFileUnit_,*) iTimeArray(5)
     write(iCodeInfoFileUnit_,*) iTimeArray(6)
     write(iCodeInfoFileUnit_,*) ""

     call time_real_to_int(EndTime, iTime)
     write(iCodeInfoFileUnit_,*) "#TIMEEND"
     write(iCodeInfoFileUnit_,*) iTimeArray(1)
     write(iCodeInfoFileUnit_,*) iTimeArray(2)
     write(iCodeInfoFileUnit_,*) iTimeArray(3)
     write(iCodeInfoFileUnit_,*) iTimeArray(4)
     write(iCodeInfoFileUnit_,*) iTimeArray(5)
     write(iCodeInfoFileUnit_,*) iTimeArray(6)
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
     write(iCodeInfoFileUnit_,*) PhotoElectronHeatingEfficiency_est
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#CFL"
     write(iCodeInfoFileUnit_,*) cfl
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#F107"
     write(iCodeInfoFileUnit_,*) f107
     write(iCodeInfoFileUnit_,*) f107a
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) "#THERMO"
     write(iCodeInfoFileUnit_,*) UseSolarHeating
     write(iCodeInfoFileUnit_,*) UseJouleHeating
     write(iCodeInfoFileUnit_,*) UseAuroralHeating
     write(iCodeInfoFileUnit_,*) UseNOCooling
     write(iCodeInfoFileUnit_,*) UseOCooling
     write(iCodeInfoFileUnit_,*) UseConduction
     write(iCodeInfoFileUnit_,*) UseTurbulentCond
     write(iCodeInfoFileUnit_,*) UseUpdatedTurbulentCond
     write(iCodeInfoFileUnit_,*) EddyScaling
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
     write(iCodeInfoFileUnit_,*) cAMIEFileNorth
     write(iCodeInfoFileUnit_,*) cAMIEFileSouth
     write(iCodeInfoFileUnit_,*) ""

     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""
     write(iCodeInfoFileUnit_,*) ""


     write(iCodeInfoFileUnit_,'(a,L2)') "# Resart=", dorestart
     write(iCodeInfoFileUnit_,'(4(a,f9.3))') "# Eddy coef: ", EddyDiffusionCoef, &
          " Eddy P0: ",EddyDiffusionPressure0,&
          " Eddy P1: ",EddyDiffusionPressure1,&
          " Eddy Scaling: ",EddyScaling
     write(iCodeInfoFileUnit_,'(2(a,L2))') "# Statistical Models Only: ", &
          usestatisticalmodelsonly, " Apex: ",useApex
     if (useEUVdata) then
        write(iCodeInfoFileUnit_,'(a,L2,a)') "# EUV Data: ", useEUVdata, "File: ", &
             cEUVFile
     else
        write(iCodeInfoFileUnit_,'(a,L2)') "# EUV Data: ", useEUVdata
     endif
     write(iCodeInfoFileUnit_,'(a,a15)') "# AMIE: ", cAmieFileNorth, cAmieFileSouth
     write(iCodeInfoFileUnit_,'(3(a,L2))') "# Solar Heating: ",useSolarHeating, &
          " Joule Heating: ",useJouleHeating, &
          " Auroral Heating: ", useAuroralHeating
     write(iCodeInfoFileUnit_,'(2(a,L2))') "# NO Cooling: ", useNOCooling, &
          " O Cooling: ", useOCooling
     write(iCodeInfoFileUnit_,'(3(a,L2))') "# Conduction: ",useConduction, &
          " Turbulent Conduction: ", useTurbulentCond, &
          " Updated Turbulent Conduction: ",useUpdatedTurbulentCond
     write(iCodeInfoFileUnit_,'(3(a,L2))') "# Pressure Grad: ", &
          usePressureGradient, " Ion Drag: ", useIonDrag, &
          " Neutral Drag: ", useNeutralDrag
     write(iCodeInfoFileUnit_,'(3(a,L2))') "# Viscosity: ", useViscosity,&
          " Coriolis: ", useCoriolis, " Gravity: ",useGravity 
     write(iCodeInfoFileUnit_,'(3(a,L2))') "# Ion Chemistry: ", useIonChemistry, &
          " Ion Advection: ", useIonAdvection, " Neutral Chemistry: ", &
          useNeutralChemistry

  endif

end subroutine write_code_information
