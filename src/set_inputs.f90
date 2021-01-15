!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-----------------------------------------------------------------------------
! $Id: set_inputs.f90,v 1.84 2017/10/30 14:05:36 ridley Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Comments: Routine to read  input options from the UAM.in file.
!
! AGB 3/31/13: Added flags for data assimilation to adjust F10.7
! AGB 10/23/13: Changed RCMR input to allow driving of photoelectron heating
!               efficiency
! AGB 12/20/13: Added input options for efield and ion diff. vel. limits
!-----------------------------------------------------------------------------

subroutine set_inputs

  use ModInputs
  use ModSizeGitm
  use ModGITM, only: iProc, f107_est, f107a_est, f107_msis, f107a_msis, &
       PhotoElectronHeatingEfficiency_est, EDC_est !! Ankit23May16: Added EDC_est
  use ModTime
  use ModPlanet
  use ModSatellites
  use ModRCMR

  implicit none

  integer, external :: bad_outputtype
  integer, external :: jday

  integer, dimension(7) :: iEndTime

  logical :: IsDone
  integer :: iDebugProc=0
  character (len=iCharLen_) :: cLine
  integer :: iLine, iSpecies, iSat
  integer :: i, iError, iOutputTypes
  integer, dimension(7) :: iTimeEnd

  character (len=iCharLen_)                 :: cTempLine
  character (len=iCharLen_)                 :: sIonChemistry, sNeutralChemistry
  character (len=iCharLen_), dimension(100) :: cTempLines

  real :: Vx, Bx, Bz, By, Kp, HemisphericPower, tsim_temp
  real :: EDC_est_tmp
  real*8 :: DTime

  call report("set_inputs",1)

  iError = 0
  IsDone = .false.
  iLine  = 1

  call IO_set_ap_single(10.0)

  do while (.not. IsDone)

     cLine = cInputText(iLine)

     if (cLine(1:1) == "#") then

        ! Remove anything after a space or TAB
        i=index(cLine,' '); if(i>0)cLine(i:len(cLine))=' '
        i=index(cLine,char(9)); if(i>0)cLine(i:len(cLine))=' '

        if (iDebugLevel > 3) write(*,*) "====> cLine : ",cLine(1:40)

        select case (cLine)

        case ("#TIMESTART","#STARTTIME")

           if (IsFramework) then

              if (iDebugLevel > 0) then
                 write(*,*) "----------------------------------"
                 write(*,*) "-   UAM trying to set STARTTIME  -"
                 write(*,*) "-          Ignoring              -"
                 write(*,*) "----------------------------------"
              endif

           else

              iStartTime = 0
              do i=1,6
                 call read_in_int(iStartTime(i), iError)
              enddo

              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #STARTTIME:'
                 write(*,*) 'This sets the starting time of the simulation.'
                 write(*,*) 'Even when you restart, the starttime should be'
                 write(*,*) 'to the real start time, not the restart time.'
                 write(*,*) '#STARTTIME'
                 write(*,*) 'iYear    (integer)'
                 write(*,*) 'iMonth   (integer)'
                 write(*,*) 'iDay     (integer)'
                 write(*,*) 'iHour    (integer)'
                 write(*,*) 'iMinute  (integer)'
                 write(*,*) 'iSecond  (integer)'
                 IsDone = .true.
              else

                 iTimeArray = iStartTime
                 call time_int_to_real(iStartTime, CurrentTime)
                 StartTime = CurrentTime
                 if (tSimulation > 0) then
                    CurrentTime = CurrentTime + tSimulation
                    call time_real_to_int(CurrentTime, iTimeArray)
                 endif

                 call fix_vernal_time

              endif

           endif

        case ("#TIMEEND","#ENDTIME")

           if (IsFramework) then
              if (iDebugLevel > 0) then
                 write(*,*) "--------------------------------"
                 write(*,*) "-   UAM trying to set ENDTIME  -"
                 write(*,*) "-          Ignoring            -"
                 write(*,*) "--------------------------------"
              endif
           else

              iTimeEnd = 0
              do i=1,6
                 call read_in_int(iTimeEnd(i), iError)
              enddo

              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #ENDTIME:'
                 write(*,*) 'This sets the ending time of the simulation.'
                 write(*,*) '#ENDTIME'
                 write(*,*) 'iYear    (integer)'
                 write(*,*) 'iMonth   (integer)'
                 write(*,*) 'iDay     (integer)'
                 write(*,*) 'iHour    (integer)'
                 write(*,*) 'iMinute  (integer)'
                 write(*,*) 'iSecond  (integer)'
                 IsDone = .true.
              else
                 call time_int_to_real(iTimeEnd, EndTime)
                 PauseTime = EndTime + 1.0
              endif

           endif

        case ("#PAUSETIME")
           call read_in_time(PauseTime, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #PAUSE:'
              write(*,*) 'This will set a time for the code to just pause.'
              write(*,*) 'Really, this should never be used.'
              write(*,*) '#PAUSETIME'
              write(*,*) 'iYear iMonth iDay iHour iMinute iSecond'
              IsDone = .true.
           endif

        case ("#ISTEP")
           call read_in_int(iStep, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #ISTEP:'
              write(*,*) 'This is typically only specified in a'
              write(*,*) 'restart header.  If you specify it in a start UAM.in'
              write(*,*) 'it will start the counter to whatever you specify.'
              write(*,*) '#ISTEP'
              write(*,*) 'iStep     (integer)'
              IsDone = .true.
           endif

        case ("#CPUTIMEMAX")
           call read_in_real(CPUTIMEMAX, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #CPUTIMEMAX:'
              write(*,*) 'This sets the maximum CPU time that the code should'
              write(*,*) 'run before it starts to write a restart file and end'
              write(*,*) 'the simulation.  It is very useful on systems that'
              write(*,*) 'have a queueing system and has limited time runs.'
              write(*,*) 'Typically, set it for a couple of minutes short of'
              write(*,*) 'the max wall clock, since it needs some time to write'
              write(*,*) 'the restart files.'
              write(*,*) '#CPUTIMEMAX'
              write(*,*) 'CPUTimeMax    (real)'
           endif

        case ("#STATISTICALMODELSONLY")
           call read_in_logical(UseStatisticalModelsOnly, iError)
           if (UseStatisticalModelsOnly .and. iError == 0) &
                call read_in_real(DtStatisticalModels,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #STATISTICALMODELSONLY:'
              write(*,*) 'This command will skip all pretty much all of the '
              write(*,*) 'physics of GITM, and will reinitialize the model '
              write(*,*) 'with the MSIS and IRI values at the interval set in '
              write(*,*) 'the second variable. If you want to compare a run to' 
              write(*,*) 'MSIS and IRI, you can run GITM with this command and'
              write(*,*) 'get output at exactly the same cadence and locations,'
              write(*,*) 'thereby allowing easier comparisons. The dt can be '
              write(*,*) 'set as low as needed, so you can run satellites '
              write(*,*) 'through MSIS and IRI.'
              write(*,*) '#STATISTICALMODELSONLY'
              write(*,*) 'UseStatisticalModelsOnly    (logical)'
              write(*,*) 'DtStatisticalModels         (real)'
           endif

        case ("#TSIMULATION")

           if (IsFramework) then
              if (iDebugLevel > 0) then
                 write(*,*) "------------------------------------"
                 write(*,*) "-   UAM trying to set tsimulation  -"
                 write(*,*) "-          Ignoring                -"
                 write(*,*) "------------------------------------"
              endif
           else
              call read_in_real(tSim_temp, iError)
              tSimulation = tSim_temp
              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #TSIMULATION:'
                 write(*,*) 'This is typically only specified in a'
                 write(*,*) 'restart header.'
                 write(*,*) 'It sets the offset from the starttime to the'
                 write(*,*) 'currenttime. Should really only be used with caution.'
                 write(*,*) '#TSIMULATION'
                 write(*,*) 'tsimulation    (real)'
                 IsDone = .true.
              else
                 CurrentTime = CurrentTime + tsimulation
                 call time_real_to_int(CurrentTime, iTimeArray)
              endif
           endif

        case ("#F107")
           if(RCMROutType == "F107") then
              call read_in_real(f107_msis, iError)
              call read_in_real(f107a_msis, iError)
           else
              call read_in_real(f107, iError)
              call read_in_real(f107a, iError)
           end if

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #F107:'
              write(*,*) 'Sets the F10.7 and 81 day average F10.7.  This is'
              write(*,*) 'used to set the initial altitude grid, and drive the'
              write(*,*) 'lower boundary conditions.'
              write(*,*) '#F107'
              write(*,*) 'f107  (real)'
              write(*,*) 'f107a (real - 81 day average of f107)'
              if(RCMROutType == "F107") then
                 f107_msis  = 150.0
                 f107a_msis = 150.0
              else
                 f107  = 150.0
                 f107a = 150.0
              end if
              IsDone = .true.
           endif

           ! AGB: Added switch for initializing the GITM f107 depending on
           !      the data assimilation (RCMR) flags

           if(RCMROutType /= "F107") then
              call IO_set_f107_single(f107)
              call IO_set_f107a_single(f107a)
           end if

        case ("#INITIAL")

           call read_in_logical(UseMSIS, iError)
           call read_in_logical(UseIRI, iError)
           if (.not.UseMSIS) then
              call read_in_real(TempMin,    iError)
              call read_in_real(TempMax,    iError)
              call read_in_real(TempHeight, iError)
              call read_in_real(TempWidth,  iError)
              do i=1,nSpecies
                 call read_in_real(LogNS0(i), iError)
                 if (iError == 0) LogNS0(i) = alog(LogNS0(i))
              enddo
              LogRho0 = 0.0
              do iSpecies = 1, nSpecies
                 LogRho0 = LogRho0 + &
                      exp(LogNS0(iSpecies)) * Mass(iSpecies)
              enddo
              LogRho0 = alog(LogRho0)
           endif
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #INITIAL:'
              write(*,*) 'This specifies the initial conditions and the'
              write(*,*) 'lower boundary conditions.  For Earth, we typically'
              write(*,*) 'just use MSIS and IRI for initial conditions.'
              write(*,*) 'For other planets, the vertical BCs can be set here.'
              write(*,*) '#INITIAL'
              write(*,*) 'UseMSIS        (logical)'
              write(*,*) 'UseIRI         (logical)'
              write(*,*) 'If UseMSIS is .false. then :'
              write(*,*) 'TempMin        (real, bottom temperature)'
              write(*,*) 'TempMax        (real, top initial temperature)'
              write(*,*) 'TempHeight     (real, Height of the middle of temp gradient)'
              write(*,*) 'TempWidth      (real, Width of the temperature gradient)'
              do i=1,nSpecies
                 write(*,*) 'Bottom N Density (real), species',i
              enddo
           endif

        case ("#TIDES")
           call read_in_logical(UseMSISOnly, iError)
           call read_in_logical(UseMSISTides, iError)
           call read_in_logical(UseGSWMTides, iError)
           call read_in_logical(UseWACCMTides, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #TIDES:'
              write(*,*) 'This says how to use tides.  The first one is using'
              write(*,*) 'MSIS with no tides.  The second is using MSIS with'
              write(*,*) 'full up tides. The third is using GSWM tides, while'
              write(*,*) 'the forth is for experimentation with using WACCM'
              write(*,*) 'tides.'
              write(*,*) '#TIDES'
              write(*,*) 'UseMSISOnly        (logical)'
              write(*,*) 'UseMSISTides       (logical)'
              write(*,*) 'UseGSWMTides       (logical)'
              write(*,*) 'UseWACCMTides      (logical)'
           else
              if (UseGSWMTides)  UseMSISOnly = .true.
              if (UseWACCMTides) UseMSISOnly = .true.
           endif

        case ("#MSISTIDES")
           call read_in_logical(UseMSISDiurnal, iError)
           call read_in_logical(UseMSISSemidiurnal, iError)
           call read_in_logical(UseMSISTerdiurnal, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #MSISTIDES:'
              write(*,*) 'This says how to use msis tides.  '
              write(*,*) 'The first one is using diurnal tide'
              write(*,*) 'The first one is using semi-diurnal tide'
              write(*,*) 'The first one is using terdiurnal tide'
              write(*,*) '#MSISTIDES'
              write(*,*) 'UseMSISDiurnal        (logical)'
              write(*,*) 'UseMSISSemidiurnal    (logical)'
              write(*,*) 'UseMSISTerdiurnal     (logical)'
           endif

        case ("#MSISOBC")
           call read_in_logical(UseOBCExperiment, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #MSISOBC:'
              write(*,*) '#MSISOBC'
              write(*,*) 'UseOBCExperiment        (logical)'
           endif

        !xianjing
        case("#USESECONDSINFILENAME")
           call read_in_logical(UseSecondsInFilename,iError)
            if (iError /= 0) then
              write(*,*) 'Incorrect format for #USESECONDSINFILENAME:'
              write(*,*) 'This says whether you want seconds in output file'
              write(*,*) 'name. F means no seconds in output file name.'
           else
              if (UseSecondsInFilename)  UseSecondsInFilename = .false.
           endif

        case ("#DUSTDATA")
           call read_in_logical(UseDustDistribution,iError) 
           call read_in_string(cDustFile,iError)
           call read_in_string(cConrathFile,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DUSTDATA'
              write(*,*) '#DUST'
              write(*,*) 'cDustFile'
              write(*,*) 'cConrathFile'
           endif

        case ("#OVERWRITEIONOSPHERE")
           call read_in_logical(DoOverwriteIonosphere,iError)
           call read_in_logical(DoOverwriteWithIRI,iError)
           call read_in_logical(DoOverwriteWithSami,iError)
           if (DoOverwriteWithSami) &
                call read_in_string(SamiInFile,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #OVERWRITEIONOSPHERE'
              write(*,*) '#OVERWRITEIONOSPHERE'
              write(*,*) 'DoOverwriteIonosphere'
              write(*,*) 'DoOverwriteWithIRI'
              write(*,*) 'DoOverwriteWithSami'
              write(*,*) 'SamiInFile'
           endif

        case ("#GITMBCS")
           call read_in_logical(UseGitmBCs,iError)
           call read_in_string(GitmBCsDir,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #GITMBCS'
              write(*,*) '#GITMBCS'
              write(*,*) 'UseGitmBCs'
              write(*,*) 'GitmBCsDir'
           endif

        case ("#DUST")
           call read_in_real(tautot_temp,iError)
           call read_in_real(Conrnu_temp,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DUST'
              write(*,*) '#DUST'
              write(*,*) 'TauTot'
              write(*,*) 'Conrnu'
           endif

        CASE ("#GSWMCOMP")
           DO i=1,4
              CALL read_in_logical(UseGswmComp(i), iError)
           enddo
           IF (iError /= 0) then
              write(*,*) 'Incorrect format for #GSWMCOMP:'
              write(*,*) 'If you selected to use GSWM tides above, you'
              write(*,*) 'can specify which components to use.'
              write(*,*) '#GSWMCOMP'
              write(*,*) 'GSWMdiurnal(1)        (logical)'
              write(*,*) 'GSWMdiurnal(2)        (logical)'
              write(*,*) 'GSWMsemidiurnal(1)    (logical)'
              write(*,*) 'GSWMsemidiurnal(2)    (logical)'
           ENDIF

        case ("#USEPERTURBATION")
           call read_in_logical(UsePerturbation,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #USEPERTURBATION:'
              write(*,*) ''
              write(*,*) ''
              write(*,*) '#USEPERTURBATION'
              write(*,*) 'UsePerturbation        (logical)'
              IsDone = .true.
           endif

        case ("#DAMPING")
           call read_in_logical(UseDamping,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DAMPING:'
              write(*,*) 'This is probably for damping vertical wind'
              write(*,*) 'oscillations that can occur in the lower atmosphere.'
              write(*,*) '#DAMPING'
              write(*,*) 'UseDamping        (logical)'
              IsDone = .true.
           endif

        case ("#GRAVITYWAVE")
           call read_in_logical(UseGravityWave,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #GRAVITYWAVE:'
              write(*,*) 'I dont know what this is for...'
              write(*,*) ''
              write(*,*) '#GRAVITYWAVE'
              write(*,*) 'UseGravityWave        (logical)'
              IsDone = .true.
           endif

        case ("#HPI")
           call read_in_real(HemisphericPower, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #HPI:'
              write(*,*) 'This sets the hemispheric power of the aurora.'
              write(*,*) 'Typical it ranges from 1-1000, although 20 is a'
              write(*,*) 'nominal, quiet time value.'
              write(*,*) '#HPI'
              write(*,*) 'HemisphericPower  (real)'
              IsDone = .true.
           else
              call IO_set_hpi_single(HemisphericPower)
           endif

        case ("#KP")
           call read_in_real(kp, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #KP:'
              write(*,*) 'I dont think that GITM actually uses this unless'
              write(*,*) 'the Foster electric field model is used.'
              write(*,*) '#KP'
              write(*,*) 'kp  (real)'
              IsDone = .true.
           else
              call IO_set_kp_single(kp)
           endif

        case ("#CFL")
           call read_in_real(cfl, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #cfl:'
              write(*,*) 'The CFL condition sets how close to the maximum time'
              write(*,*) 'step that GITM will take.  1.0 is the maximum value.'
              write(*,*) 'A value of about 0.75 is typical.  If instabilities'
              write(*,*) 'form, a lower value is probably needed.'
              write(*,*) '#CFL'
              write(*,*) 'cfl  (real)'
              IsDone = .true.
           endif

        case ("#FIXEDDT")
           call read_in_real(fixeddt, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #FIXEDDT:'
              write(*,*) 'If you would like to force GITM to take a fixed dt'
              write(*,*) 'you can use this.  It will try to take that fixed dt,'
              write(*,*) 'unless the CFL condition is violated.'
              write(*,*) '#FIXEDDT'
              write(*,*) 'FixedDt  (real)'
              IsDone = .true.
           endif

        case ("#SOLARWIND")
           call read_in_real(bx, iError)
           call read_in_real(by, iError)
           call read_in_real(bz, iError)
           call read_in_real(vx, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #SOLARWIND:'
              write(*,*) 'This sets the driving conditions for the high-latitude'
              write(*,*) 'electric field models.  This is static for the whole'
              write(*,*) 'run, though.  It is better to use the MHD\_INDICES'
              write(*,*) 'command to have dynamic driving conditions.'
              write(*,*) '#SOLARWIND'
              write(*,*) 'bx  (real)'
              write(*,*) 'by  (real)'
              write(*,*) 'bz  (real)'
              write(*,*) 'vx  (real)'
              IsDone = .true.
           else
              call IO_set_imf_by_single(by)
              call IO_set_imf_bz_single(bz)
              call IO_set_sw_v_single(abs(vx))
           endif

        case ("#MHD_INDICES")

           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #MHD_INDICES:'
              write(*,*) 'Use this for dynamic IMF and solar wind conditions.'
              write(*,*) 'The exact format of the file is discussed further'
              write(*,*) 'in the manual.'
              write(*,*) '#MHD_INDICES'
              write(*,*) 'filename  (string)'
              IsDone = .true.
           else

              cTempLines(2) = cTempLine
              cTempLines(3) = " "
              cTempLines(4) = "#END"

              iError = 0
              call IO_set_inputs(cTempLines)

              call read_MHDIMF_Indices_new(iError, CurrentTime, EndTime)
              if (iError /= 0) then 
                 write(*,*) "read indices was NOT successful (imf file)"
                 IsDone = .true.
              else
                 UseVariableInputs = .true.
              endif

           endif

        case ("#AURORAMODS")
           call read_in_logical(NormalizeAuroraToHP, iError)
           call read_in_real(AuroralHeightFactor, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #AURORAMODS'
              write(*,*) 'This is for modifying the aurora a bit.  The'
              write(*,*) 'NormalizeAuroraToHP variable calculates the '
              write(*,*) 'modeled hemispheric power and then normalizes it'
              write(*,*) 'the hemispheric power read in.  AuroralHeightFactor'
              write(*,*) 'adjusts the pressure to that the ionization height'
              write(*,*) 'can be at the right height.  A value of 1 is the'
              write(*,*) 'default height.  Lowering it to 0.4 changes the'
              write(*,*) 'height to be consistent with FLIP (but that was a'
              write(*,*) 'while ago that I checked this.'
              write(*,*) ''
              write(*,*) '#AURORAMODS'
              write(*,*) 'NormalizeAuroraToHP     (logical)'
              write(*,*) 'AuroralHeightFactor     (real)'
           endif

        case ("#NEWELLAURORA")
           call read_in_logical(UseNewellAurora, iError)
           call read_in_logical(UseNewellAveraged, iError)
           call read_in_logical(UseNewellMono, iError)
           call read_in_logical(UseNewellWave, iError)
           call read_in_logical(DoNewellRemoveSpikes, iError)
           call read_in_logical(DoNewellAverage, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #NEWELLAURORA'
              write(*,*) 'This is for using Pat Newells aurora (Ovation).'
              write(*,*) ''
              write(*,*) '#NEWELLAURORA'
              write(*,*) 'UseNewellAurora   (logical)'
              write(*,*) 'UseNewellAveraged (logical)'
              write(*,*) 'UseNewellMono (logical)'
              write(*,*) 'UseNewellWave (logical)'
              write(*,*) 'UseNewellRemoveSpikes (logical)'
              write(*,*) 'UseNewellAverage      (logical)'
              IsDone = .true.
           endif

        case ("#OVATIONSME")
           call read_in_logical(UseOvationSME, iError)
           call read_in_logical(UseOvationSMEMono, iError)
           call read_in_logical(UseOvationSMEWave, iError)
           call read_in_logical(UseOvationSMEIon, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #OVATIONSME'
              write(*,*) 'This is for using Betsy Michells aurora (OvationSME).'
              write(*,*) ''
              write(*,*) '#OVATIONSME'
              write(*,*) 'UseOvationSME     (logical)'
              write(*,*) 'UseOvationSMEMono (logical)'
              write(*,*) 'UseOvationSMEWave (logical)'
              write(*,*) 'UseOvationSMEIon  (logical)'
              IsDone = .true.
           endif

        case ("#AEMODEL")
           call read_in_logical(UseAeModel, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #AEMODEL'
              write(*,*) 'This is for using Dongjies aurora.'
              write(*,*) ''
              write(*,*) '#AEMODEL'
              write(*,*) 'UseAeModel        (logical)'
              IsDone = .true.
           endif

        case ("#FANGENERGY")
           call read_in_logical(UseFangEnergyDeposition, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #FANGENERGY'
              write(*,*) 'This is for using Dongjies aurora.'
              write(*,*) ''
              write(*,*) '#FANGENERGY'
              write(*,*) 'UseFangEnergyDeposition        (logical)'
              IsDone = .true.
           endif

        case ("#USECUSP")
           call read_in_logical(UseCusp, iError)
           call read_in_real(CuspAveE, iError)
           call read_in_real(CuspEFlux, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #USECUSP'
              write(*,*) 'This is for specifying a cusp.'
          write(*,*) ''
              write(*,*) '#USECUSP'
              write(*,*) 'UseCusp        (logical)'
              write(*,*) 'CuspAveE       (real)'
              write(*,*) 'CuspEFlux      (real)'
              IsDone = .true.
           endif

        case ("#AMIEFILES")
           call read_in_string(cAMIEFileNorth, iError)
           call read_in_string(cAMIEFileSouth, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #AMIEFILES:'
              write(*,*) ''
              write(*,*) '#AMIEFILES'
              write(*,*) 'cAMIEFileNorth  (string)'
              write(*,*) 'cAMIEFileSouth  (string)'
              IsDone = .true.
           endif

        case ("#USEREGIONALAMIE")
           call read_in_logical(UseRegionalAMIE, iError)
           call read_in_logical(UseTwoAMIEPotentials, iError)
           call read_in_time(AMIETimeStart, iError)
           call read_in_time(AMIETimeEnd, iError)
           call read_in_real(AMIELonStart, iError)
           call read_in_real(AMIELonEnd, iError)
           call read_in_real(AMIELatStart, iError)
           call read_in_real(AMIELatEnd, iError)
           call read_in_real(AMIEBoundaryWidth, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #USEREGIONALAMIE:'
              write(*,*) 'This is to set up a local region with specified potential'
              write(*,*) 'from AMIE files. Use Weimer potential elsewhere.'
              write(*,*) 'AMIEBoundaryWidth is padded outside of the region with '
              write(*,*) 'the geographic lon and lat boundaries set below.'
              write(*,*) ''
              write(*,*) '#USEREGIONALAMIE'
              write(*,*) 'UseRegionalAMIE      (logical)'
              write(*,*) 'UseTwoAMIEPotentials (logical)'
              write(*,*) 'AMIETimeStart        (yyyy mm dd hh mm ss)'
              write(*,*) 'AMIETimeEnd          (yyyy mm dd hh mm ss)'
              write(*,*) 'AMIELonStart         (real)'
              write(*,*) 'AMIELonEnd           (real)'
              write(*,*) 'AMIELatStart         (real)'
              write(*,*) 'AMIELatEnd           (real)'
              write(*,*) 'AMIEBoundaryWidth    (real)' 
              IsDone = .true.
           endif

        case ("#LIMITER")
           call read_in_string(TypeLimiter, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #LIMITER:'
              write(*,*) 'The limiter is quite important.  It is a value'
              write(*,*) 'between 1.0 and 2.0, with 1.0 being more diffuse and'
              write(*,*) 'robust, and 2.0 being less diffuse, but less robust.'
              write(*,*) '#LIMITER'
              write(*,*) 'TypeLimiter  (string)'
              write(*,*) 'BetaLimiter  (real between 1.0-minmod and 2.0-mc)'
              IsDone = .true.
           endif

           call read_in_real(BetaLimiter, iError)
           if (iError /= 0) then
              BetaLimiter = 1.6
              if (TypeLimiter == "minmod") BetaLimiter = 1.0
              if (iProc == 0) then
                 write(*,*) "You can now set the limiter value yourself!"
                 write(*,*) '#LIMITER'
                 write(*,*) 'TypeLimiter  (string)'
                 write(*,*) 'BetaLimiter  (real between 1.0-minmod and 2.0-mc)'
              endif
              iError = 0
           else
              if (BetaLimiter < 1.0) BetaLimiter = 1.0
              if (BetaLimiter > 2.0) BetaLimiter = 2.0
              if (iDebugLevel > 2) &
                   write(*,*) "===>Beta Limiter set to ",BetaLimiter
           endif

        case ("#NANCHECK")
           call read_in_logical(DoCheckForNans, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #NANCHECK:'
              write(*,*) 'This will turn on all of the NaN checks in the code!'
              write(*,*) '#NANCHECK'
              write(*,*) 'DoCheckForNans (logical)'
              IsDone = .true.
           endif

        case ("#DEBUG")
           call read_in_int(iDebugLevel, iError)
           call read_in_int(iDebugProc, iError)
           call read_in_real(DtReport, iError)
           call read_in_logical(UseBarriers, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DEBUG:'
              write(*,*) 'This will set how much information the code screams'
              write(*,*) 'at you - set to 0 to get minimal, set to 10 to get'
              write(*,*) 'EVERYTHING.  You can also change which processor is'
              write(*,*) 'shouting the information - PE 0 is the first one.'
              write(*,*) 'If you set the iDebugLevel to 0, you can set the dt'
              write(*,*) 'of the reporting.  If you set it to a big value,'
              write(*,*) 'you wont get very many outputs.  If you set it to a'
              write(*,*) 'tiny value, you will get a LOT of outputs.'
              write(*,*) 'UseBarriers will force the code to sync up a LOT.'
              write(*,*) '#DEBUG'
              write(*,*) 'iDebugLevel (integer)'
              write(*,*) 'iDebugProc  (integer)'
              write(*,*) 'DtReport    (real)'
              write(*,*) 'UseBarriers (logical)'
              IsDone = .true.
           endif

        case ("#IONLIMITS")
           call read_in_real(MaxVParallel, iError)
           call read_in_real(MaxEField, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #IONLIMITS:'
              write(*,*) ''
              write(*,*) '#IONLIMITS'
              write(*,*) "MaxVParallel   (real, default=100 m/s)"
              write(*,*) "MaxEField      (real, default=0.1 V/m)"
              MaxVParallel = 100.0
              MaxEField    = 0.1
              IsDone       = .true.
           end if

        case ("#PHOTOELECTRON")
           if(RCMROutType /= 'PHOTOELECTRON') then
              call read_in_real(PhotoElectronHeatingEfficiency, iError)
              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #PHOTOELECTRON:'
                 write(*,*) ''
                 write(*,*) '#PHOTOELECTRON'
                 write(*,*) "PhotoElectronHeatingEfficiency   (real)"
                 IsDone = .true.
              end if
           end if

        case ("#NEUTRALHEATING")
           call read_in_real(NeutralHeatingEfficiency, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #NEUTRALHEATING:'
              write(*,*) ''
              write(*,*) '#NEUTRALHEATING'
              write(*,*) "NeutralHeatingEfficiency   (real)"
              IsDone = .true.
           end if

        case ("#THERMO")
           call read_in_logical(UseSolarHeating, iError)
           call read_in_logical(UseJouleHeating, iError)
           call read_in_logical(UseAuroralHeating, iError)
           call read_in_logical(UseNOCooling, iError)
           call read_in_logical(UseOCooling, iError)
           call read_in_logical(UseConduction, iError)
           call read_in_logical(UseTurbulentCond, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #THERMO:'
              write(*,*) ''
              write(*,*) '#THERMO'
              write(*,*) "UseSolarHeating   (logical)"
              write(*,*) "UseJouleHeating   (logical)"
              write(*,*) "UseAuroralHeating (logical)"
              write(*,*) "UseNOCooling      (logical)"
              write(*,*) "UseOCooling       (logical)"
              write(*,*) "UseConduction     (logical)"
              write(*,*) "UseTurbulentCond  (logical)"
              IsDone = .true.
           endif

        case ("#THERMALDIFFUSION")
           call read_in_real(KappaTemp0, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #THERMALDIFFUSION:'
              write(*,*) ''
              write(*,*) '#THERMALDIFFUSION'
              write(*,*) "KappaTemp0    (thermal conductivity, real)"
           endif

        case ("#THERMALCONDUCTION")
           call read_in_real(ThermalConduction_AO2, iError)
           call read_in_real(ThermalConduction_AO,  iError)
           call read_in_real(ThermalConduction_s,   iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #THERMALCONDUCTION:'
              write(*,*) ''
              write(*,*) '#THERMALCONDUCTION'
              write(*,*) "ThermalConduction_AO2 (Conduction A(O2): 3.6e-4, real)"
              write(*,*) "ThermalConduction_AO  (Conduction A(O): 5.6e-4, real)"
              write(*,*) "ThermalConduction_s   (Conduction s: 0.75, real)"
           endif

        case ("#VERTICALSOURCES")
           call read_in_real(MaximumVerticalVelocity, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #VERTICALSOURCES:'
              write(*,*) ''
              write(*,*) '#VERTICALSOURCES'
              write(*,*) "MaximumVerticalVelocity      (real)"
           endif

        case ("#AUSMSOLVER")
           call read_in_logical(UseAUSMSolver, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #VERTICALSOURCES:'
              write(*,*) ''
              write(*,*) '#AUSMSOLVER'
              write(*,*) "Use AUSM Solver      (logical)"
           endif

        case ("#DIFFUSION")
           call read_in_logical(UseDiffusion, iError)
           if (UseDiffusion .and. iError == 0) then
              call read_in_real(EddyDiffusionCoef,iError)
              call read_in_real(EddyDiffusionPressure0,iError)
              call read_in_real(EddyDiffusionPressure1,iError)
           endif

           if (EddyDiffusionPressure0 < EddyDiffusionPressure1) then
              iError = 1
           endif
             
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DIFFUSION:'
              write(*,*) ''
              write(*,*) "If you use eddy diffusion, you must specify two pressure"
              write(*,*) "levels - under the first, the eddy diffusion is constant."
              write(*,*) "Between the first and the second, there is a linear drop-off."
              write(*,*) "Therefore The first pressure must be larger than the second!"
              write(*,*) '#DIFFUSION'
              write(*,*) "UseDiffusion (logical)"
              write(*,*) "EddyDiffusionCoef (real)"
              write(*,*) "EddyDiffusionPressure0 (real)"
              write(*,*) "EddyDiffusionPressure1 (real)"
           endif

        case ("#FORCING")
           call read_in_logical(UsePressureGradient, iError)
           call read_in_logical(UseIonDrag, iError)
           call read_in_logical(UseNeutralFriction, iError)
           call read_in_logical(UseViscosity, iError)
           call read_in_logical(UseCoriolis, iError)
           call read_in_logical(UseGravity, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #THERMO:'
              write(*,*) ''
              write(*,*) '#FORCING'
              write(*,*) "UsePressureGradient (logical)"
              write(*,*) "UseIonDrag          (logical)"
              write(*,*) "UseNeutralFriction  (logical)"
              write(*,*) "UseViscosity        (logical)"
              write(*,*) "UseCoriolis         (logical)"
              write(*,*) "UseGravity          (logical)"
              IsDone = .true.
           endif

        case ("#MODIFIYPLANET")
           call read_in_real(RotationPeriodInput, iError)
           call read_in_real(DaysPerYearInput, iError)
           call read_in_real(PlanetTiltInput, iError)
           OmegaBodyInput = 2.0*pi/RotationPeriodInput
           HoursPerDayInput = RotationPeriodInput/3600.0
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #MODIFYPLANET:'
              write(*,*) ''
              write(*,*) '#MODIFYPLANET'
              write(*,*) "RotationPeriodInput        (real)"
              write(*,*) "DaysPerYearInput           (real)"
              write(*,*) "DaysPerYearInput           (real)"
           endif
        
        case ("#USEIMPLICITIONMOMENTUM")
           call read_in_logical(UseImplicitFieldAlignedMomentum, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #USEIMPLICITIONMOMENTUM:'
              write(*,*) ''
              write(*,*) '#USEIMPLICITIONMOMENTUM'
              write(*,*) "UseImplicitFieldAlignedMomentum      (logical)"
           endif

        case ("#USEIMPROVEDIONADVECTION")
           call read_in_logical(UseImprovedIonAdvection, iError)
           call read_in_logical(UseNighttimeIonBCs, iError)
           call read_in_real(MinTEC, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #USEIMPROVEDIONADVECTION:'
              write(*,*) ''
              write(*,*) '#USEIMPROVEDIONADVECTION'
              write(*,*) "UseImprovedIonAdvection      (logical)"
              write(*,*) "UseNighttimeIonBCs           (logical)"
           endif

        case ("#USETESTVISCOSITY")
           call read_in_logical(UseTestViscosity, iError)
           call read_in_real(TestViscosityFactor, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #USETESTVISCOSITY:'
              write(*,*) ''
              write(*,*) '#USETESTVISCOSITY'
              write(*,*) "UseTestViscosity         (logical)"
              write(*,*) "TestViscosityFactor      (real)"
           endif

        case ("#DYNAMO")
           call read_in_logical(UseDynamo, iError)
           if (UseDynamo) then
              call read_in_real(DynamoHighLatBoundary, iError)
              call read_in_int(nItersMax, iError)
              call read_in_real(MaxResidual, iError)
              call read_in_logical(IncludeCowling, iError)
              call read_in_real(DynamoLonAverage, iError)
           endif
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DYNAMO:'
              write(*,*) ''
              write(*,*) '#DYNAMO'
              write(*,*) "UseDynamo              (logical)"
              write(*,*) "DynamoHighLatBoundary  (real)"
              write(*,*) "nItersMax              (integer)"
              write(*,*) "MaxResidual            (V,real)"
              write(*,*) "IncludeCowling         (logical)"
              write(*,*) "DynamoLonAverage       (real)"
           endif

        case ("#IONFORCING")
           call read_in_logical(UseExB, iError)
           call read_in_logical(UseIonPressureGradient, iError)
           call read_in_logical(UseIonGravity, iError)
           call read_in_logical(UseNeutralDrag, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #IONFORCING:'
              write(*,*) ''
              write(*,*) '#IONFORCING'
              write(*,*) "UseExB                 (logical)"
              write(*,*) "UseIonPressureGradient (logical)"
              write(*,*) "UseIonGravity          (logical)"
              write(*,*) "UseNeutralDrag         (logical)"
              IsDone = .true.
           endif

        case ("#CHEMISTRY")

           call read_in_logical(UseIonChemistry, iError)
           call read_in_logical(UseIonAdvection, iError)
           call read_in_logical(UseNeutralChemistry, iError)

!           call read_in_string(sNeutralChemistry, iError)
!           call read_in_string(sIonChemistry, iError)
!
!           iInputIonChemType = -1
!           iInputNeutralChemType = -1
!
!           do i = 1, nChemTypes_
!              if (sNeutralChemistry == sChemType(i)) iInputNeutralChemType = i
!              if (sIonChemistry == sChemType(i)) iInputIonChemType = i
!           enddo
!
!           if (iInputNeutralChemType < 1) then
!              write(*,*) "Error in #CHEMISTRY for Neutrals"
!              write(*,*) "Input type : ", sNeutralChemistry
!              write(*,*) "Acceptable types :"
!              do i = 1, nChemTypes_
!                 write(*,*) sChemType(i)
!              enddo
!              IsDone = .true.
!           endif
!
!           if (iInputIonChemType < 1) then
!              write(*,*) "Error in #CHEMISTRY for Ions"
!              write(*,*) "Input type : ", sIonChemistry
!              write(*,*) "Acceptable types :"
!              do i = 1, nChemTypes_
!                 write(*,*) sChemType(i)
!              enddo
!              IsDone = .true.
!           endif

        case ("#FIXTILT")

           call read_in_logical(IsFixedTilt, iError)

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #FIXTILT:'
              write(*,*) ''
              write(*,*) '#FIXTILT'
              write(*,*) 'IsFixedTilt   (logical)'
              IsDone = .true.
           endif

        case ("#DIPOLE")

           call read_in_real(MagneticPoleRotation, iError)
           call read_in_real(MagneticPoleTilt    , iError)
           call read_in_real(xDipoleCenter       , iError)
           call read_in_real(yDipoleCenter       , iError)
           call read_in_real(zDipoleCenter       , iError)

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DIPOLE:'
              write(*,*) ''
              write(*,*) '#DIPOLE'
              write(*,*) 'MagneticPoleRotation   (real)'
              write(*,*) 'MagneticPoleTilt       (real)'
              write(*,*) 'xDipoleCenter          (real)'
              write(*,*) 'yDipoleCenter          (real)'
              write(*,*) 'zDipoleCenter          (real)'
              IsDone = .true.

           else

              MagneticPoleRotation = MagneticPoleRotation * pi / 180.0
              MagneticPoleTilt     = MagneticPoleTilt     * pi / 180.0
              xDipoleCenter = xDipoleCenter * 1000.0
              yDipoleCenter = yDipoleCenter * 1000.0
              zDipoleCenter = zDipoleCenter * 1000.0

           endif

        case ("#APEX")

           if (IsFramework .and. UseApex) then
              if (iDebugLevel >= 0) then
                 write(*,*) "---------------------------------------"
                 write(*,*) "-  While using Framework, you can not -"
                 write(*,*) "-   use APEX coordinates, sorry.      -"
                 write(*,*) "-          Ignoring                   -"
                 write(*,*) "---------------------------------------"
              endif
              UseApex = .false.
           else

              call read_in_logical(UseApex, iError)
              
              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #APEX:'
                 write(*,*) ''
                 write(*,*) '#APEX'
                 write(*,*) 'UseApex (logical)'
                 write(*,*) '        Sets whether to use a realistic magnetic'
                 write(*,*) '        field (T) or a dipole (F)'
                 IsDone = .true.
              endif

           endif


        case ("#ALTITUDE")
           call read_in_real(AltMin, iError)
           call read_in_real(AltMax, iError)
           call read_in_logical(UseStretchedAltitude, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #ALTITUDE'
              write(*,*) 'For Earth, the AltMin is the only variable used here.'
              write(*,*) 'The altitudes are set to 0.3 times the scale height'
              write(*,*) 'reported by MSIS, at the equator for the specified'
              write(*,*) 'F107 and F107a values.'
              write(*,*) '#ALTITUDE'
              write(*,*) 'AltMin                (real, km)'
              write(*,*) 'AltMax                (real, km)'
              write(*,*) 'UseStretchedAltitude  (logical)'
           else
              AltMin = AltMin * 1000.0
              AltMax = AltMax * 1000.0
           endif

        case ("#GRID")
           if (nLats == 1 .and. nLons == 1) Is1D = .true.
           call read_in_int(nBlocksLon, iError)
           call read_in_int(nBlocksLat, iError)
           call read_in_real(LatStart, iError)
           call read_in_real(LatEnd, iError)
           call read_in_real(LonStart, iError)
           if (iError == 0) then 
              call read_in_real(LonEnd, iError)
              if (iError /= 0) then
                 write(*,*) "You can now model part of the globe in"
                 write(*,*) "longitude.  Include a LonEnd after the LonStart"
                 write(*,*) "variable in #GRID."
                 LonEnd = LonStart
                 iError = 0
              endif
           endif

           if (nLats > 1 .and. LatEnd-LatStart < 1) iError=1

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #GRID:'
              write(*,*) ''
              write(*,*) 'If LatStart and LatEnd are set to < -90 and'
              write(*,*) '> 90, respectively, then GITM does a whole'
              write(*,*) 'sphere.  If not, it models between the two.'
              write(*,*) 'If you want to do 1-D, set nLons=1, nLats=1 in'
              write(*,*) 'ModSizeGitm.f90, then recompile, then set LatStart'
              write(*,*) 'and LonStart to the point on the Globe you want'
              write(*,*) 'to model.'
              write(*,*) '#GRID'
              write(*,*) 'nBlocksLon   (integer)'
              write(*,*) 'nBlocksLat   (integer)'
              write(*,*) 'LatStart     (real)'
              write(*,*) 'LatEnd       (real)'
              write(*,*) 'LonStart     (real)'
              write(*,*) 'LonEnd       (real)'
              IsDone = .true.
           else

              if (LatStart <= -90.0 .and. LatEnd >= 90.0) then
                 IsFullSphere = .true.
              else
                 IsFullSphere = .false.
                 LatStart = LatStart * pi / 180.0
                 LatEnd   = LatEnd * pi / 180.0
              endif

              LonStart = LonStart * pi / 180.0
              LonEnd   = LonEnd * pi / 180.0

           endif

        case ("#NEWSTRETCH")
           NewStretchedGrid = .true.
           call read_in_real(ConcentrationLatitude, iError)
           call read_in_real(StretchWidth, iError)
           call read_in_real(StretchingPercentage, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #NEWSTRETCH:'
              write(*,*) ''
              write(*,*) '#NEWSTRETCH'
              write(*,*) 'Poleward Edge of Stretch Region (real, degrees)'
              write(*,*) 'StretchWidth  (real, 1.0-20.0, deg)'
              write(*,*) 'StretchingPercentage  (real, 0-1)'
              write(*,*) 'Example (auroral zone):'
              write(*,*) '#NEWSTRETCH'
              write(*,*) '65.0 ! location of minimum grid spacing'
              write(*,*) '5.0	 ! Width of stretched region'
              write(*,*) '0.6	 ! Amount of stretch 0 (none) to 1 (lots)'
              IsDone = .true.
           endif

        case ("#STRETCH")
           call read_in_real(ConcentrationLatitude, iError)
           call read_in_real(StretchingPercentage, iError)
           call read_in_real(StretchingFactor, iError)

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #STRETCH:'
              write(*,*) ''
              write(*,*) '#STRETCH'
              write(*,*) 'ConcentrationLatitude (real, degrees)'
              write(*,*) 'StretchingPercentage  (real, 0-1)'
              write(*,*) 'StretchingFactor      (real)'
              write(*,*) 'Example (no stretching):'
              write(*,*) '#STRETCH'
              write(*,*) '65.0 ! location of minimum grid spacing'
              write(*,*) '0.0	 ! Amount of stretch 0 (none) to 1 (lots)'
              write(*,*) '1.0  '// &
                   '! More control of stretch ( > 1 stretch less '//&
                   '< 1 stretch more)'
              IsDone = .true.
           endif

        case("#TOPOGRAPHY")
           call read_in_logical(UseTopography, iError)
            call read_in_real(AltMinUniform,iError)
           if(iError /= 0) then
              write(*,*) 'Incorrect format for #TOPOGRAPHY:'
              write(*,*) '#TOPOGRAPHY'
              write(*,*) 'UseTopography (logical)'
              write(*,*) 'AltMinUniform (real)'
              IsDone = .true.
           endif
           AltMinUniform=AltMinUniform*1000.0

        case ("#RESTART")
           call read_in_logical(DoRestart, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #RESTART:'
              write(*,*) ''
              write(*,*) '#RESTART'
              write(*,*) 'DoRestart (logical)'
              IsDone = .true.
           endif

        case ("#PLOTTIMECHANGE")
           call read_in_time(PlotTimeChangeStart, iError)
           call read_in_time(PlotTimeChangeEnd,   iError)
           do iOutputTypes=1,nOutputTypes
              call read_in_real(PlotTimeChangeDt(iOutputTypes), iError)
           enddo
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #PLOTTIMECHANGE:'
              write(*,*) 'This allows you to change the output cadence of the '
              write(*,*) 'files for a limited time.  If you have an event then '
              write(*,*) 'you can output much more often during that event.'
              write(*,*) '#PLOTTIMECHANGE'
              write(*,*) 'yyyy mm dd hh mm ss ms (start)'
              write(*,*) 'yyyy mm dd hh mm ss ms (end)'
              do iOutputTypes=1,nOutputTypes
                 write(*,*) 'dt (real)'
              enddo
              IsDone = .true.
           endif

        case ("#APPENDFILES")
           call read_in_logical(DoAppendFiles, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #APPENDFILES:'
              write(*,*) 'For satellite files, you can have one single file'
              write(*,*) 'per satellite, instead of one for every output.'
              write(*,*) 'This makes GITM output significantly less files.'
              write(*,*) 'It only works for satellite files now.'
              write(*,*) '#APPENDFILES'
              write(*,*) 'DoAppendFiles    (logical)'
           endif

        case ("#CCMCFILENAME")
           call read_in_logical(UseCCMCFileName, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #CCMCFILENAME:'
              write(*,*) 'Typicaly file is named (e.g.) 1DALL_yymmdd_hhmmss.bin'
              write(*,*) 'With this it will be named'
              write(*,*) '1DALL_GITM_yyyy-mm-ddThh-mm-ss.bin'
              write(*,*) '#CCMCFILENAME'
              write(*,*) 'UseCCMCFileName    (logical)'
           endif

        case ("#SAVEPLOTS", "#SAVEPLOT")
           call read_in_real(DtRestart, iError)
           call read_in_int(nOutputTypes, iError)
           if (nOutputTypes > nMaxOutputTypes) then
              write(*,*) "Sorry, nOutputTypes is Limited to ",nMaxOutputTypes
              iError = 1
           else
              do iOutputTypes=1,nOutputTypes
                 call read_in_string(OutputType(iOutputTypes), iError)
                 call read_in_real(DtPlot(iOutputTypes), iError)
              enddo
              iError = bad_outputtype()
              if (iError /= 0) then
                 write(*,*) "Error in #SAVEPLOTS"
                 write(*,*) "Bad output type : ",OutputType(iError)
              endif
              DtPlotSave = DtPlot
           endif
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #SAVEPLOT'
              write(*,*) 'The DtRestart variable sets the time in between '
              write(*,*) 'writing full restart files to the UA/restartOUT '
              write(*,*) 'directory.\\'
              write(*,*) 'This sets the output files.  The most common type'
              write(*,*) 'is 3DALL, which outputs all primary state variables.'
              write(*,*) 'Types include : 3DALL, 3DNEU, 3DION, 3DTHM, 3DCHM, '
              write(*,*) '3DUSR, 3DGLO, 2DGEL, 2DMEL, 2DUSR, 1DALL, 1DGLO, '
              write(*,*) '1DTHM, 1DNEW, 1DCHM, 1DCMS, 1DUSR. DtPlot sets the'
              write(*,*) 'frequency of output'
              write(*,*) '#SAVEPLOT'
              write(*,*) 'DtRestart (real, seconds)'
              write(*,*) 'nOutputTypes  (integer)'
              write(*,*) 'Outputtype (string, 3D, 2D, ION, NEUTRAL, ...)'
              write(*,*) 'DtPlot    (real, seconds)'
              IsDone = .true.
           endif

        case ("#SAVEHIMEPLOT")
           call read_in_real(HIMEPlotLonStart, iError)
           call read_in_real(HIMEPlotLonEnd, iError)
           call read_in_real(HIMEPlotLatStart, iError)
           call read_in_real(HIMEPlotLatEnd, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #SAVEHIMEPLOT'
              write(*,*) ''
              write(*,*) '#SAVEHIMEPLOT'
              write(*,*) 'HIMEPlotLonStart (real)'
              write(*,*) 'HIMEPlotLonEnd (real)'
              write(*,*) 'HIMEPlotLatStart (real)'
              write(*,*) 'HIMEPlotLatEnd (real)'
              IsDone = .true.
           endif

        case ("#SATELLITES")
           call read_in_int(nSats, iError)
           if (nSats > nMaxSats) then
              iError = 1
              write(*,*) "Too many satellites requested!!"
           else if(RCMRFlag .and. nRCMRSat > nSats) then
              iError = 1
              write (*,*) "Assimilation satellite index out of range"
           else
              do iSat=1,nSats
                 call read_in_string(cSatFileName(iSat), iError)
                 call read_in_string(SatOutputType(iSat), iError)       !!! Xing 
                 call read_in_real(SatDtPlot(iSat), iError)
                 iSatCurrentIndex(iSat) = 0
              enddo
           endif
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #SATELLITES'
              write(*,*) ''
              write(*,*) '#SATELLITES'
              write(*,*) 'nSats     (integer - max = ',nMaxSats,')'
              write(*,*) 'SatFile1  (string)'
              write(*,*) 'SatOutputtype1  (string, 0DUSR or 1DUSR)'
              write(*,*) 'DtPlot1   (real, seconds)'
              write(*,*) 'etc...'
              IsDone = .true.
           else
              call read_satellites(iError)
              if (iError /= 0) IsDone = .true.
           endif

        case("#RCMR")
           ! Read in the data assimilation type
           call read_in_string(RCMRInType, iError)
           ! Read in the output driver type
           call read_in_string(RCMROutType, iError)

           if(iError == 0) then
              if(RCMROutType == 'F107') then
                 RCMRFlag = .true.
                 call read_in_real(f107_est, iError)

                 if(iError /= 0) then
                    write (*,*) 'Unrecognizable initial F10.7 estimate'
                    write (*,*) 'running with 100.0, which is only appropriate '
                    write (*,*) 'if the actual F10.7 is not close to 100'
                    f107_est = 100
                    iError   = 0
                 end if

                 f107a_est = f107_est
                 call IO_set_f107_single(f107_est)
                 call IO_set_f107a_single(f107a_est)

              else if(RCMROutType == 'PHOTOELECTRON') then
                 RCMRFlag = .true.
                 call read_in_real(PhotoElectronHeatingEfficiency_est, iError)

                 if(iError /= 0) then
                    write (*,*) 'Unrecognizable initial photoelectron heating '
                    write (*,*) 'efficiency.  Running with 0.1, which is only '
                    write (*,*) 'appropraite if the actual efficiency is not '
                    write (*,*) 'close to 0.1'
                    PhotoElectronHeatingEfficiency_est = 0.1
                    iError = 0
                 end if

                 if(RCMROutType == 'PHOTOELECTRON') then
                    PhotoElectronHeatingEfficiency = PhotoElectronHeatingEfficiency_est
                 end if
              
                 !!ANKIT: TEC read
              else if(RCMROutType == 'EDC') then
                 RCMRFlag = .true.
                 call read_in_real(EDC_est_tmp, iError)
                 EDC_est(1,1) = EDC_est_tmp
                 ! call read_in_real(EDC_est, iError)
                 
                 if(iError /= 0) then
                    write (*,*) 'Unrecognizable initial eddy diffusion '
                    write (*,*) 'coefficient.  Running with 1750, which is only '
                    write (*,*) 'appropraite if the actual coefficient is not '
                    write (*,*) 'close to 1750'
                    EDC_est = 1750
                    iError = 0
                 end if
                 
                 if(RCMROutType == 'EDC') then
                    EddyDiffusionCoef = EDC_est(1,1)
                 end if
              end if
              ! Quality check of the satellite indexes will occur when
              ! reading the satellite block.
              call read_in_int(nRCMRSat, iError)
              do iSat=1,nRCMRSat
                 call read_in_int(RCMRSat(iSat), iError)
              end do
              call read_in_string(RCMRrun, iError)    !Ankit23May16: Truth or ID
              call read_in_int(Nc, iError)            !Ankit23May16: Controller Order.
              call read_in_int(RegZ, iError)          !Ankit23May16: Regressor Variable. 1=z, 0=y.
              call read_in_real(lambda1, iError)      !Ankit23May16: Forgetting factor. 
              call read_in_real(W_Rtheta, iError)     !Ankit23May16: Rtheta. Or inv(P0).
              call read_in_real(Dts, iError)          !Ankit23May16: Time step used with RCMR
              call read_in_real(Measure_Dts, iError)  !Ankit23May16: RCMR executes every Measure_Dts*Dts*steps timesteps.
              call read_in_int(C_on, iError)          !Ankit23May16: RCMR switches on after C_on steps.
              call read_in_int(uFiltLength, iError)   !Ankit23May16: Length of averaging filter.
           end if
           
           if(iError /= 0) then
              write (*,*) 'Incorrect format for #RCMR'
              write (*,*) ''
              write (*,*) '#RCMR'
              write (*,*) 'Input data type to assimilate (RHO/VTEC)'
              write (*,*) 'Output variable to drive (F107/PHOTOELECTRON)'
              write (*,*) 'Initial output estimate'
              write (*,*) 'Number satellites to assimilate'
              write (*,*) '1st satellite index to assimilate'
              write (*,*) '2nd satellite index to assimilate'
              write (*,*) 'etc...'
              write (*,*) ''
              write (*,*) 'If driving F10.7, use #F10.7 input block to set'
              write (*,*) 'the MSIS F10.7 and F10.7a to realistic values'
              IsDone = .true.
           else
              if(RCMROutType == 'EDC') then
                 call read_tec_locations              
              end if
              call alloc_rcmr
              call init_markov_matrix
              call init_rcmr
           end if

        case ("#DART") !alexey
           call read_in_int(useDART, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DART:'
              write(*,*) ''
              write(*,*) '#DART'
              write(*,*) "useDART (integer, {default 0=no}, 1=master ensemble member, 2=slave ens.)"
           endif

        case ("#ELECTRODYNAMICS")
           call read_in_real(dTPotential, iError)
           call read_in_real(dTAurora, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #ELECTRODYNAMICS'
              write(*,*) 'Sets the time for updating the high-latitude'
              write(*,*) '(and low-latitude) electrodynamic drivers, such as'
              write(*,*) 'the potential and the aurora.'
              write(*,*) '#ELECTRODYNAMICS'
              write(*,*) 'DtPotential (real, seconds)'
              write(*,*) 'DtAurora    (real, seconds)'
              IsDone = .true.
           endif

        case ("#INPUTTIMEDELAY")
           call read_in_real(TimeDelayHighLat, iError)
           call read_in_real(TimeDelayEUV, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #INPUTTIMEDELAY'
              write(*,*) 'Sets the time delay for the high latitude drivers'
              write(*,*) 'and solar EUV inputs.'
              write(*,*) '#INPUTTIMEDELAY'
              write(*,*) 'TimeDelayHighLat (real, seconds)'
              write(*,*) 'TimeDelayEUV     (real, seconds)'
              IsDone = .true.
           endif

        case ("#LTERadiation")
           call read_in_real(DtLTERadiation, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #LTERadiation:'
              write(*,*) ''
              write(*,*) '#LTERadiation'
              write(*,*) 'DtLTERadiation (real)'
              IsDone = .true.
           endif

        case ("#IONPRECIPITATION")
           call read_in_logical(UseIonPrecipitation, iError)
           call read_in_string(IonIonizationFilename, iError)
           call read_in_string(IonHeatingRateFilename, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #IONPRECIPITATION'
              write(*,*) 'This is really highly specific.  You dont want this.'
              write(*,*) '#IONPRECIPITATION'
              write(*,*) 'UseIonPrecipitation     (logical)'
              write(*,*) 'IonIonizationFilename   (string)'
              write(*,*) 'IonHeatingRateFilename  (string)'
              IsDone = .true.
           endif

        case ("#LOGFILE")
           call read_in_real(dTLogFile, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #LOGFILE'
              write(*,*) 'You really want a log file.  They are very important.'
              write(*,*) 'It is output in UA/data/log*.dat.'
              write(*,*) 'You can output the log file at whatever frequency'
              write(*,*) 'you would like, but if you set dt to some very small'
              write(*,*) 'value, you will get an output every iteration, which'
              write(*,*) 'is probably a good thing.'
              write(*,*) '#LOGFILE'
              write(*,*) 'DtLogFile   (real, seconds)'
              IsDone = .true.
           endif

        case ("#EUVMODEL")
           call read_in_logical(UseEUVAC, iError)
           call read_in_logical(UseTobiska, iError)
           call read_in_logical(UseAboveHigh, iError)
           call read_in_logical(UseBelowLow, iError)

           if ( .not. UseEUVAC     .and. &
                .not. UseTobiska   .and. &
                .not. UseAboveHigh .and. &
                .not. UseBelowLow) UseRidleyEUV = .true.

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #EUVMODEL'
              write(*,*) 'If EUVAC   = true, use EUVAC Model'
              write(*,*) 'If Tobiska = true, use Tobiska91 Model'
              write(*,*) 'If both are true, average them together'
              write(*,*) 'If UseAboveHigh, then extend the spectrum to the longer wavelengths'
              write(*,*) 'If UseBelowLow, then extend the spectrum to the shorter wavelengths'
              write(*,*) '#EUV_DATA'
              write(*,*) 'UseEUVData            (logical)'
              write(*,*) 'cEUVFile              (string)'
           endif

        case ("#EUV_DATA")
           call read_in_logical(UseEUVData, iError)
           call read_in_string(cEUVFile, iError)

           if (UseEUVData) call Set_Euv(iError, CurrentTime, EndTime)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #EUV_DATA'
              write(*,*) 'This is for a FISM or some other solar spectrum file.'
              write(*,*) '#EUV_DATA'
              write(*,*) 'UseEUVData            (logical)'
              write(*,*) 'cEUVFile              (string)'
           endif

        case ("#ECLIPSE")
           IncludeEclipse = .true.
           call read_in_time(EclipseStartTime, iError)
           call read_in_time(EclipseEndTime, iError)
           call read_in_real(EclipseStartY, iError)
           call read_in_real(EclipseStartZ, iError)
           call read_in_real(EclipseEndY, iError)
           call read_in_real(EclipseEndZ, iError)
           call read_in_real(EclipsePeak, iError)
           call read_in_real(EclipseMaxDistance, iError)
           call read_in_real(EclipseExpAmp, iError)
           call read_in_real(EclipseExpWidth, iError)

           EclipseStartY = EclipseStartY * 1000.0
           EclipseStartZ = EclipseStartZ * 1000.0
           EclipseEndY = EclipseEndY * 1000.0
           EclipseEndZ = EclipseEndZ * 1000.0
           EclipseMaxDistance = EclipseMaxDistance * 1000.0
           EclipseExpWidth = EclipseExpWidth * 1000.0
           
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #ECLIPSE'
              write(*,*) '#ECLIPSE'
              write(*,*) 'start y m d h m s            (int x 6)'
              write(*,*) 'end   y m d h m s            (int x 6)'
              write(*,*) 'start-y-gse                  (real)'
              write(*,*) 'start-z-gse                  (real)'
              write(*,*) 'end-y-gse                    (real)'
              write(*,*) 'end-z-gse                    (real)'
              write(*,*) 'peak                         (real)'
              write(*,*) 'Max Distance                 (real)'
              write(*,*) 'ExpAmp                       (real)'
              write(*,*) 'ExpWidth                     (real)'
           endif

        case ("#GLOW")
           call read_in_logical(UseGlow, iError) 
           call read_in_real(dTGlow, iError)

        case ("#NGDC_INDICES")
           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           cTempLines(2) = cTempLine
           cTempLines(3) = " "
           cTempLines(4) = "#END"

           call IO_set_inputs(cTempLines)
           call read_NGDC_Indices_new(iError,CurrentTime,EndTime)

           if (iError /= 0) then 
              write(*,*) "read indices was NOT successful (NOAA file)"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif
        
        case ("#OMNIWEB_AP_INDICES")
           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           write(*,*) cTempLine
           cTempLines(2) = cTempLine
           cTempLines(3) = " "
           cTempLines(4) = "#END"

           call IO_set_inputs(cTempLines)
           call read_OMNIWEB_Ap_Indices_new(iError,CurrentTime,EndTime)

           if (iError /= 0) then
              write(*,*) "read indices was NOT successful (OMNIWEB Ap file)"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif

        case ("#SME_INDICES")
           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           cTempLines(2) = cTempLine
           call read_in_string(cTempLine, iError)
           cTempLines(3) = cTempLine
           cTempLines(4) = " "
           cTempLines(5) = "#END"

           call IO_set_inputs(cTempLines)
           call read_sme(iError,CurrentTime+TimeDelayHighLat,EndTime+TimeDelayHighLat)

           if (iError /= 0) then
              write(*,*) "read indices was NOT successful (SME file)"
              IsDone = .true.
           endif

           call read_al_onset_list(iError,CurrentTime+TimeDelayHighLat,EndTime+TimeDelayHighLat)

           if (iError /= 0) then
              write(*,*) "read indices was NOT successful (onset file)"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif

        case ("#ACE_DATA")
           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           cTempLines(2) = cTempLine
           call read_in_string(cTempLine, iError)
           cTempLines(3) = cTempLine
           cTempLines(4) = " "
           cTempLines(5) = "#END"
           call IO_set_inputs(cTempLines)
           call read_ACE_data(iError)

           if (iError /= 0) then 
              write(*,*) "read indices was NOT successful (ACE SC file)"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif

        case ("#SWPC_INDICES")
           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           cTempLines(2) = cTempLine
           call read_in_string(cTempLine, iError)
           cTempLines(3) = cTempLine
           cTempLines(4) = " "
           cTempLines(5) = "#END"

           call IO_set_inputs(cTempLines)
           call read_SWPC_Indices(iError)

           if (iError /= 0) then 
              write(*,*) "read indices was NOT successful (SWPC file)"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif

        case ("#NOAAHPI_INDICES")
           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           cTempLines(2) = cTempLine
           cTempLines(3) = " "
           cTempLines(4) = "#END"

           call IO_set_inputs(cTempLines)
           call read_NOAAHPI_Indices_new(iError,StartTime,EndTime)

           if (iError /= 0) then 
              write(*,*) "read indices was NOT successful (NOAA HPI file)"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif

        case ("#END")
           IsDone = .true.

        end select

        if (iError /= 0) IsDone = .true.

    endif

    iLine = iLine + 1

    if (iLine >= nInputLines) IsDone = .true.

  enddo

  if (iDebugProc >= 0 .and. iProc /= iDebugProc) then
     iDebugLevel = -1
  endif

  if (iError /= 0) then
     call stop_gitm("Must Stop!!")
  endif

  RestartTime = CurrentTime

!  KappaTemp0 = 3.6e-4

contains

  subroutine read_in_int(variable, iError)
    integer, intent(out) :: variable
    integer :: iError
    if (iError == 0) then 
       iline = iline + 1
       read(cInputText(iline),*,iostat=iError) variable
    endif
  end subroutine read_in_int

  subroutine read_in_time(variable, iError)
    real(Real8_), intent(out) :: variable
    integer :: iError
    integer :: iTimeReadIn(7)
    if (iError == 0) then 
       iline = iline + 1
       read(cInputText(iline),*,iostat=iError) iTimeReadIn(1:6)
       iTimeReadIn(7) = 0
       call time_int_to_real(iTimeReadIn, variable)
    endif
  end subroutine read_in_time

  subroutine read_in_logical(variable, iError)
    logical, intent(out) :: variable
    integer :: iError
    if (iError == 0) then
       iline = iline + 1
       read(cInputText(iline),*,iostat=iError) variable
    endif
  end subroutine read_in_logical

  subroutine read_in_string(variable, iError)
    character (len=iCharLen_), intent(out) :: variable
    integer :: iError
    if (iError == 0) then 
       iline = iline + 1
       variable = cInputText(iline)
       ! Remove anything after a space or TAB
       i=index(variable,' '); if(i>0)variable(i:len(cLine))=' '
       i=index(variable,char(9)); if(i>0)variable(i:len(cLine))=' '
    endif
  end subroutine read_in_string

  subroutine read_in_real(variable, iError)
    real, intent(out) :: variable
    integer :: iError
    if (iError == 0) then
       iline = iline + 1
       read(cInputText(iline),*, iostat=iError) variable
    endif
  end subroutine read_in_real

end subroutine set_inputs


subroutine fix_vernal_time

  use ModTime
  use ModPlanet
  use ModInputs, only: DaysPerYearInput, RotationPeriodInput

  implicit none

  real*8 :: DTime
  integer :: n
  integer, external :: jday

  DTime = CurrentTime - VernalTime

  if (DTime < 0.0) then
     n = 0
     do while (CurrentTime < VernalTime)
        VernalTime = VernalTime-int(DaysPerYearInput)*RotationPeriodInput
        n = n + 1
        if (floor(n*(DaysPerYearInput - int(DaysPerYearInput))) > 0) then
           VernalTime = VernalTime - &
                floor(n*(DaysPerYearInput - int(DaysPerYearInput))) * &
                RotationPeriodInput
           n = 0
        endif
     enddo
     DTime = CurrentTime - VernalTime
  else
     n = 0
     do while (DTime > DaysPerYearInput*RotationPeriodInput)
        VernalTime = VernalTime+int(DaysPerYearInput)*RotationPeriodInput
        n = n + 1
        if (floor(n*(DaysPerYearInput - int(DaysPerYearInput))) > 0) then
           VernalTime = VernalTime + &
                floor(n*(DaysPerYearInput - int(DaysPerYearInput))) * &
                RotationPeriodInput
           n = 0
        endif
        DTime = CurrentTime - VernalTime
     enddo
  endif
  iDay  = DTime / RotationPeriodInput
  uTime = (DTime / RotationPeriodInput - iDay) * RotationPeriodInput
  iJulianDay = jday( &
       iTimeArray(1), iTimeArray(2), iTimeArray(3)) 

end subroutine fix_vernal_time

