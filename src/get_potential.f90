! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine init_get_potential

  use ModGITM
  use ModTime
  use ModIndicesInterfaces
  use ModInputs
  use ModUserGITM
  use ModNewell
  use ModOvationSME
  use ModAeAuroralModel
  use ModEIE_Interface, only: EIEr3_HaveLats, EIEr3_HaveMLTs

  implicit none

  character (len=iCharLen_), dimension(100) :: Lines
  character (len=iCharLen_) :: TimeLine
  real    :: bz

  logical :: IsFirstTime = .true.
  integer :: iError

  integer :: nAmieLats, nAmieMlts, nAmieBlocks

  iError = 0

  if (.not.IsFirstTime .or. IsFramework) return
  call report("init_get_potential",2)

  IsFirstTime = .false.

  if (UseNewellAurora) then
     call init_newell
     UseIMF = .true.
  endif

  if (UseOvationSME) then
     call read_ovationsm_files
  endif

  if (UseAeModel) then
     call read_ae_model_files(iError)
  endif

  call report("AMIE vs Weimer",4)

  !!! Xing Meng Nov 2018 added UseRegionalAMIE to set up a local region
  !!! with the potential from AMIE files and Weimer potential elsewhere
  if (index(cAMIEFileNorth,"none") > 0 .or. UseRegionalAMIE) then

     Lines(1) = "#BACKGROUND"
     Lines(2) = "EIE/"

     UseHPI = .true.
     if (UseNewellAurora .or. UseOvationSME .or. UseAeModel) UseHPI = .false.

     call get_IMF_Bz(CurrentTime+TimeDelayHighLat, bz, iError)
     
     call IO_SetIMFBz(bz)
     if (iError /= 0) then
        write(*,*) "Can not find IMF Bz."
        !write(*,*) "Setting potential to Millstone HPI."
        !Lines(3) = "millstone_hpi"    ! Change to "zero" if you want
        call stop_gitm("must stop! Check your IMF file!!!")
     else
        Lines(3) = PotentialModel    ! Change to "zero" if you want
        UseIMF = .true.
     endif
     Lines(4) = AuroralModel
     Lines(5) = "idontknow"
     Lines(6) = ""

  else

     if (index(cAMIEFileNorth,"mhd") > 0) then

        Lines(1) = "#MHDFILE"
        Lines(2) = cAMIEFileSouth
        write(TimeLine,'(i4)') iTimeArray(1)
        Lines(3) = TimeLine
        write(TimeLine,'(i2)') iTimeArray(2)
        Lines(4) = TimeLine
        write(TimeLine,'(i2)') iTimeArray(3)
        Lines(5) = TimeLine
        Lines(6) = ""

        UseIMF = .false.
           
     else

        Lines(1) = "#AMIEFILES"
        Lines(2) = cAMIEFileNorth
        Lines(3) = cAMIEFileSouth
        Lines(4) = ""
        Lines(5) = ""
        Lines(6) = ""
        
        UseIMF = .false. 

     endif

  endif

  Lines(7) = "#DEBUG"
  Lines(8) = "0"
  Lines(9) = "0"
  Lines(10) = ""

  if (IsFixedTilt) then
     Lines(11) = "#FIXTILT"
     Lines(12) = "T"
     Lines(13) = ""
     Lines(14) = "#END"
  else
     Lines(11) = "#END"
  endif

  call report("EIE_stuff",4)

  call EIE_set_inputs(Lines)

  call EIE_Initialize(iError)
  if (iError /= 0) then
     write(*,*) &
          "Code Error in IE_Initialize called from get_potential.f90"
     write(*,*) "Error : ",iError
     call stop_gitm("Stopping in get_potential")
  endif

  if (UseRegionalAMIE) then
     !!! read AMIE files
     call AMIE_SetFileName(cAMIEFileNorth)
     call readAMIEOutput(2, .false., iError)
     if (iError /= 0) then
        write(*,*) &
             "Code Error in readAMIEOutput called from get_potential.f90"
        write(*,*) "Error : ",iError
        call stop_gitm("Stopping in get_potential")
     endif

     if (index(cAMIEFileSouth,'mirror') > 0) then
        call AMIE_SetFileName(cAMIEFileNorth)
        call readAMIEOutput(1, .true., iError)
     else
        call AMIE_SetFileName(cAMIEFileSouth)
        call readAMIEOutput(1, .false., iError)
     endif

     call AMIE_GetnLats(nAmieLats)
     call AMIE_GetnMLTs(nAmieMlts)
     nAmieBlocks = 2

     call EIE_InitGrid(nAmieLats, nAmieMlts, nAmieBlocks, iError)
     if (iError /= 0) then
        write(*,*) &
             "Code Error in EIE_InitGrid called from get_potential.f90"
        write(*,*) "Error : ",iError
        call stop_gitm("Stopping in get_potential")
     endif

     call AMIE_GetLats(nAmieMlts,nAmieLats,nAmieBlocks, &
          EIEr3_HaveLats,iError)     

     call AMIE_GetMLTs(nAmieMlts,nAmieLats,nAmieBlocks, &
          EIEr3_HaveMLTs,iError)

  endif


end subroutine init_get_potential

!--------------------------------------------------------------------
! set_indices
!--------------------------------------------------------------------

subroutine set_indices

  use ModGITM
  use ModTime
  use ModIndicesInterfaces
  use ModInputs
  use ModUserGITM
  use ModNewell

  implicit none

  real :: temp, by, bz, vx, den
  integer :: iError
  iError = 0

  if (UseIMF) then

     call read_NOAAHPI_Indices_new( &
          iError, &
          CurrentTime+TimeDelayHighLat, &
          EndTime+TimeDelayHighLat)
     call read_MHDIMF_Indices_new( &
          iError, &
          CurrentTime+TimeDelayHighLat, &
          EndTime+TimeDelayHighLat)
     call get_IMF_Bz(CurrentTime+TimeDelayHighLat, bz, iError)
     if (bz < -50.0) bz = -50.0
     if (bz >  50.0) bz =  50.0
     call IO_SetIMFBz(bz)

     if (iError /= 0) then
        write(*,*) &
             "Code Error in get_IMF_Bz called from get_potential.f90"
        write(*,*) "Code : ",iError
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 1) write(*,*) "==> IMF Bz : ",bz

     call get_IMF_By(CurrentTime+TimeDelayHighLat, by, iError)
     if (by < -50.0) by = -50.0
     if (by >  50.0) by =  50.0
     call IO_SetIMFBy(by)

     if (iError /= 0) then
        write(*,*) &
             "Code Error in get_IMF_By called from get_potential.f90"
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 1) write(*,*) "==> IMF By : ",by

     call get_SW_V(CurrentTime+TimeDelayHighLat, vx, iError)
     if (vx < -2000.0) vx = -2000.0
     if (vx >  2000.0) vx =  2000.0
     call IO_SetSWV(vx)

     if (iError /= 0) then
        write(*,*) "Code Error in get_sw_v called from get_potential.f90"
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 1) write(*,*) "==> Solar Wind Velocity : ",vx

     call get_SW_N(CurrentTime+TimeDelayHighLat, den, iError)
     if (iError /= 0) den = 5.0
     call IO_SetSWN(den)

!!!     call get_kp(CurrentTime+TimeDelayHighLat, temp, iError)
!!!     call IO_Setkp(temp)
!!!
!!!     if (iError /= 0) then
!!!        write(*,*) "Code Error in get_kp called from get_potential.f90"
!!!        call stop_gitm("Stopping in get_potential")
!!!     endif

     if (UseNewellAurora) call calc_dfdt(by, bz, vx)

  endif

  if (UseHPI) then

     call get_HPI(CurrentTime+TimeDelayHighLat, temp, iError)
     call IO_SetHPI(temp)

     if (iError /= 0) then
        write(*,*) "Code Error in get_hpi called from get_potential.f90"
        call stop_gitm("Stopping in get_potential")
     endif

  endif

  if (index(cAMIEFileNorth,"none") <= 0 .and. &
     index(cAMIEFileNorth,"SPS") <= 0 .and. (.not. UseRegionalAMIE)) then
     ! when UseRegionalAMIE, only get AMIE values during the specified
     ! time interval.
     if (iDebugLevel > 1) &
          write(*,*) "==> Reading AMIE values for time :",CurrentTime
     call get_AMIE_values(CurrentTime+TimeDelayHighLat)
  endif

end subroutine set_indices

!--------------------------------------------------------------------
! get_potential
!--------------------------------------------------------------------

subroutine get_potential(iBlock)

  use ModGITM
  use ModTime
  use ModIndicesInterfaces
  use ModInputs
  use ModUserGITM
  use ModNewell
  use ModOvationSME, only: run_ovationsme
  use ModAeAuroralModel, only: run_ae_model
  use ModEIE_Interface, only: UAl_UseGridBasedEIE
  use ModMpi

  implicit none

  integer, intent(in) :: iBlock

  integer :: iError, iLat, iLon, iAlt, iPot, nPot, iDir, nDir=1
  logical :: IsFirstTime = .true.
  logical :: IsFirstPotential(nBlocksMax) = .true.
  logical :: IsFirstAurora(nBlocksMax) = .true.
  real    :: mP, dis, TempWeight
  real    ::  LocalSumDiffPot, MeanDiffPot

  real, dimension(-1:nLons+2, -1:nLats+2) :: TempPotential2d
  real, dimension(-1:nLons+2, -1:nLats+2, 2) :: TempPotential, AMIEPotential
  real, dimension(-1:nLons+2, -1:nLats+2) :: Grid, dynamo, SubMLats, SubMLons
  real, dimension(-1:nLons+2, -1:nLats+2) :: lats, mlts, EFlux
  real :: by, bz, CuspLat, CuspMlt

  call start_timing("get_potential")
  call report("get_potential",2)

  iError = 0

  if (index(cPlanet,"Earth") == 0) then 

     potential = 0.0
     ElectronAverageEnergy = 0.1
     ElectronEnergyFlux = 0.0001
     return

  endif

  if (floor((tSimulation-dt)/DtPotential) /= &
       floor((tsimulation)/DtPotential) .or. IsFirstPotential(iBlock)) then

     call report("Setting up IE Grid",1)

     call init_get_potential
     call UA_SetnMLTs(nLons+4)
     call UA_SetnLats(nLats+4)
     if (.not. IsFramework) then 
        call IO_SetTime(CurrentTime)
        call set_indices
     endif
     call UA_SetNorth

     call report("Getting Potential",1)

     Potential(:,:,:,iBlock) = 0.0

     do iAlt=-1,nAlts+2

        call UA_SetGrid(                    &
             MLT(-1:nLons+2,-1:nLats+2,iAlt), &
             MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), iError)

        if (iError /= 0) then
           write(*,*) "Error in routine get_potential (UA_SetGrid):"
           write(*,*) iError
           call stop_gitm("Stopping in get_potential")
        endif

        if (iDebugLevel > 1 .and. iAlt == 1) &
             write(*,*) "==> Getting IE potential"

        TempPotential = 0.0
        TempPotential2d = 0.0

        call UA_GetPotential(TempPotential2d, iError)
        TempPotential(:,:,1) = TempPotential2d

        if (iError /= 0) then
           write(*,*) "Error in get_potential (UA_GetPotential):"
           write(*,*) iError
           TempPotential = 0.0
!           call stop_gitm("Stopping in get_potential")
        endif

        nDir = 1
        if (UseRegionalAMIE .and. &
             CurrentTime >= AMIETimeStart .and. CurrentTime <= AMIETimeEnd) then

           if (iDebugLevel > 1) &
                write(*,*) "==> Reading AMIE values for time :",CurrentTime
           call get_AMIE_values(CurrentTime+TimeDelayHighLat)
           
           if (iDebugLevel > 1) write(*,*) "==> Getting AMIE Potential"

           AMIEPotential = 0.0
           UAl_UseGridBasedEIE = .true.

           call UA_SetGrid(                    &
             MLT(-1:nLons+2,-1:nLats+2,iAlt), &
             MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), iError)           

           call UA_GetPotential(AMIEPotential(:,:,1), iError)
           if (iError /= 0) then
              write(*,*) "Error in get_potential (UA_GetPotential) for AMIE:"
              write(*,*) iError
              call stop_gitm("Stopping in get_potential")
           endif

           if (UseTwoAMIEPotentials) then
              if (iDebugLevel > 1) &
                   write(*,*) "==> Reading AMIE PotentialY for time :",CurrentTime
              call get_AMIE_PotentialY(CurrentTime+TimeDelayHighLat)
!              the following setgrid seems unnecessary
!              call UA_SetGrid(                    &
!                   MLT(-1:nLons+2,-1:nLats+2,iAlt), &
!                   MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), iError)
              call UA_GetPotential(AMIEPotential(:,:,2), iError)
              if (iError /= 0) then
                 write(*,*) "Error in get_potential (UA_GetPotential) for AMIE PotentialY:"
                 write(*,*) iError
                 call stop_gitm("Stopping in get_potential")
              endif
              TempPotential(:,:,2) = TempPotential(:,:,1)
              ! two directions for two-potentials
              nDir = 2
           else
              nDir = 1
           endif

           ! Two iterations for two-potentials, one iteration for one-potential
           do iDir = 1, nDir

              ! AMIE potential inside a geographic rectangular; Weimer potential outside.
              ! Spatially weighted AMIE and Weimer potentials at
              ! the rectangular boundaries, within AMIEBoundaryWidth.
              LocalSumDiffPot = 0.0
              iPot = 0
              do iLat = -1,nLats+2
                 do iLon = -1,nLons+2
                    
                    ! find the mean of the difference between AMIE potential and 
                    ! Weimer potential within the region
                    if (MLatitude(iLon,iLat,iAlt,iBlock) >= AMIELatStart &
                         .and. MLatitude(iLon,iLat,iAlt,iBlock) <= AMIELatEnd &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) >= AMIELonStart  &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) <= AMIELonEnd) then
                       LocalSumDiffPot = LocalSumDiffPot + &
                            AMIEPotential(iLon,iLat,iDir) - TempPotential(iLon,iLat,iDir)
                       iPot = iPot + 1
                    endif
                 enddo
              enddo
              call MPI_ALLREDUCE(LocalSumDiffPot, MeanDiffPot, 1, &
                   MPI_REAL, MPI_SUM, iCommGITM, iError)
              call MPI_ALLREDUCE(iPot, nPot, 1, &
                   MPI_INTEGER, MPI_SUM, iCommGITM, iError)
              MeanDiffPot = MeanDiffPot/nPot
              if (iDebugLevel > 1 .and. iProc == 0 .and. iAlt == 10) &
                   write(*,*)  'iDir, nPot, MeanDiffPot = ', iDir, nPot, MeanDiffPot

              do iLat = -1,nLats+2
                 do iLon = -1,nLons+2

                    ! adjust AMIE potential according to the mean of potential difference
                    AMIEPotential(iLon,iLat,iDir) = AMIEPotential(iLon,iLat,iDir) - MeanDiffPot

                    ! fill in the rectangular with AMIE potential
                    if (MLatitude(iLon,iLat,iAlt,iBlock) >= AMIELatStart &
                         .and. MLatitude(iLon,iLat,iAlt,iBlock) <= AMIELatEnd &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) >= AMIELonStart  &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) <= AMIELonEnd) &
                         TempPotential(iLon,iLat,iDir) = AMIEPotential(iLon,iLat,iDir)
                     
                    ! left and right boundaries, excluding corners
                    if (MLatitude(iLon,iLat,iAlt,iBlock) > AMIELatStart &
                         .and. MLatitude(iLon,iLat,iAlt,iBlock) < AMIELatEnd &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) > AMIELonStart-AMIEBoundaryWidth &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) < AMIELonStart) then
                       TempWeight = (AMIELonStart-MLongitude(iLon,iLat,iAlt,iBlock))/AMIEBoundaryWidth
                       TempPotential(iLon,iLat,iDir) = TempWeight*TempPotential(iLon,iLat,iDir) + &
                            (1-TempWeight)*AMIEPotential(iLon,iLat,iDir)
                    endif

                    if (MLatitude(iLon,iLat,iAlt,iBlock) > AMIELatStart &
                         .and. MLatitude(iLon,iLat,iAlt,iBlock) < AMIELatEnd &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) > AMIELonEnd &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) < AMIELonEnd+AMIEBoundaryWidth) then
                       TempWeight = (MLongitude(iLon,iLat,iAlt,iBlock)-AMIELonEnd)/AMIEBoundaryWidth
                       TempPotential(iLon,iLat,iDir) = TempWeight*TempPotential(iLon,iLat,iDir) + &
                            (1-TempWeight)*AMIEPotential(iLon,iLat,iDir)
                    endif

                    ! low and high lat boundaries
                    if (MLongitude(iLon,iLat,iAlt,iBlock) > AMIELonStart-AMIEBoundaryWidth &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) < AMIELonEnd+AMIEBoundaryWidth &
                         .and. MLatitude(iLon,iLat,iAlt,iBlock) > AMIELatStart-AMIEBoundaryWidth &
                         .and. MLatitude(iLon,iLat,iAlt,iBlock) < AMIELatStart) then
                       TempWeight = (AMIELatStart-MLatitude(iLon,iLat,iAlt,iBlock))/AMIEBoundaryWidth
                       TempPotential(iLon,iLat,iDir) = TempWeight*TempPotential(iLon,iLat,iDir) + &
                            (1-TempWeight)*AMIEPotential(iLon,iLat,iDir)
                    endif
                    if (MLongitude(iLon,iLat,iAlt,iBlock) > AMIELonStart-AMIEBoundaryWidth &
                         .and. MLongitude(iLon,iLat,iAlt,iBlock) < AMIELonEnd+AMIEBoundaryWidth &
                         .and. MLatitude(iLon,iLat,iAlt,iBlock) > AMIELatEnd &
                         .and. MLatitude(iLon,iLat,iAlt,iBlock) < AMIELatEnd+AMIEBoundaryWidth) then
                       TempWeight = (MLatitude(iLon,iLat,iAlt,iBlock)-AMIELatEnd)/AMIEBoundaryWidth
                       TempPotential(iLon,iLat,iDir) = TempWeight*TempPotential(iLon,iLat,iDir) + &
                            (1-TempWeight)*AMIEPotential(iLon,iLat,iDir)
                    endif
                 enddo
              enddo
             
           enddo
           ! set back to false for getting Weimer potential at the next iteration
           UAl_UseGridBasedEIE = .false.
        endif

        if (UseDynamo .and. .not. Is1D) then
           dynamo = 0.0
           call get_dynamo_potential( &
                MLongitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), &
                 MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), dynamo)

           do iDir = 1, nDir
              do iLon = -1,nLons+2
                 do iLat = -1,nLats+2
                    if (abs(MLatitude(iLon, iLat, iAlt, iBlock)) < DynamoHighLatBoundary) then
                       dis= (DynamoHighLatBoundary - &
                            abs(MLatitude(iLon, iLat, iAlt, iBlock)))/20.0
                       if (dis > 1.0) then
                          TempPotential(iLon,iLat,iDir) = dynamo(iLon,iLat)
                       else
                          TempPotential(iLon,iLat,iDir) = &
                               (1.0-dis) * TempPotential(iLon,iLat,iDir) + &
                               dis * dynamo(iLon,iLat)
                       endif
                    endif
                 enddo
              enddo
           enddo

        endif

        Potential(:,:,iAlt,iBlock) = TempPotential(:,:,1)
        if (UseTwoAMIEPotentials) then
           PotentialY(:,:,iAlt,iBlock) = TempPotential(:,:,2)
        else
           PotentialY(:,:,iAlt,iBlock) = Potential(:,:,iAlt,iBlock)
        endif

        !----------------------------------------------
        ! Another example of user output

        if (iAlt == 1) then 

           UserData2d(1:nLons,1:nLats,1,1,iBlock) = &
                TempPotential(1:nLons,1:nLats,1)/1000.0
        endif

     enddo

     IsFirstPotential(iBlock) = .false.

  endif

  ! -----------------------------------------------------
  ! Now get the aurora.
  ! This assumes that the field lines are basically
  ! vertical starting at the top of the model.
  ! -----------------------------------------------------
  
  if (floor((tSimulation-dt)/DtAurora) /= &
       floor((tsimulation)/DtAurora) .or. IsFirstAurora(iBlock)) then

     call report("Getting Aurora",1)

     iAlt = nAlts+1

     if (UseNewellAurora) then
        call run_newell(iBlock)
     elseif (UseOvationSME) then 
        call run_ovationsme(StartTime, CurrentTime, iBlock)
     elseif (UseAeModel) then
        call run_ae_model(CurrentTime, iBlock)
     else

        call UA_SetGrid(                    &
             MLT(-1:nLons+2,-1:nLats+2,iAlt), &
             MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), iError)

        if (iError /= 0) then
           write(*,*) "Error in routine get_potential (UA_SetGrid):"
           write(*,*) iError
           call stop_gitm("Stopping in get_potential")
        endif

        call UA_GetAveE(ElectronAverageEnergy, iError)
        if (iError /= 0) then
           write(*,*) "Error in get_potential (UA_GetAveE):"
           write(*,*) iError
           ElectronAverageEnergy = 1.0
        endif

        ! Sometimes, in AMIE, things get messed up in the
        ! Average energy, so go through and fix some of these.
        
        do iLat=-1,nLats+2
           do iLon=-1,nLons+2
              if (ElectronAverageEnergy(iLon,iLat) < 0.0) then
                 ElectronAverageEnergy(iLon,iLat) = 0.1
                 write(*,*) "ave e i,j Negative : ",iLon,iLat,&
                      ElectronAverageEnergy(iLon,iLat)
              endif
              if (ElectronAverageEnergy(iLon,iLat) > 100.0) then
                 write(*,*) "ave e i,j Positive : ",iLon,iLat,&
                      ElectronAverageEnergy(iLon,iLat)
                 ElectronAverageEnergy(iLon,iLat) = 0.1
              endif
           enddo
        enddo

        call UA_GetEFlux(ElectronEnergyFlux, iError)
        if (iError /= 0) then
           write(*,*) "Error in get_potential (UA_GetEFlux):"
           write(*,*) iError
           ElectronEnergyFlux = 0.1
        endif

        ! -----------------------------------------------------
        ! Get Ion Precipitation if desired
        ! -----------------------------------------------------

        if (UseIonPrecipitation) then

           call UA_GetIonAveE(IonAverageEnergy, iError)
           if (iError /= 0) then
              write(*,*) "Error in get_potential (UA_GetAveE):"
              write(*,*) iError
              IonAverageEnergy = 1.0
           endif

           call UA_GetIonEFlux(IonEnergyFlux, iError)
           if (iError /= 0) then
              write(*,*) "Error in get_potential (UA_GetEFlux):"
              write(*,*) iError
              IonEnergyFlux = 0.1
           endif

        endif
        
     endif

    if (UseCusp) then

       lats = abs(MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock))

       if (maxval(lats) > 50) then

          mlts = mod(MLT(-1:nLons+2,-1:nLats+2,iAlt)+24.0,24.0)

          call get_IMF_Bz(CurrentTime+TimeDelayHighLat, bz, iError)
          call get_IMF_By(CurrentTime+TimeDelayHighLat, by, iError)

          ! If we are in the southern hemisphere, reverse by:
          if (lats(nLons/2, nLats/2) < 0.0) by = -by

          if (bz > 0) then
             ! Newell et al., 1988:
             CuspLat = 77.2 + 0.11 * bz
             ! Asai et al., Earth Planets Space, 2005:
             CuspMlt = 11.755 + 0.169 * by
          else
             ! Asai et al., Earth Planets Space, 2005:
             CuspMlt = 11.949 + 0.0826 * by
             ! Zhang et al., JGR, 2005:
             if (Bz > -10) then
                CuspLat = 77.2 + 1.1 * bz
             else
                CuspLat = 21.7 * exp(0.1 * bz) + 58.2
             endif
          endif

          EFlux = CuspEFlux * &
               exp(-abs(lats - CuspLat)/CuspLatHalfWidth) *  &
               exp(-abs(mlts - CuspMlt)/CuspMltHalfWidth)

          do iLat=-1,nLats+2
             do iLon=-1,nLons+2
                if (EFlux(iLon,iLat) > 0.1) then
                   ElectronEnergyFlux(iLon,iLat) = EFlux(iLon,iLat)
                   ElectronAverageEnergy(iLon,iLat) = CuspAveE
                endif
             enddo
          enddo

       endif

    endif

     if (iDebugLevel > 2) &
          write(*,*) "==> Max, electron_ave_ene : ", &
          maxval(ElectronAverageEnergy), &
          maxval(ElectronEnergyFlux)

     IsFirstAurora(iBlock) = .false.

  endif

  if (iDebugLevel > 1) &
       write(*,*) "==> Min, Max, CPC Potential : ", &
       minval(Potential(:,:,:,iBlock))/1000.0, &
       maxval(Potential(:,:,:,iBlock))/1000.0, &
       (maxval(Potential(:,:,:,iBlock))-minval(Potential(:,:,:,iBlock)))/1000.0

  call end_timing("get_potential")

end subroutine get_potential


subroutine get_dynamo_potential(lons, lats, pot)

  use ModGITM
  use ModInputs, only: iDebugLevel, DynamoHighLatBoundary
  use ModElectrodynamics

  implicit none

  real, dimension(-1:nLons+2,-1:nLats+2), intent(in) :: lons
  real, dimension(-1:nLons+2,-1:nLats+2), intent(in) :: lats
  real, dimension(-1:nLons+2,-1:nLats+2), intent(out) :: pot

  integer :: iLon, iLat, iL, iM
  real    :: dM, dL, LatIn, LonIn, iError

  logical :: IsFound

  pot = 0.0

  IsFound = .false.

  do iLon = -1, nLons+2 
     do iLat = -1, nLats+2 

        LatIn = lats(iLon, iLat)
        LonIn = mod(lons(iLon, iLat)+360.0,360.0)

        IsFound = .false.

        if (abs(LatIn) <= MagLatMC(nMagLons,nMagLats)) then

           iM = 1
           do while (iM < nMagLons+1)
              iL = 1
              do while (iL < nMagLats)

                 !\
                 ! Check to see if the point is within the current cell
                 !/
                 
                 if ( LatIn <  MagLatMC(iM,iL+1) .and. &
                      LatIn >= MagLatMC(iM,iL) .and. &
                      LonIn <  MagLonMC(iM+1,iL) .and. &
                      LonIn >= MagLonMC(iM,iL)) then

                    dM=  (LonIn            -MagLonMC(iM,iL)) / &
                         (MagLonMC(iM+1,iL)-MagLonMC(iM,iL))

                    dL=  (LatIn        -MagLatMC(iM,iL))/&
                         (MagLatMC(iM,iL+1)-MagLatMC(iM,iL))

                    pot(iLon, iLat) = &
                         (1.0 - dM) * (1.0 - dL) * DynamoPotentialMC(iM  , iL  ) + &
                         (1.0 - dM) * (      dL) * DynamoPotentialMC(iM  , iL+1) + &
                         (      dM) * (      dL) * DynamoPotentialMC(iM+1, iL+1) + &
                         (      dM) * (1.0 - dL) * DynamoPotentialMC(iM+1, iL  )
                    

                    iL = nMagLats
                    iM = nMagLons

                    IsFound = .true.

                 endif

                 iL = iL + 1

              enddo

              iM = iM + 1

           enddo

           ! Check for near pole
           if (.not.IsFound) then 

              if (LatIn >  88.0) then
                 IsFound = .true.
                 pot(iLon, iLat) = sum(DynamoPotentialMC(:, nMagLats)) / (nMagLons+1)
              endif

              if (LatIn < -88.0) then
                 IsFound = .true.
                 pot(iLon, iLat) = sum(DynamoPotentialMC(:, 0)) / (nMagLons+1)
              endif

              write(*,*) "Inside the low latitude, but can't find the point!"
              write(*,*) LatIn, LonIn

           endif

        endif

        if (.not.IsFound) then 
           if (abs(LatIn) < MagLatMC(nMagLons, nMagLats)) &
                write(*,*) "=====> Could not find point : ", &
                LatIn, LonIn, DynamoHighLatBoundary, &
                MagLatMC(nMagLons, nMagLats)
           pot(iLon,iLat) = 0.0
        endif

     enddo

  enddo

end subroutine get_dynamo_potential
