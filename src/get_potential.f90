!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine init_get_potential

  use ModGITM
  use ModTime
  use ModIndicesInterfaces
  use ModInputs
  use ModUserGITM
  use ModNewell
  use ModOvationSME

  implicit none

  character (len=100), dimension(100) :: Lines
  character (len=100) :: TimeLine
  real    :: bz

!  real    :: dynamo(-1:nLons+2,-1:nLats+2)
  logical :: IsFirstTime = .true.
  integer :: iError
  iError = 0

  if (.not.IsFirstTime .or. IsFramework) return

  IsFirstTime = .false.

  if (UseNewellAurora) then
     call init_newell
     UseIMF = .true.
  endif

  if (UseOvationSME) then
     call read_ovationsm_files
  endif

  if (index(cAMIEFileNorth,"none") > 0) then

     Lines(1) = "#BACKGROUND"
     Lines(2) = "EIE/"

     UseHPI = .true.
     if (UseNewellAurora .or. UseOvationSME) UseHPI = .false.

     call get_IMF_Bz(CurrentTime+TimeDelayHighLat, bz, iError)

     call IO_SetIMFBz(bz)
     if (iError /= 0) then
!        write(*,*) "Can not find IMF Bz."
!        write(*,*) "Setting potential to Millstone HPI."
        Lines(3) = "millstone_hpi"    ! Change to "zero" if you want
     else
!        write(*,*) "Setting potential to ",PotentialModel
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

  call EIE_set_inputs(Lines)

  call EIE_Initialize(iError)

  if (iError /= 0) then
     write(*,*) &
          "Code Error in IE_Initialize called from get_potential.f90"
     write(*,*) "Error : ",iError
     call stop_gitm("Stopping in get_potential")
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

     call read_NOAAHPI_Indices_new(iError, CurrentTime+TimeDelayHighLat, EndTime+TimeDelayHighLat)
     call read_MHDIMF_Indices_new(iError, CurrentTime+TimeDelayHighLat, EndTime+TimeDelayHighLat)
     call get_IMF_Bz(CurrentTime+TimeDelayHighLat, bz, iError)
     if (bz < -20.0) bz = -20.0
     if (bz >  20.0) bz =  20.0
     call IO_SetIMFBz(bz)

     if (iError /= 0) then
        write(*,*) &
             "Code Error in get_IMF_Bz called from get_potential.f90"
        write(*,*) "Code : ",iError
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 1) write(*,*) "==> IMF Bz : ",bz

     call get_IMF_By(CurrentTime+TimeDelayHighLat, by, iError)
     if (by < -20.0) by = -20.0
     if (by >  20.0) by =  20.0
     call IO_SetIMFBy(by)

     if (iError /= 0) then
        write(*,*) &
             "Code Error in get_IMF_By called from get_potential.f90"
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 1) write(*,*) "==> IMF By : ",by

     call get_SW_V(CurrentTime+TimeDelayHighLat, vx, iError)
     if (vx < -800.0) vx = -800.0
     if (vx >  800.0) vx =  800.0
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

!  if (index(cAMIEFileNorth,"none") <= 0 .and. &
!     index(cAMIEFileNorth,"SPS") <= 0 .and. iBlock == 1) then 
  if (index(cAMIEFileNorth,"none") <= 0 .and. &
     index(cAMIEFileNorth,"SPS") <= 0) then 
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

  implicit none

  integer, intent(in) :: iBlock

  integer :: iError, iLat, iLon, iAlt
  logical :: IsFirstTime = .true.
  logical :: IsFirstPotential(nBlocksMax) = .true.
  logical :: IsFirstAurora(nBlocksMax) = .true.
  real    :: mP, dis

  real, dimension(-1:nLons+2,-1:nLats+2) :: TempPotential, Grid, dynamo
  real, dimension(-1:nLons+2,-1:nLats+2) :: SubMLats, SubMLons


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

     if (iDebugLevel > 1) write(*,*) "==> Setting up IE Grid"

     call init_get_potential
     call UA_SetnMLTs(nLons+4)
     call UA_SetnLats(nLats+4)
     if (.not. IsFramework) then 
        call IO_SetTime(CurrentTime)
        call set_indices
     endif
     call UA_SetNorth

     if (iDebugLevel > 1) write(*,*) "==> Getting Potential"

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

        call UA_GetPotential(TempPotential, iError)

        if (iError /= 0) then
           write(*,*) "Error in get_potential (UA_GetPotential):"
           write(*,*) iError
           TempPotential = 0.0
!           call stop_gitm("Stopping in get_potential")
        endif

        if (UseDynamo .and. .not. Is1D) then
           dynamo = 0.0
           call get_dynamo_potential( &
                MLongitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), &
                 MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), dynamo)

           do iLon = -1,nLons+2
              do iLat = -1,nLats+2 
!!                 if (abs(MLatitude(iLon, iLat, iAlt, iBlock)) < DynamoHighLatBoundary) then
!!                    TempPotential(iLon,iLat) = TempPotential(iLon,iLat) + dynamo(iLon,iLat)
                 if (abs(MLatitude(iLon, iLat, iAlt, iBlock)) < DynamoHighLatBoundary) then
                    dis= (DynamoHighLatBoundary - &
                          abs(MLatitude(iLon, iLat, iAlt, iBlock)))/20.0
                    if (dis > 1.0) then
                       TempPotential(iLon,iLat) = dynamo(iLon,iLat)
                    else
                       TempPotential(iLon,iLat) = &
                            (1.0-dis) * TempPotential(iLon,iLat) + &
                            dis * dynamo(iLon,iLat)
                    endif
                 endif
              enddo
           enddo

        endif

        Potential(:,:,iAlt,iBlock) = TempPotential

        !----------------------------------------------
        ! Another example of user output

        if (iAlt == 1) then 

           UserData2d(1:nLons,1:nLats,1,1,iBlock) = &
                TempPotential(1:nLons,1:nLats)/1000.0
        endif

     enddo

     IsFirstPotential(iBlock) = .false.

  endif

  if (floor((tSimulation-dt)/DtAurora) /= &
       floor((tsimulation)/DtAurora) .or. IsFirstAurora(iBlock)) then

     if (iDebugLevel > 1) write(*,*) "==> Getting Aurora"

     iAlt = nAlts+1

     if (UseNewellAurora) then
        call run_newell(iBlock)
     elseif (UseOvationSME) then 
        call run_ovationsme(StartTime, CurrentTime, iBlock)
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
!           call stop_gitm("Stopping in get_potential")
           ElectronAverageEnergy = 1.0
        endif

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
!           call stop_gitm("Stopping in get_potential")
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
