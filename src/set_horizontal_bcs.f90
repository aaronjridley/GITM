! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine set_horizontal_bcs(iBlock)

  use ModSizeGitm
  use ModPlanet, only : nSpecies, nIons
  use ModInputs, only : LatStart, LatEnd, LonStart, LonEnd, GitmBCsDir, UseGitmBCs
  use ModReadGitm3d
  use ModGITM
  use ModTime, only : CurrentTime
  
  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iIon, iSpecies, iLat, iLon, iErr, iDir
  logical :: IsFirstTime = .true.

  character (len=nGitmVarCharLength), allocatable :: vars(:)
  ! --> declared in ModReadGitm3d --> real, allocatable :: GitmFileData(:,:)
  real, allocatable :: lonsBCs(:),latsBCs(:),altsBCs(:)
  integer :: nVars, nPoints, iPoint, iVar
  
  logical :: IsWestBC = .false.
  logical :: IsEastBC = .false.
  logical :: IsSouthBC = .false.
  logical :: IsNorthBC = .false.
  
  call report("set_horizontal_bcs",2)

  if (IsFirstTime .and. UseGitmBCs) then

     iErr = 0
     
     call GitmSetDir(GitmBCsDir)
     call GetGitmFileList(iErr)
     if (iErr /= 0) then
        call con_stop("Error in trying to read GITM 3D Filelist in horizontal bcs")
     endif

     call GetGitmGeneralHeaderInfo(iErr)
     if (iErr /= 0) then
        call con_stop("Error in reading gitm header information in horizontal bcs")
     endif

     call GitmGetnVars(nVars)
     if (iErr /= 0) then
        call con_stop("Error in getting number of variables in horizontal bcs")
     endif
    
     allocate(vars(nVars))
     call GitmGetVars(vars)

     ! Need to calculate the number of boundary points.

     nPoints = 0
     
     ! Western Boundary
     if (LonStart /= LonEnd .and. &
          minval(Longitude(:,iBlock)) < LonStart) then
        nPoints = nPoints + (nAlts+4) * (nLats+4) * 2
        IsWestBC = .true.
     endif
     
     ! Eastern Boundary
     if (LonStart /= LonEnd .and. &
          maxval(Longitude(:,iBlock)) > LonEnd) then
        nPoints = nPoints + (nAlts+4) * (nLats+4) * 2
        IsEastBC = .true.
     endif

     ! South Boundary
     if ( LatStart > -pi/2.0 .and. &
          minval(Latitude(:,iBlock)) < LatStart) then
        nPoints = nPoints + (nAlts+4) * (nLons+4) * 2
        IsSouthBC = .true.
     endif

     ! North Boundary
     if ( LatStart > -pi/2.0 .and. &
          maxval(Latitude(:,iBlocK)) > LatEnd) then
        nPoints = nPoints + (nAlts+4) * (nLons+4) * 2
        IsNorthBC = .true.
     endif

     call GitmSetnPointsToGet(nPoints)
     
     allocate(GitmFileData(nPoints, nVars))
     allocate(lonsBCs(nPoints))
     allocate(latsBCs(nPoints))
     allocate(altsBCs(nPoints))

     iPoint = 1
     
     if (IsWestBC) then
        do iLon = 0,-1,-1
           do iLat = -1, nLats+2
              do iAlt = -1, nAlts+2
                 lonsBCs(iPoint) = longitude(iLon,iBlock)
                 latsBCs(iPoint) = latitude(iLat,iBlock)
                 altsBCs(iPoint) = Altitude_GB(iLon,iLat,iAlt,iBlock)
                 iPoint = iPoint + 1
              enddo
           enddo
        enddo
     endif

     if (IsEastBC) then
        do iLon = nLons+1,nLons+2
           do iLat = -1, nLats+2
              do iAlt = -1, nAlts+2
                 lonsBCs(iPoint) = longitude(iLon,iBlock)
                 latsBCs(iPoint) = latitude(iLat,iBlock)
                 altsBCs(iPoint) = Altitude_GB(iLon,iLat,iAlt,iBlock)
                 iPoint = iPoint + 1
              enddo
           enddo
        enddo
     endif

     if (IsSouthBC) then
        do iLon = -1,nLons+2
           do iLat = 0,-1,-1
              do iAlt = -1, nAlts+2
                 lonsBCs(iPoint) = longitude(iLon,iBlock)
                 latsBCs(iPoint) = latitude(iLat,iBlock)
                 altsBCs(iPoint) = Altitude_GB(iLon,iLat,iAlt,iBlock)
                 iPoint = iPoint + 1
              enddo
           enddo
        enddo
     endif

     if (IsNorthBC) then
        do iLon = -1,nLons+2
           do iLat = nLats+1,nLats+2
              do iAlt = -1, nAlts+2
                 lonsBCs(iPoint) = longitude(iLon,iBlock)
                 latsBCs(iPoint) = latitude(iLat,iBlock)
                 altsBCs(iPoint) = Altitude_GB(iLon,iLat,iAlt,iBlock)
                 iPoint = iPoint + 1
              enddo
           enddo
        enddo
     endif

     if (iPoint > 1) then 

        ! Convert to degrees and km
        lonsBCs = lonsBCs * 360.0 / twopi
        latsBCs = latsBCs * 360.0 / twopi
        altsBCs = altsBCs / 1000.0
     
        call GitmSetGrid(lonsBCs,latsBCs,altsBCs)

     endif

     IsFirstTime = .false.

  endif
  
  if (UseGitmBCs) then 

     call GitmUpdateTime(CurrentTime, iErr)
     if (iErr == 0) then 
        call GitmGetData(GitmFileData)
     else
        write(*,*) 'Error in getting GITM data in set_horizontal_BCs!'
        call CON_stop('Must Stop!')
     endif
          
     iPoint = 1

     if (IsWestBC) then
        do iLon = 0,-1,-1
           do iLat = -1, nLats+2
              do iAlt = -1, nAlts+2
                 call set_horizontal_bcs_1point(iLon,iLat,iAlt,iBlock,iPoint)
                 iPoint = iPoint + 1
              enddo
           enddo
        enddo
     endif
              
     if (IsEastBC) then
        do iLon = nLons+1,nLons+2
           do iLat = -1, nLats+2
              do iAlt = -1, nAlts+2
                 call set_horizontal_bcs_1point(iLon,iLat,iAlt,iBlock,iPoint)
                 iPoint = iPoint + 1
              enddo
           enddo
        enddo
     endif
              
     if (IsSouthBC) then
        do iLon = -1,nLons+2
           do iLat = 0,-1,-1
              do iAlt = -1, nAlts+2
                 call set_horizontal_bcs_1point(iLon,iLat,iAlt,iBlock,iPoint)
                 iPoint = iPoint + 1
              enddo
           enddo
        enddo
     endif

     if (IsNorthBC) then
        do iLon = -1,nLons+2
           do iLat = nLats+1,nLats+2
              do iAlt = -1, nAlts+2
                 call set_horizontal_bcs_1point(iLon,iLat,iAlt,iBlock,iPoint)
                 iPoint = iPoint + 1
              enddo
           enddo
        enddo
     endif

  else

     ! Western Boundary

     if (LonStart /= LonEnd .and. &
          minval(Longitude(:,iBlock)) < LonStart) then

        do iAlt = -1, nAlts+2
           do iLon = 0,-1,-1

              do iSpecies = 1, nSpecies
                 VerticalVelocity(iLon,:,iAlt,iSpecies,iBlock) = &
                      VerticalVelocity(iLon+1,:,iAlt,iSpecies,iBlock)
                 nDensityS(iLon,:,iAlt,iSpecies,iBlock) = &
                      nDensityS(iLon+1,:,iAlt,iSpecies,iBlock)
              enddo
              
              Rho(iLon, :,iAlt,iBlock) = &
                   Rho(iLon+1,:,iAlt,iBlock)
              Temperature( iLon,:,iAlt,iBlock)= &
                   Temperature(iLon+1,:,iAlt,iBlock)
           iTemperature(iLon,:,iAlt,iBlock) = &
                iTemperature(iLon+1,:,iAlt,iBlock)
           eTemperature(iLon,:,iAlt,iBlock) = &
                eTemperature(iLon+1,:,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(iLon,:,iAlt,iIon,iBlock) = &
                   IDensityS(iLon+1,:,iAlt,iIon,iBlock)
           enddo

           Velocity(iLon,:,iAlt,iNorth_,iBlock) = &
                Velocity(iLon+1,:,iAlt,iNorth_,iBlock)
           Velocity(iLon,:,iAlt,iEast_,iBlock) = &
                Velocity(iLon+1,:,iAlt,iEast_,iBlock)
           Velocity(iLon,:,iAlt,iUp_,iBlock) = &
                Velocity(iLon+1,:,iAlt,iUp_,iBlock)

           iVelocity(iLon,:,iAlt,iNorth_,iBlock) = &
                iVelocity(iLon+1,:,iAlt,iNorth_,iBlock)
           iVelocity(iLon,:,iAlt,iEast_,iBlock) = &
                iVelocity(iLon+1,:,iAlt,iEast_,iBlock)
           iVelocity(iLon,:,iAlt,iUp_,iBlock) = &
                iVelocity(iLon+1,:,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif

  if (LonStart /= LonEnd .and. &
       maxval(Longitude(:,iBlock)) > LonEnd) then

     do iAlt = -1, nAlts+2
        do iLon = nLons+1,nLons+2

           do iSpecies = 1, nSpecies
              VerticalVelocity(iLon,:,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(iLon-1,:,iAlt,iSpecies,iBlock)
              nDensityS(iLon,:,iAlt,iSpecies,iBlock) = &
                   nDensityS(iLon-1,:,iAlt,iSpecies,iBlock)
!              nDensityS(iLon,:,iAlt,iSpecies,iBlock) = &
!                   2*nDensityS(iLon-1,:,iAlt,iSpecies,iBlock) - &
!                   nDensityS(iLon-2,:,iAlt,iSpecies,iBlock)
           enddo

           Rho(iLon, :,iAlt,iBlock) = &
                Rho(iLon-1,:,iAlt,iBlock)
           Temperature( iLon,:,iAlt,iBlock)= &
                Temperature(iLon-1,:,iAlt,iBlock)
           iTemperature(iLon,:,iAlt,iBlock) = &
                iTemperature(iLon-1,:,iAlt,iBlock)
           eTemperature(iLon,:,iAlt,iBlock) = &
                eTemperature(iLon-1,:,iAlt,iBlock)

!           Rho(iLon, :,iAlt,iBlock) = &
!                2*Rho(iLon-1,:,iAlt,iBlock) - &
!                Rho(iLon-2,:,iAlt,iBlock)
!           Temperature( iLon,:,iAlt,iBlock)= &
!                2*Temperature(iLon-1,:,iAlt,iBlock) - &
!                Temperature(iLon-2,:,iAlt,iBlock)
!           iTemperature(iLon,:,iAlt,iBlock) = &
!                2*iTemperature(iLon-1,:,iAlt,iBlock) - &
!                iTemperature(iLon-2,:,iAlt,iBlock)
!           eTemperature(iLon,:,iAlt,iBlock) = &
!                2*eTemperature(iLon-1,:,iAlt,iBlock) - &
!                eTemperature(iLon-2,:,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(iLon,:,iAlt,iIon,iBlock) = &
                   IDensityS(iLon-1,:,iAlt,iIon,iBlock)
!              IDensityS(iLon,:,iAlt,iIon,iBlock) = &
!                   2*IDensityS(iLon-1,:,iAlt,iIon,iBlock) - &
!                   IDensityS(iLon-2,:,iAlt,iIon,iBlock)
           enddo

           Velocity(iLon,:,iAlt,iNorth_,iBlock) = &
                Velocity(iLon-1,:,iAlt,iNorth_,iBlock)
           Velocity(iLon,:,iAlt,iEast_,iBlock) = &
                Velocity(iLon-1,:,iAlt,iEast_,iBlock)
           Velocity(iLon,:,iAlt,iUp_,iBlock) = &
                Velocity(iLon-1,:,iAlt,iUp_,iBlock)

           iVelocity(iLon,:,iAlt,iNorth_,iBlock) = &
                iVelocity(iLon-1,:,iAlt,iNorth_,iBlock)
           iVelocity(iLon,:,iAlt,iEast_,iBlock) = &
                iVelocity(iLon-1,:,iAlt,iEast_,iBlock)
           iVelocity(iLon,:,iAlt,iUp_,iBlock) = &
                iVelocity(iLon-1,:,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif


  ! Southern Boundary

  if ( LatStart > -pi/2.0 .and. &
       minval(Latitude(:,iBlock)) < LatStart) then

     do iAlt = -1, nAlts+2
        do iLat = 0,-1,-1

           do iSpecies = 1, nSpecies
              VerticalVelocity(:,iLat,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(:,iLat+1,iAlt,iSpecies,iBlock)
              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
                   nDensityS(:,iLat+1,iAlt,iSpecies,iBlock)
!              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
!                   2*nDensityS(:,iLat+1,iAlt,iSpecies,iBlock)- &
!                   nDensityS(:,iLat+2,iAlt,iSpecies,iBlock)
           enddo

           Rho(:, iLat,iAlt,iBlock) = &
                Rho(:,iLat+1,iAlt,iBlock)
           Temperature( :, iLat,iAlt,iBlock)= &
                Temperature(:,iLat+1,iAlt,iBlock)
!           Rho(:, iLat,iAlt,iBlock) = &
!                2*Rho(:,iLat+1,iAlt,iBlock) - &
!                Rho(:,iLat+2,iAlt,iBlock)
!           Temperature( :, iLat,iAlt,iBlock)= &
!                2*Temperature(:,iLat+1,iAlt,iBlock)- &
!                Temperature(:,iLat+2,iAlt,iBlock)
           iTemperature(:, iLat,iAlt,iBlock)=iTemperature(:,iLat+1,iAlt,iBlock)
           eTemperature(:, iLat,iAlt,iBlock)=eTemperature(:,iLat+1,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(:,iLat,iAlt,iIon,iBlock) = &
                   IDensityS(:,iLat+1,iAlt,iIon,iBlock)
           enddo

           Velocity(:,iLat,iAlt,iNorth_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iNorth_,iBlock)
           Velocity(:,iLat,iAlt,iEast_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iEast_,iBlock)
           Velocity(:,iLat,iAlt,iUp_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iUp_,iBlock)

           iVelocity(:,iLat,iAlt,iNorth_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iNorth_,iBlock)
           iVelocity(:,iLat,iAlt,iEast_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iEast_,iBlock)
           iVelocity(:,iLat,iAlt,iUp_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif

  if (LatEnd < pi/2.0 .and. &
       maxval(Latitude(:,iBlocK)) > LatEnd) then

     do iAlt = -1, nAlts+2
        do iLat = nLats+1, nLats+2

           do iSpecies = 1, nSpecies
              VerticalVelocity(:,iLat,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(:,iLat-1,iAlt,iSpecies,iBlock)
              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
                   nDensityS(:,iLat-1,iAlt,iSpecies,iBlock)
!              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
!                   2*nDensityS(:,iLat-1,iAlt,iSpecies,iBlock) - &
!                   nDensityS(:,iLat-2,iAlt,iSpecies,iBlock)
           enddo

           Rho(:, iLat,iAlt,iBlock) = &
                Rho(:,iLat-1,iAlt,iBlock)
           Temperature( :, iLat,iAlt,iBlock)= &
                Temperature(:,iLat-1,iAlt,iBlock)
!           Rho(:, iLat,iAlt,iBlock) = &
!                2*Rho(:,iLat-1,iAlt,iBlock) - &
!                Rho(:,iLat-2,iAlt,iBlock)
!           Temperature( :, iLat,iAlt,iBlock)= &
!                2*Temperature(:,iLat-1,iAlt,iBlock) - &
!                Temperature(:,iLat-2,iAlt,iBlock)
           iTemperature(:, iLat,iAlt,iBlock)=iTemperature(:,iLat-1,iAlt,iBlock)
           eTemperature(:, iLat,iAlt,iBlock)=eTemperature(:,iLat-1,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(:,iLat,iAlt,iIon,iBlock) = &
                   IDensityS(:,iLat-1,iAlt,iIon,iBlock)
           enddo

           Velocity(:,iLat,iAlt,iNorth_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iNorth_,iBlock)
           Velocity(:,iLat,iAlt,iEast_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iEast_,iBlock)
           Velocity(:,iLat,iAlt,iUp_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iUp_,iBlock)

           iVelocity(:,iLat,iAlt,iNorth_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iNorth_,iBlock)
           iVelocity(:,iLat,iAlt,iEast_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iEast_,iBlock)
           iVelocity(:,iLat,iAlt,iUp_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif

endif
  
end subroutine set_horizontal_bcs

subroutine set_horizontal_bcs_1point(iLon,iLat,iAlt,iBlock,iPoint)

  use ModPlanet
  use ModGITM
  use ModReadGitm3d
  
  implicit none
  
  integer, intent(in) :: iLon,iLat,iAlt,iBlock,iPoint
  integer :: iDir,iSpecies
  
  Rho(iLon, iLat, iAlt, iBlock) = GitmFileData(iPoint,iRho_)

  do iSpecies = 1, nSpeciesTotal
     NDensityS(iLon,iLat,iAlt,iSpecies,iBlock) = &
          GitmFileData(iPoint,iSpecies+iNeutralStart_-1)
!     if (iProc == 0) then
!        write(*,*) 'ndensity : ',iSpecies, GitmFileData(iPoint,iSpecies+iNeutralStart_-1)
!     endif
  enddo
  do iSpecies = 1, nIons
     IDensityS(iLon,iLat,iAlt,iSpecies,iBlock) = &
          GitmFileData(iPoint,iSpecies+iIonStart_-1)
!     if (iProc == 0) then
!        write(*,*) 'idensity : ',iSpecies, GitmFileData(iPoint,iSpecies+iIonStart_-1)
!     endif
  enddo

  ! Need to calculate tempunit to get normalized temperature
  MeanMajorMass(iLon,iLat,iAlt) = 0.0
  do iSpecies = 1, nSpeciesTotal
     MeanMajorMass(iLon,iLat,iAlt) = &
          MeanMajorMass(iLon,iLat,iAlt) + &
          Mass(iSpecies) * NDensityS(iLon,iLat,iAlt,iSpecies,iBlock) / &
          sum(NDensityS(iLon,iLat,iAlt,:,iBlock))
  enddo
  TempUnit(iLon,iLat,iAlt) = MeanMajorMass(iLon,iLat,iAlt) / Boltzmanns_Constant
  Temperature(iLon,iLat,iAlt,iBlock) = &
       GitmFileData(iPoint,iTn_) / TempUnit(iLon,iLat,iAlt)

  eTemperature(iLon,iLat,iAlt,iBlock) = GitmFileData(iPoint,iTe_)
  iTemperature(iLon,iLat,iAlt,iBlock) = GitmFileData(iPoint,iTe_+1)
  ! Velocities
  do iDir = 1, 3
     Velocity(iLon,iLat,iAlt,iDir,iBlock) = &
          GitmFileData(iPoint,iVn_+iDir-1)
     iVelocity(iLon,iLat,iAlt,iDir,iBlock) = &
          GitmFileData(iPoint,iVi_+iDir-1)
  enddo
  ! Vertical Velocities
  do iSpecies = 1, nSpecies
     VerticalVelocity(iLon,iLat,iAlt,iSpecies,iBlock) = &
          GitmFileData(iPoint,iVn_+3+iSpecies-1)
  enddo

end subroutine set_horizontal_bcs_1point
  
