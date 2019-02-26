
!---------------------------------------------------------------------------
! Routines for making GITM use another model of the ionosphere
!---------------------------------------------------------------------------

subroutine overwrite_ionosphere

  use ModInputs
  use ModSizeGitm
  
  implicit none

  integer :: iBlock

  write(*,*) 'Overwriting Ionosphere!'
  
  if (DoOverwriteWithIRI) then
     ! This doesn't include the velocities
     call init_iri
  else
     do iBlock = 1, nBlocks
        call ionosphere_overwrite_sami(iBlock)
     enddo
  endif
  
end subroutine overwrite_ionosphere

!---------------------------------------------------------------------------
! Routine for overwriting GITM data with SAMI data
!---------------------------------------------------------------------------

subroutine ionosphere_overwrite_sami(iBlock)

  use ModReadSami3D

  use ModSizeGitm
  use ModPlanet, only : nSpecies, nIons
  use ModInputs
  use ModReadSami3d
  use ModGITM
  use ModTime, only : CurrentTime
  
  implicit none

  save
  
  integer, intent(in) :: iBlock
  character (len=nSamiVarCharLength), allocatable :: SamiVars(:)
  real, allocatable :: GitmLons(:),GitmLats(:),GitmAlts(:)
  real, allocatable :: SamiFileData(:,:)
  logical :: IsFirstTime = .true.
  
  integer :: nPoints, nVars, iErr
  integer :: iPoint, iLon, iLat, iAlt

  call report("ionosphere_overwrite_sami",2)

  if (IsFirstTime) then

     iErr = 0

     call SamiReadInputFile(SamiInFile)
     call SamiGetnVars(nVars)

     nPoints = (nLons+4) * (nLats+4) * (nAlts)
     call SamiSetnPointsToGet(nPoints)
     
     allocate(SamiVars(nVars))
     allocate(SamiFileData(nPoints,nVars))

     call SamiGetVars(SamiVars)

     allocate(GitmLons(nPoints))
     allocate(GitmLats(nPoints))
     allocate(GitmAlts(nPoints))

     iPoint = 1

     do iLon = -1,nLons+2
        do iLat = -1, nLats+2
           do iAlt = 1, nAlts
              GitmLons(iPoint) = longitude(iLon,iBlock)
              GitmLats(iPoint) = latitude(iLat,iBlock)
              GitmAlts(iPoint) = Altitude_GB(iLon,iLat,iAlt,iBlock)
              iPoint = iPoint + 1
           enddo
        enddo
     enddo

     if (iPoint > 1) then 

        ! Convert to degrees and km
        GitmLons = GitmLons * 360.0 / twopi
        GitmLats = GitmLats * 360.0 / twopi
        GitmAlts = GitmAlts / 1000.0
     
        call SamiSetGrid(GitmLons,GitmLats,GitmAlts)

     endif
     
     IsFirstTime = .false.
     
  endif

  call SamiUpdateTime(CurrentTime, iErr)
  if (iErr == 0) then 
     call SamiGetData(SamiFileData)
  endif

  iPoint = 1

  do iLon = -1,nLons+2
     do iLat = -1, nLats+2
        do iAlt = 1, nAlts
           !write(*,*) iLon,iLat,iAlt,iPoint
           IDensityS(iLon,iLat,iAlt,iBlock,iO_4SP_) = SamiFileData(iPoint, iSami_Op_)*1e6
           IDensityS(iLon,iLat,iAlt,iBlock,iO2P_) = SamiFileData(iPoint, iSami_O2p_)*1e6
           IDensityS(iLon,iLat,iAlt,iBlock,iNOP_) = SamiFileData(iPoint, iSami_NOp_)*1e6
           IDensityS(iLon,iLat,iAlt,iBlock,iHP_) = SamiFileData(iPoint, iSami_Hp_)*1e6
           IDensityS(iLon,iLat,iAlt,iBlock,iHeP_) = SamiFileData(iPoint, iSami_Hep_)*1e6

           IDensityS(iLon,iLat,iAlt,iBlock,ie_) = &
                sum(IDensityS(iLon,iLat,iAlt,iBlock,1:nIons-1))

           eTemperature(iLon,iLat,iAlt,iBlock) = SamiFileData(iPoint, iSami_Te_)
           iTemperature(iLon,iLat,iAlt,iBlock) = SamiFileData(iPoint, iSami_Ti_)
           
           iPoint = iPoint + 1
        enddo
     enddo
  enddo


  
end subroutine ionosphere_overwrite_sami
