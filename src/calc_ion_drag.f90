! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_ion_drag(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  implicit none

  integer, intent(in) :: iBlock
  integer :: iErr, iSpecies, iDir

  real :: tmp(nLons, nLats, nAlts)
  real :: RhoI(nLons, nLats, nAlts)
  
  if (iDebugLevel > 4) write(*,*) "=====> Ion Drag", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM, iErr)
  
  RhoI = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock) * &
       MeanIonMass(1:nLons,1:nLats,1:nAlts)

  !
  ! This method should work, but seems to be way, way too strong, and I don't know why:
  !

!  if (UseIonDrag) then
!
!     do iSpecies = 1, nSpecies
!
!        Weight = Mass(iSpecies) * NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) / &
!             Rho(1:nLons,1:nLats,1:nAlts,iBlock)
!
!        do iIon = 1, nIons-1
!
!           ! Generalize to use all collisions, AJR - Nov. 3, 2019
!
!           tmp = &
!                IonCollisions(1:nLons,1:nLats,1:nAlts,iIon,iSpecies) * &
!                IDensityS(1:nLons,1:nLats,1:nAlts,iIon,iBlock) * MassI(iIon) / &
!                (NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) * Mass(iSpecies))
!
!           do iDir = 1, 3
!              IonDrag(:,:,:,iDir) = IonDrag(:,:,:,iDir) + &
!                   tmp * &
!                   (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
!                   Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock)) * &
!                   Weight
!           enddo
!
!           VerticalIonDrag(:,:,:,iSpecies) = VerticalIonDrag(:,:,:,iSpecies) + &
!                tmp * &
!                (IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) - &
!                VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock))
!
!        enddo
!     enddo
!
!  endif

  tmp = Collisions(1:nLons, 1:nLats, 1:nAlts, iVIN_)*&
       RhoI/Rho(1:nLons, 1:nLats, 1:nAlts, iBlock)

  do iDir = 1, 3
     IonDrag(:,:,:,iDir) = tmp * &
          (IVelocity(1:nLons, 1:nLats, 1:nAlts, iDir, iBlock) - &
          Velocity(1:nLons, 1:nLats, 1:nAlts, iDir, iBlock))
  enddo

  ! F_iondrag = rho_i/rho * Vis * (Ui-Un)
  ! where Vis = Vin *(Ns/N)  

  do iSpecies = 1, nSpecies

     tmp = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
          RhoI / &
          (Mass(iSpecies) * &
          NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock)) * &
          (NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) / &
          NDensity(1:nLons,1:nLats,1:nAlts,iBlock))

     VerticalIonDrag(:,:,:,iSpecies) = tmp * &
          (IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) - &
          VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock))

  enddo
  
  
end subroutine calc_ion_drag
  
