!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!\
! ------------------------------------------------------------
! calc_timestep
! ------------------------------------------------------------
!/

subroutine calc_timestep_horizontal

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi
  implicit none

  logical :: DoWarning = .true.

  real :: DtCond, DtLocal, DtHorizontal, DtEnd

  real :: cSound_H(0:nLons+1,0:nLats+1)

  integer :: iBlock, iError, iLat, iLon, iAlt

  call report("calc_timestep_horizontal",2)

  DtLocal = 1.0e32

  do iBlock = 1, nBlocks

     do iAlt = 1,nAlts

        ! Calculate maximum propagation speeds for the horizontal directions
        cSound_H = sqrt(Gamma(0:nLons+1,0:nLats+1,iAlt,iBlock) * &
                   Temperature(0:nLons+1,0:nLats+1,iAlt,iBlock))

        cMax_GDB(:,:,iAlt,iNorth_,iBlock) = &
             abs(Velocity(0:nLons+1,0:nLats+1,iAlt,iNorth_,iBlock)) + cSound_H

        cMax_GDB(:,:,iAlt,iEast_,iBlock) = &
             abs(Velocity(0:nLons+1,0:nLats+1,iAlt,iEast_,iBlock)) + cSound_H

        ! Find stability limit on the time step
        do iLat = 1,nLats
           do iLon = 1,nLons
           
              DtLocal = min(DtLocal, Cfl / ( &
                   cMax_GDB(iLon, iLat, iAlt, iEast_,  iBlock) / &
                   dLonDist_GB(iLon,iLat,iAlt,iBlock) + &
                   cMax_GDB(iLon, iLat, iAlt, iNorth_, iBlock) / &
                   dLatDist_GB(iLon,iLat,iAlt,iBlock)))

              if (UseIonAdvection) &
                   DtLocal = min(DtLocal, Cfl / ( &
                   (abs(IVelocity(iLon, iLat, iAlt, iEast_,  iBlock))+ &
                   cSound_H(iLon, iLat))/ &
                   dLonDist_GB(iLon,iLat,iAlt,iBlock) + &
                   (abs(IVelocity(iLon, iLat, iAlt, iNorth_, iBlock))+ &
                   cSound_H(iLon, iLat))/ &
                   dLatDist_GB(iLon,iLat,iAlt,iBlock))/2.0)

           end do
        end do

     end do

  end do

  DtEnd = EndTime - CurrentTime
  DtLocal = min(DtLocal, DtEnd)

  call MPI_AllREDUCE(DtLocal, DtHorizontal, 1, MPI_REAL, MPI_MIN, &
       iCommGITM, iError)

  if (iDebugLevel > 2) &
       write(*,*) "===> DtVertical, DtHorizontal : ", Dt, DtHorizontal

  Dt = min(Dt, DtHorizontal)

end subroutine calc_timestep_horizontal

!==========================================================================

subroutine calc_timestep_vertical

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi
  implicit none

  logical :: DoWarning = .true.

  real :: DtCond, DtLocal, DtVertical, DtEnd

  real :: cm(1:nAlts)

  integer :: iBlock, iError, iLat, iLon, iAlt
  integer :: iLatSave, iLonSave

  call report("calc_timestep_vertical",2)

  DtLocal = 1.0e32

  iLatSave = 0
  iLonSave = 0

  do iBlock = 1, nBlocks

     do iLon = 1,nLons
        do iLat = 1,nLats

           do iAlt = 0, nAlts+1
              cMax_GDB(iLon,iLat,iAlt,iUp_,iBlock) = &
                   maxval(abs(VerticalVelocity(iLon,iLat,iAlt,:,iBlock))) + &
                   sqrt(Gamma(iLon,iLat,iAlt,iBlock) * Temperature(iLon,iLat,iAlt,iBlock))
           enddo

           DtLocal = min(DtLocal, &
                Cfl / &
                maxval(cMax_GDB(iLon,iLat,1:nAlts,iUp_,iBlock)&
                /      dAlt_GB(iLon,iLat,1:nAlts,iBlock)))

           if (UseIonAdvection) then

              cm = abs(IVelocity(iLon,iLat,1:nAlts,iUp_,iBlock)) + &
                sqrt(Gamma(iLon,iLat,1:nAlts,iBlock) * Temperature(iLon,iLat,1:nAlts,iBlock))

              DtLocal = min(DtLocal, &
                   Cfl / &
                   maxval(cm/dAlt_GB(iLon,iLat,1:nAlts,iBlock)))

           endif

        enddo
     enddo

  enddo

  DtEnd = EndTime - CurrentTime
  dTLocal = min(DtLocal, DtEnd)

  call MPI_AllREDUCE(dTLocal, DtVertical, 1, MPI_REAL, MPI_MIN, &
       iCommGITM, iError)

  Dt = min(Dt, DtVertical)

  if (iDebugLevel > 2) &
       write(*,*) "===> DtVertical : ", Dt

  if (dt < cfl/100.0 .and. dt < 0.99*DtEnd) then
     write(*,*) "Dt too slow!!!", dt
     call stop_gitm("Stopping in calc_timestep_vertical")
  endif

end subroutine calc_timestep_vertical

