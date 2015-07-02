!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine advance

  use ModConstants
  use ModGITM
  use ModTime
  use ModEUV
  use ModInputs

  implicit none

  integer :: iBlock, iAlt, iLat, iLon

  real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2) :: &
       NewLogRho, NewVel, NewTemp

  real, dimension(1:nLons,1:nLats,1:nAlts) ::           &
       GradLogRho, DivVel, GradTemp,                    &
       DiffLogRho, DiffVel, DiffTemp, GradTmp, DiffTmp

  real, dimension(1:nLons,1:nLats,1:nAlts+1) ::           &
       Conduction

  call start_timing("advance")

  do iBlock = 1, nBlocks

     call calc_physics(iBlock)
     call chapman_integrals(iBlock)
     call calc_rates(iBlock)
     call euv_ionization_heat(iBlock)

     if (floor((tSimulation-dt)/DtPlot) /= floor((tsimulation)/DtPlot)) &
          call output("data/",".dat",iBlock)

     call rusanov(LogRho(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock), &
          GradLogRho, DiffLogRho, iBlock)
     call rusanov(Temperature(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock), &
          GradTemp,   DiffTemp, iBlock)
     call rusanov(Velocity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iUp_,iBlock),&
          DivVel,     DiffVel, iBlock)

     do iLon = 1, nLons
        do iLat = 1, nLats
           do iAlt = 1, nAlts

              NewLogRho(iLon,iLat,iAlt) = LogRho(iLon,iLat,iAlt,iBlock) &
                   - Dt * &
                   (DivVel(iLon,iLat,iAlt) + &
                   Velocity(iLon,iLat,iAlt,iUp_,iBlock) * &
                   GradLogRho(iLon,iLat,iAlt) ) &
                   + Dt * DiffLogRho(iLon,iLat,iAlt)

              NewVel(iLon,iLat,iAlt)    = Velocity(iLon,iLat,iAlt,iUp_,iBlock)&
                   - Dt * &
                   (Velocity(iLon,iLat,iAlt,iUp_,iBlock)* &
                   DivVel(iLon,iLat,iAlt) &
                   + GradTemp(iLon,iLat,iAlt) &
                   + Temperature(iLon,iLat,iAlt,iBlock)* &
                   GradLogRho(iLon,iLat,iAlt) &
                   - Gravity(iAlt)) &
                   + Dt * DiffVel(iLon,iLat,iAlt)

              NewTemp(iLon,iLat,iAlt)   = Temperature(iLon,iLat,iAlt,iBlock) &
                   - Dt * &
                   (Velocity(iLon,iLat,iAlt,iUp_,iBlock)*&
                   GradTemp(iLon,iLat,iAlt) &
                   + (Gamma(iLon,iLat,iAlt,iBlock) - 1.0) * Temperature(iLon,iLat,iAlt,iBlock)* &
                   DivVel(iLon,iLat,iAlt)) &
                   + Dt * DiffTemp(iLon,iLat,iAlt)

           enddo
        enddo
     enddo

     iLon = 3
     iLat = 3
     iAlt = 1

!     if (iBlock == 1 .and. iProc == 0) then
!        write(*,*) "newlogrho : ",NewLogRho(3,3,1), &
!             minval(velocity(3,3,:,iUp_,iBlock)), &
!             maxval(velocity(3,3,:,iUp_,iBlock)), &
!             minval(Temperature(3,3,:,iBlock)), &
!             maxval(Temperature(3,3,:,iBlock))
!     endif

     do iAlt = 1,nAlts
        NewTemp(1:nLons,1:nLats,iAlt) = NewTemp(1:nLons,1:nLats,iAlt) + &
             dt * EuvHeating(1:nLons,1:nLats,iAlt,iBlock)
     enddo

     ! Conduction (1/rho)*d(kappa dT/dZ)/dZ
     do iAlt = 1,nAlts+1
        Conduction(:,:,iAlt) = KappaTemp(:,:,iAlt,iBlock) * &
             (Temperature(1:nLons,1:nLats,iAlt,iBlock) - &
              Temperature(1:nLons,1:nLats,iAlt-1,iBlock)) / &
              dAlt(iAlt)
     end do
     do iAlt = 1,nAlts

        NewTemp(1:nLons,1:nLats,iAlt)   = NewTemp(1:nLons,1:nLats,iAlt) + &
             Dt * &
             (Conduction(1:nLons,1:nLats,iAlt+1) - &
              Conduction(1:nLons,1:nLats,iAlt)) / &
              (dAlt(iAlt)*exp(LogRho(1:nLons,1:nLats,iAlt,iBlock)))

     end do

     LogRho(1:nLons,1:nLats,1:nAlts,iBlock) = &
          NewLogRho(1:nLons,1:nLats,1:nAlts)
     Velocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) = &
          NewVel(1:nLons,1:nLats,1:nAlts)
     Temperature(1:nLons,1:nLats,1:nAlts,iBlock) = &
          NewTemp(1:nLons,1:nLats,1:nAlts)

!     write(*,*) "logrho : ", minval(NewLogRho(1:nLons,1:nLats,1:nAlts)), &
!          maxval(NewLogRho(1:nLons,1:nLats,1:nAlts))

  enddo

!  if (iProc == 0) then
!     write(*,*) "mm(velocity) : ",&
!          minval(Velocity(1:nLons,1:nLats,1:nAlts,iUp_,:)),&
!          maxval(Velocity(1:nLons,1:nLats,1:nAlts,iUp_,:))
!  endif

  tSimulation = tSimulation + dt
  CurrentTime = CurrentTime + dt

  call time_real_to_int(CurrentTime, iTimeArray)
  iJulianDay = jday(iTimeArray(1), iTimeArray(2), iTimeArray(3)) 
  utime = iTimeArray(4)*3600.0 + iTimeArray(5)*60.0 + &
       iTimeArray(6) + iTimeArray(7)/1000.0

  call end_timing("advance")

end subroutine advance
