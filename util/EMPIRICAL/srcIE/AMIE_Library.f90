!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------

subroutine AMIE_SetFileName(cFileNameIn)
  use ModAMIE_Interface
  implicit none
  character (len=iCharLenIE_), intent(in) :: cFileNameIn
  AMIE_FileName = cFileNameIn
end subroutine AMIE_SetFileName

!----------------------------------------------------------------------

subroutine AMIE_GetFileName(cFileNameOut)
  use ModAMIE_Interface
  implicit none
  character (len=iCharLenIE_), intent(out) :: cFileNameOut
  cFileNameOut = AMIE_FileName
end subroutine AMIE_GetFileName

!----------------------------------------------------------------------

subroutine AMIE_GetnTimes(nTimesOut)
  use ModAMIE_Interface
  implicit none
  integer, intent(out) :: nTimesOut
  nTimesOut = AMIE_nTimes
end subroutine AMIE_GetnTimes

!----------------------------------------------------------------------

subroutine AMIE_GetnMLTs(nMLTsOut)
  use ModAMIE_Interface
  implicit none
  integer, intent(out) :: nMLTsOut
  nMLTsOut = AMIE_nMLTs
end subroutine AMIE_GetnMLTs

!----------------------------------------------------------------------

subroutine AMIE_GetnLats(nLatsOut)
  use ModAMIE_Interface
  implicit none
  integer, intent(out) :: nLatsOut
  nLatsOut = AMIE_nLats
end subroutine AMIE_GetnLats

!----------------------------------------------------------------------

subroutine AMIE_GetLats(IEi_nMLTs, IEi_nLats, IEi_nBLKs, LatsOut)

  use ModAMIE_Interface

  implicit none

  integer, intent(in) :: IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMlts,IEi_nLats,IEi_nBLKs), intent(out) :: LatsOut
  integer :: i,j

  do i=1,IEi_nMLTs
     do j=IEi_nLats,1,-1
        LatsOut(i,j,AMIE_North_) = AMIE_Lats(IEi_nLats-j+1)
     enddo
     LatsOut(i,1:IEi_nLats,AMIE_South_) = -AMIE_Lats(1:IEi_nLats)
  enddo

end subroutine AMIE_GetLats

!----------------------------------------------------------------------

subroutine AMIE_GetMLTs(IEi_nMLTs, IEi_nLats, IEi_nBLKs, MLTsOut)

  use ModAMIE_Interface

  implicit none

  integer, intent(in) :: IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMlts,IEi_nLats,IEi_nBLKs), intent(out) :: MLTsOut
  integer :: j

  do j=1,IEi_nLats
     MLTsOut(1:IEi_nMLTs,j,1) = AMIE_MLTs(1:IEi_nMLTs)
     MLTsOut(1:IEi_nMLTs,j,2) = AMIE_MLTs(1:IEi_nMLTs)
  enddo

end subroutine AMIE_GetMLTs

!----------------------------------------------------------------------
! If there is an error, report it and crash
!----------------------------------------------------------------------

subroutine report_error_and_crash(iError, StringSource)

  use ModErrors

  implicit none

  integer, intent(in) :: iError
  character (len=*) :: StringSource

  if (iError /= 0) then
     write(*,*) "Error in AMIE_Library!"
     write(*,*) "Source : ", StringSource
     write(*,*) cErrorCodes(iError)
     stop
  endif

end subroutine report_error_and_crash

!----------------------------------------------------------------------
! New method!
! East/West potential(y)
!----------------------------------------------------------------------

subroutine AMIE_GetPotentialY_New(TimeIn, Method, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, PotentialYOut, iError)

  use ModAMIE_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: PotentialYOut
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  AMIE_Value = AMIE_PotentialY
  
  call AMIE_GetValue_New(TimeIn, Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs, &
       ValueOut, iError)

  call report_error_and_crash(iError, "AMIE_GetPotentialY_New")

  PotentialYOut = ValueOut

end subroutine AMIE_GetPotentialY_New

!----------------------------------------------------------------------
! New method!
! Overall potential
!----------------------------------------------------------------------

subroutine AMIE_GetPotential_New(TimeIn, Method, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, PotentialOut, iError)

  use ModAMIE_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: PotentialOut
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  AMIE_Value = AMIE_potential
  
  call AMIE_GetValue_New(TimeIn, Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs, &
       ValueOut, iError)

  call report_error_and_crash(iError, "AMIE_GetPotential_New")

  PotentialOut = ValueOut

end subroutine AMIE_GetPotential_New

!----------------------------------------------------------------------
! New method!
! Electron Energy Flux
!----------------------------------------------------------------------

subroutine AMIE_GetEFlux_New(TimeIn, Method, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, EFluxOut, iError)

  use ModAMIE_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: EFluxOut
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  AMIE_Value = AMIE_eflux
  
  call AMIE_GetValue_New(TimeIn, Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs, &
       ValueOut, iError)

  call report_error_and_crash(iError, "AMIE_GetEFlux_New")

  EFluxOut = ValueOut

end subroutine AMIE_GetEFlux_New

!----------------------------------------------------------------------
! New method!
! Electron Average Energy
!----------------------------------------------------------------------

subroutine AMIE_GetAveE_New(TimeIn, Method, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, AveEOut, iError)

  use ModAMIE_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: AveEOut
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  AMIE_Value = AMIE_avee
  
  call AMIE_GetValue_New(TimeIn, Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs, &
       ValueOut, iError)

  call report_error_and_crash(iError, "AMIE_GetAveE_New")

  AveEOut = ValueOut

end subroutine AMIE_GetAveE_New

!----------------------------------------------------------------------
! New method!
! Ion Energy Flux
!----------------------------------------------------------------------

subroutine AMIE_GetIonEFlux_New(TimeIn, Method, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, IonEFluxOut, iError)

  use ModAMIE_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: IonEFluxOut
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  AMIE_Value = AMIE_IonEflux
  
  call AMIE_GetValue_New(TimeIn, Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs, &
       ValueOut, iError)

  call report_error_and_crash(iError, "AMIE_GetIonEFlux_New")

  IonEFluxOut = ValueOut

end subroutine AMIE_GetIonEFlux_New

!----------------------------------------------------------------------
! New method!
! Ion Average Energy
!----------------------------------------------------------------------

subroutine AMIE_GetIonAveE_New(TimeIn, Method, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, IonAveEOut, iError)

  use ModAMIE_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: IonAveEOut
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  AMIE_Value = AMIE_IonAveE

  call AMIE_GetValue_New(TimeIn, Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs, &
       ValueOut, iError)

  call report_error_and_crash(iError, "AMIE_GetIonAveE_New")

  IonAveEOut = ValueOut

end subroutine AMIE_GetIonAveE_New

!----------------------------------------------------------------------
! "East/West" Potential
!----------------------------------------------------------------------

subroutine AMIE_GetPotentialY(TimeIn, Method, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, PotentialOut, iError)

  use ModAMIE_Interface
  use ModErrors

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: PotentialOut
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  call AMIE_GetValue(TimeIn, Method, potentialy_, &
       IEi_nMLTs, IEi_nLats, IEi_nBLKs, ValueOut, iError)

  if (iError /= 0) then
     write(*,*) "Error in routine AMIE_GetPotential:"
     write(*,*) cErrorCodes(iError)
     stop
  else
     PotentialOut = ValueOut
!     write(*,*) 'AMIE_GetPotentialY PotentialOut = ', PotentialOut
  endif

end subroutine AMIE_GetPotentialY

!----------------------------------------------------------------------
! Get General Values
!----------------------------------------------------------------------

subroutine AMIE_GetValue(TimeIn, Method, iValue, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, ValueOut, iError)

  use ModErrors
  use ModAMIE_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs, iValue
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: ValueOut
  integer, intent(out) :: iError

  integer :: iTime, i, j, iLat, iBLK
  logical :: IsDone
  real*8  :: dT, VerySmall = 1.0e-6


  iError = 0

  do iBLK = AMIE_South_, AMIE_North_

     IsDone = .false.
     iTime = 1

     do while (.not. IsDone)
        if (TimeIn - AMIE_Time(iTime,iBLK) < VerySmall) IsDone = .true.
        if ((iTime == AMIE_nTimes) .and. (.not.IsDone)) then
           iTime = iTime + 1 
           IsDone = .true.
        endif
        iTime = iTime + 1
     enddo

     if (iTime <= AMIE_nTimes+1) then

        iTime = iTime - 1

        if (iTime == 1) then

           ! If we are before the start time, allow users to extrapolate
           ! up to 5 dT.

           dT = AMIE_Time(2,iBLK) - AMIE_Time(1,iBLK)
           if (TimeIn + 5*dt < AMIE_Time(1,iBLK)) then
              ValueOut = -1.0e32
              iError = ecBeforeStartTime_
              return
           endif
        endif
     else
        dT = AMIE_Time(2,iBLK) - AMIE_Time(1,iBLK)

        ! If we are after the end time, allow users to extrapolate
        ! up to 5 dT.

        if (TimeIn - 5*dt < AMIE_Time(AMIE_nTimes,iBLK)) then
           iTime = AMIE_nTimes
        else
           ValueOut = -1.0e32
           iError = ecAfterEndTime_
           return
        endif
     endif

     if (Method == AMIE_After_) then
        if (iBLK == AMIE_South_) then

           if (iValue == potential_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   AMIE_Potential(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           endif

           if (iValue == potentialy_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   AMIE_PotentialY(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           endif

           if (iValue == eflux_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   AMIE_EFlux(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           endif

           if (iValue == avee_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   AMIE_AveE(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           endif

           ! Reverse the North block of AMIE data for now...
        else
           do iLat = AMIE_nLats,1,-1

              if (iValue == potential_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      AMIE_Potential(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
              endif

              if (iValue == potentialy_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      AMIE_PotentialY(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
              endif

              if (iValue == eflux_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      AMIE_EFlux(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
              endif

              if (iValue == avee_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      AMIE_AveE(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
              endif

           enddo
        endif
     endif

     if (Method == AMIE_Closest_) then
        if (iTime > 1) then
           if (abs(TimeIn-AMIE_Time(iTime,iBLK)) > &
               abs(TimeIn-AMIE_Time(iTime-1,iBLK))) &
               iTime = iTime - 1
        endif
        if (iBLK == AMIE_South_) then
           if (iValue == potential_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   AMIE_Potential(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           endif
           if (iValue == potentialy_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   AMIE_PotentialY(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           endif
           if (iValue == eflux_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   AMIE_EFlux(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           endif
           if (iValue == avee_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   AMIE_AveE(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           endif
        else
           ! Reverse the North block of AMIE data for now...
           do iLat = AMIE_nLats,1,-1
              if (iValue == potential_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      AMIE_Potential(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
              endif
              if (iValue == potentialy_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      AMIE_PotentialY(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
              endif
              if (iValue == eflux_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      AMIE_EFlux(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
              endif
              if (iValue == avee_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      AMIE_AveE(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
              endif
           enddo
        endif
     endif

     if (Method == AMIE_Interpolate_) then
        ! This will do extrapolation if it is before the first time
        if (iTime == 1) iTime = iTime + 1
        ! dT is the percentage of the way away from the current point
        dT = (AMIE_Time(iTime,iBLK) - TimeIn) / &
             (AMIE_Time(iTime,iBLK) - AMIE_Time(iTime-1,iBLK))

        ! Use 1-dT for the selected point, since dt = 0 if you are exactly
        ! on the selected point
        if (iBLK == AMIE_South_) then
           if (iValue == potential_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   (1.0 - dt)*AMIE_Potential(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)+&
                           dt*AMIE_Potential(1:AMIE_nMLTs, 1:AMIE_nLats,iTime-1,iBLK)
           endif
           if (iValue == potentialy_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   (1.0 - dt)*AMIE_PotentialY(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)+&
                           dt*AMIE_PotentialY(1:AMIE_nMLTs, 1:AMIE_nLats,iTime-1,iBLK)
           endif
           if (iValue == eflux_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   (1.0 - dt)*AMIE_EFlux(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)+&
                           dt*AMIE_EFlux(1:AMIE_nMLTs, 1:AMIE_nLats,iTime-1,iBLK)
           endif
           if (iValue == avee_) then
              ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                   (1.0 - dt)*AMIE_AveE(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)+&
                           dt*AMIE_AveE(1:AMIE_nMLTs, 1:AMIE_nLats,iTime-1,iBLK)
           endif
        else
           ! Reverse the 2nd block of AMIE data for now...
           do iLat = AMIE_nLats,1,-1
              if (iValue == potential_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      (1.0 - dt)*AMIE_Potential(1:AMIE_nMLTs,AMIE_nLats-iLat+1,&
                                         iTime,iBLK) + &
                     dt*AMIE_Potential(1:AMIE_nMLTs, AMIE_nLats-iLat+1,iTime-1,iBLK)
              endif
              if (iValue == potentialy_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      (1.0 - dt)*AMIE_PotentialY(1:AMIE_nMLTs,AMIE_nLats-iLat+1,&
                                         iTime,iBLK) + &
                     dt*AMIE_PotentialY(1:AMIE_nMLTs, AMIE_nLats-iLat+1,iTime-1,iBLK)
              endif
              if (iValue == eflux_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      (1.0 - dt)*AMIE_EFlux(1:AMIE_nMLTs,AMIE_nLats-iLat+1,&
                                         iTime,iBLK) + &
                     dt*AMIE_EFlux(1:AMIE_nMLTs, AMIE_nLats-iLat+1,iTime-1,iBLK)
              endif
              if (iValue == avee_) then
                 ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                      (1.0 - dt)*AMIE_AveE(1:AMIE_nMLTs,AMIE_nLats-iLat+1,&
                                         iTime,iBLK) + &
                     dt*AMIE_AveE(1:AMIE_nMLTs, AMIE_nLats-iLat+1,iTime-1,iBLK)
              endif
           enddo
        endif
     endif

  enddo

end subroutine AMIE_GetValue

!----------------------------------------------------------------------
! New Method!!!
! Get General Values
!----------------------------------------------------------------------

subroutine AMIE_GetValue_New(TimeIn, Method, &
     IEi_nMLTs, IEi_nLats, IEi_nBLKs, ValueOut, iError)

  use ModErrors
  use ModAMIE_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
  real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: ValueOut
  integer, intent(out) :: iError

  integer :: iTime, i, j, iLat, iBLK, iL
  logical :: IsDone
  real*8  :: dT, VerySmall = 1.0e-6

  iError = 0

  do iBLK = AMIE_South_, AMIE_North_

     IsDone = .false.
     iTime = 1

     do while (.not. IsDone)
        if (TimeIn - AMIE_Time(iTime,iBLK) < VerySmall) IsDone = .true.
        if ((iTime == AMIE_nTimes) .and. (.not.IsDone)) then
           iTime = iTime + 1 
           IsDone = .true.
        endif
        iTime = iTime + 1
     enddo

     if (iTime <= AMIE_nTimes+1) then

        iTime = iTime - 1

        if (iTime == 1) then

           ! If we are before the start time, allow users to extrapolate
           ! up to 5 dT.

           dT = AMIE_Time(2,iBLK) - AMIE_Time(1,iBLK)
           if (TimeIn + 5*dt < AMIE_Time(1,iBLK)) then
              ValueOut = -1.0e32
              iError = ecBeforeStartTime_
              return
           endif
        endif
     else
        dT = AMIE_Time(2,iBLK) - AMIE_Time(1,iBLK)

        ! If we are after the end time, allow users to extrapolate
        ! up to 5 dT.

        if (TimeIn - 5*dt < AMIE_Time(AMIE_nTimes,iBLK)) then
           iTime = AMIE_nTimes
        else
           ValueOut = -1.0e32
           iError = ecAfterEndTime_
           return
        endif
     endif

     if ( Method == AMIE_After_ .or. &
          Method == AMIE_Closest_) then

        if (Method == AMIE_Closest_) then
           if (iTime > 1) then
              if (abs(TimeIn-AMIE_Time(iTime,iBLK)) > &
                   abs(TimeIn-AMIE_Time(iTime-1,iBLK))) &
                   iTime = iTime - 1
           endif
        endif

        if (iBLK == AMIE_South_) then
           ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                AMIE_Value(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
        else
           ! Reverse the North block of AMIE data for now...
           do iLat = AMIE_nLats,1,-1
              ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                   AMIE_Value(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
           enddo
        endif

     endif

     if (Method == AMIE_Interpolate_) then
        ! This will do extrapolation if it is before the first time
        if (iTime == 1) iTime = iTime + 1
        ! dT is the percentage of the way away from the current point
        dT = (AMIE_Time(iTime,iBLK) - TimeIn) / &
             (AMIE_Time(iTime,iBLK) - AMIE_Time(iTime-1,iBLK))

        ! Use 1-dT for the selected point, since dt = 0 if you are exactly
        ! on the selected point
        if (iBLK == AMIE_South_) then
           ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                (1.0 - dt)*AMIE_Value(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)+&
                        dt*AMIE_Value(1:AMIE_nMLTs, 1:AMIE_nLats,iTime-1,iBLK)
        else
           ! Reverse the 2nd block of AMIE data for now...
           do iLat = AMIE_nLats,1,-1
              iL = AMIE_nLats-iLat+1
              ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  &
                   (1.0 - dt)*AMIE_Value(1:AMIE_nMLTs, iL, iTime, iBLK) + &
                           dt*AMIE_Value(1:AMIE_nMLTs, iL, iTime-1, iBLK)
           enddo
        endif
     endif

  enddo

end subroutine AMIE_GetValue_New

!----------------------------------------------------------------------
! Get All Values:
! - potential
! - e- average energy
! - e- total energy flux
! - ion average energy
! - ion total energy flux
!----------------------------------------------------------------------

subroutine get_AMIE_values(rtime)

  use ModEIE_Interface

  real*8, intent(in) :: rtime
  integer :: iError

  call AMIE_GetPotential_New(rtime, EIE_Interpolate_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, &
       EIEr3_HavePotential, iError)

  call AMIE_GetAveE_New(rtime, EIE_Closest_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HaveAveE, iError)

  call AMIE_GetEFlux_New(rtime, EIE_Closest_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HaveEFlux, iError)

  call AMIE_GetIonAveE_New(rtime, EIE_Closest_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HaveIonAveE, iError)

  call AMIE_GetIonEFlux_New(rtime, EIE_Closest_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HaveIonEFlux, iError)

end subroutine get_AMIE_values

!----------------------------------------------------------------------
! Get the secondary potential (this is the potential for the E/W direction)
!----------------------------------------------------------------------

subroutine get_AMIE_PotentialY(rtime)

  use ModEIE_Interface

  real*8, intent(in) :: rtime
  integer :: iError

  call AMIE_GetPotentialY(rtime, EIE_Interpolate_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, &
       EIEr3_HavePotential, iError)

end subroutine get_AMIE_PotentialY
