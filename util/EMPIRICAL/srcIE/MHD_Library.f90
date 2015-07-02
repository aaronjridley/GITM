!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------

subroutine MHD_SetFileName(cFileNameIn)
  use ModMHD_Interface
  implicit none
  character (len=100), intent(in) :: cFileNameIn
  MHD_FileName = cFileNameIn
end subroutine MHD_SetFileName

!----------------------------------------------------------------------

subroutine MHD_SetYear(iYearIn)
  use ModMHD_Interface
  implicit none
  integer, intent(in) :: iYearIn
  MHD_Year = iYearIn
end subroutine MHD_SetYear

!----------------------------------------------------------------------

subroutine MHD_SetMonth(iMonthIn)
  use ModMHD_Interface
  implicit none
  integer, intent(in) :: iMonthIn
  MHD_Month = iMonthIn
end subroutine MHD_SetMonth

!----------------------------------------------------------------------

subroutine MHD_SetDay(iDayIn)
  use ModMHD_Interface
  implicit none
  integer, intent(in) :: iDayIn
  MHD_Day = iDayIn
end subroutine MHD_SetDay

!----------------------------------------------------------------------

subroutine MHD_GetFileName(cFileNameOut)
  use ModMHD_Interface
  implicit none
  character (len=100), intent(out) :: cFileNameOut
  cFileNameOut = MHD_FileName
end subroutine MHD_GetFileName

!----------------------------------------------------------------------

subroutine MHD_GetnTimes(nTimesOut)
  use ModMHD_Interface
  implicit none
  integer, intent(out) :: nTimesOut
  nTimesOut = MHD_nTimes
end subroutine MHD_GetnTimes

!----------------------------------------------------------------------

subroutine MHD_GetnMLTs(nMLTsOut)
  use ModMHD_Interface
  implicit none
  integer, intent(out) :: nMLTsOut
  nMLTsOut = MHD_nMLTs
end subroutine MHD_GetnMLTs

!----------------------------------------------------------------------

subroutine MHD_GetnLats(nLatsOut)
  use ModMHD_Interface
  implicit none
  integer, intent(out) :: nLatsOut
  nLatsOut = MHD_nLats
end subroutine MHD_GetnLats

!----------------------------------------------------------------------

subroutine MHD_GetLats(EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, LatsOut)

  use ModMHD_Interface

  implicit none

  integer, intent(in) :: EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs
  real, dimension(EIEi_nMlts,EIEi_nLats,EIEi_nBLKs), intent(out) :: LatsOut
  integer :: i,j

  do i=1,EIEi_nMLTs
     do j=EIEi_nLats,1,-1
        LatsOut(i,j,MHD_North_) = MHD_Lats(EIEi_nLats-j+1)
     enddo
     LatsOut(i,1:EIEi_nLats,MHD_South_) = -MHD_Lats(1:EIEi_nLats)
  enddo

end subroutine MHD_GetLats

!----------------------------------------------------------------------

subroutine MHD_GetMLTs(EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, MLTsOut)

  use ModMHD_Interface

  implicit none

  integer, intent(in) :: EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs
  real, dimension(EIEi_nMlts,EIEi_nLats,EIEi_nBLKs), intent(out) :: MLTsOut
  integer :: j

  do j=1,EIEi_nLats
     MLTsOut(1:EIEi_nMLTs,j,1) = MHD_MLTs(1:EIEi_nMLTs)
     MLTsOut(1:EIEi_nMLTs,j,2) = MHD_MLTs(1:EIEi_nMLTs)
  enddo

end subroutine MHD_GetMLTs

!----------------------------------------------------------------------

subroutine MHD_GetPotential(TimeIn, Method, &
     EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, PotentialOut, iError)

  use ModMHD_Interface
  use ModErrors

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs
  real, dimension(EIEi_nMLTs,EIEi_nLats,EIEi_nBLKs), intent(out) :: PotentialOut
  real, dimension(EIEi_nMLTs,EIEi_nLats,EIEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  MHD_Value = MHD_Potential

  call MHD_GetValue(TimeIn, Method, &
       EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, ValueOut, iError)

  if (iError /= 0) then
     write(*,*) "Error in routine MHD_GetPotential:"
     write(*,*) cErrorCodes(iError)
     stop
  else
     PotentialOut = ValueOut
  endif

end subroutine MHD_GetPotential

!----------------------------------------------------------------------

subroutine MHD_GetEFlux(TimeIn, Method, &
     EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, EFluxOut, iError)

  use ModMHD_Interface
  use ModErrors

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs
  real, dimension(EIEi_nMLTs,EIEi_nLats,EIEi_nBLKs), intent(out) :: EFluxOut
  real, dimension(EIEi_nMLTs,EIEi_nLats,EIEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  MHD_Value = MHD_EFlux

  call MHD_GetValue(TimeIn, Method, &
       EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, ValueOut, iError)

  if (iError /= 0) then
     write(*,*) "Error in routine MHD_GetEFlux:"
     write(*,*) cErrorCodes(iError)
     stop
  else
     EFluxOut = ValueOut
  endif

end subroutine MHD_GetEFlux

!----------------------------------------------------------------------

subroutine MHD_GetAveE(TimeIn, Method, &
     EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, AveEOut, iError)

  use ModMHD_Interface
  use ModErrors

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs
  real, dimension(EIEi_nMLTs,EIEi_nLats, EIEi_nBLKs), intent(out) :: AveEOut
  real, dimension(EIEi_nMLTs,EIEi_nLats, EIEi_nBLKs)              :: ValueOut
  integer, intent(out) :: iError

  MHD_Value = MHD_AveE
  call MHD_GetValue(TimeIn, Method, &
       EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, ValueOut, iError)

  if (iError /= 0) then
     write(*,*) "Error in routine MHD_GetAveE:"
     write(*,*) cErrorCodes(iError)
     stop
  else
     AveEOut = ValueOut
  endif

end subroutine MHD_GetAveE

!----------------------------------------------------------------------

subroutine MHD_GetValue(TimeIn, Method, &
     EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs, ValueOut, iError)

  use ModErrors
  use ModMHD_Interface

  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method, EIEi_nMLTs, EIEi_nLats, EIEi_nBLKs
  real, dimension(EIEi_nMLTs,EIEi_nLats,EIEi_nBLKs), intent(out) :: ValueOut
  integer, intent(out) :: iError

  integer :: iTime, i, j, iLat, iBLK
  logical :: IsDone
  real*8  :: dT, VerySmall = 1.0e-6

  iError = 0

  do iBLK = MHD_South_, MHD_North_

     IsDone = .false.
     iTime = 1

     do while (.not. IsDone)
        if (TimeIn - MHD_Time(iTime,iBLK) < VerySmall) IsDone = .true.
        if ((iTime == MHD_nTimes) .and. (.not.IsDone)) then
           iTime = iTime + 1 
           IsDone = .true.
        endif
        iTime = iTime + 1
     enddo

     if (iTime <= MHD_nTimes+1) then

        iTime = iTime - 1

        if (iTime == 1) then

           ! If we are before the start time, allow users to extrapolate
           ! up to 1 dT.

           dT = MHD_Time(2,iBLK) - MHD_Time(1,iBLK)
           if (TimeIn + dt < MHD_Time(1,iBLK)) then
              ValueOut = -1.0e32
              iError = ecBeforeStartTime_
              return
           endif
        endif
     else
        dT = MHD_Time(2,iBLK) - MHD_Time(1,iBLK)

        ! If we are after the end time, allow users to extrapolate
        ! up to 1 dT.

        if (TimeIn - dt < MHD_Time(MHD_nTimes,iBLK)) then
           iTime = MHD_nTimes
        else
           ValueOut = -1.0e32
           iError = ecAfterEndTime_
           return
        endif
     endif

     if (Method == MHD_After_) then
        if (iBLK == MHD_South_) then
           ValueOut(1:MHD_nMLTs, 1:MHD_nLats,iBLK) =  &
                MHD_Value(1:MHD_nMLTs, 1:MHD_nLats,iTime,iBLK)
           ! Reverse the North block of MHD data for now...
        else
           do iLat = MHD_nLats,1,-1
              ValueOut(1:MHD_nMLTs, iLat,iBLK) =  &
                   MHD_Value(1:MHD_nMLTs, MHD_nLats - iLat + 1,iTime,iBLK)
           enddo
        endif
     endif

     if (Method == MHD_Closest_) then
        if (iTime > 1) then
           if (abs(TimeIn-MHD_Time(iTime,iBLK)) > &
               abs(TimeIn-MHD_Time(iTime-1,iBLK))) &
               iTime = iTime - 1
        endif
        if (iBLK == MHD_South_) then
           ValueOut(1:MHD_nMLTs, 1:MHD_nLats,iBLK) =  &
                MHD_Value(1:MHD_nMLTs, 1:MHD_nLats,iTime,iBLK)
        else
           ! Reverse the North block of MHD data for now...
           do iLat = MHD_nLats,1,-1
              ValueOut(1:MHD_nMLTs, iLat,iBLK) =  &
                   MHD_Value(1:MHD_nMLTs, MHD_nLats - iLat + 1,iTime,iBLK)
           enddo
        endif
     endif

     if (Method == MHD_Interpolate_) then
        ! This will do extrapolation if it is before the first time
        if (iTime == 1) iTime = iTime + 1
        ! dT is the percentage of the way away from the current point
        dT = (MHD_Time(iTime,iBLK) - TimeIn) / &
             (MHD_Time(iTime,iBLK) - MHD_Time(iTime-1,iBLK))
        ! Use 1-dT for the selected point, since dt = 0 if you are exactly
        ! on the selected point
        if (iBLK == MHD_South_) then
           ValueOut(1:MHD_nMLTs, 1:MHD_nLats,iBLK) =  &
                (1.0 - dt)*MHD_Value(1:MHD_nMLTs, 1:MHD_nLats,iTime,iBLK)+&
                       dt*MHD_Value(1:MHD_nMLTs, 1:MHD_nLats,iTime-1,iBLK)
        else
           ! Reverse the 2nd block of MHD data for now...
           do iLat = MHD_nLats,1,-1
              ValueOut(1:MHD_nMLTs, iLat,iBLK) =  &
                   (1.0 - dt)*MHD_Value(1:MHD_nMLTs,MHD_nLats-iLat+1,&
                                         iTime,iBLK) + &
                   dt*MHD_Value(1:MHD_nMLTs, MHD_nLats-iLat+1,iTime-1,iBLK)
           enddo
        endif
     endif

  enddo

end subroutine MHD_GetValue

subroutine get_MHD_values(rtime)

  use ModEIE_Interface

  real*8, intent(in) :: rtime
  integer :: iError

  call MHD_GetPotential(rtime, EIE_Interpolate_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HavePotential, iError)

  call MHD_GetAveE(rtime, EIE_Closest_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HaveAveE, iError)

  call MHD_GetEFlux(rtime, EIE_Closest_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HaveEFlux, iError)

end subroutine get_MHD_values
