!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!----------------------------------------------------------------------

subroutine IO_SetnMLTs(nMLTsIn)
  use ModEIE_Interface
  implicit none
  integer, intent(in) :: nMLTsIn
  IOi_NeednMLTs = nMLTsIn
end subroutine IO_SetnMLTs

!----------------------------------------------------------------------

subroutine IO_SetnLats(nLatsIn)
  use ModEIE_Interface
  implicit none
  integer, intent(in) :: nLatsIn
  IOi_NeednLats = nLatsIn
end subroutine IO_SetnLats

!----------------------------------------------------------------------

subroutine IO_SetTime(TimeIn)
!  use ModKind
  use ModEIE_Interface
  implicit none
!  Integer, parameter :: Real8_ = selected_real_kind(14,200)
  real (kind=Real8_), intent(in) :: TimeIn
  IOd_NeedTime = TimeIn
end subroutine IO_SetTime

!----------------------------------------------------------------------

subroutine IO_SetIMFBz(IMFBzIn)
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: IMFBzIn
  IOr_NeedIMFBz = IMFBzIn
end subroutine IO_SetIMFBz

!----------------------------------------------------------------------

subroutine IO_SetIMFBy(IMFByIn)
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: IMFByIn
  IOr_NeedIMFBy = IMFByIn
end subroutine IO_SetIMFBy

!----------------------------------------------------------------------

subroutine IO_SetSWV(SWVIn)
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: SWVIn
  IOr_NeedSWV = SWVIn
end subroutine IO_SetSWV

!----------------------------------------------------------------------

subroutine IO_SetSWN(SWNIn)
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: SWNIn
  IOr_NeedSWN = SWNIn
end subroutine IO_SetSWN

!----------------------------------------------------------------------

subroutine IO_SetHPI(HPIIn)
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: HPIIn
  IOr_NeedHPI = HPIIn
  IOr_NeedHPINorm = 2.09 * ALOG(HPIIn) * 1.0475
end subroutine IO_SetHPI

!----------------------------------------------------------------------

subroutine IO_SetHPINorm(HPINormIn)
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: HPINormIn
  IOr_NeedHPINorm = HPINormIn
end subroutine IO_SetHPINorm

!----------------------------------------------------------------------

subroutine IO_SetKp(KpIn)
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: KpIn
  IOr_NeedKp = KpIn
end subroutine IO_SetKp

!----------------------------------------------------------------------

subroutine IO_SetNorth
!  use ModKind
  use ModEIE_Interface
  implicit none
  IOl_IsNorth = .true.
end subroutine IO_SetNorth

!----------------------------------------------------------------------

subroutine IO_SetSouth
!  use ModKind
  use ModEIE_Interface
  implicit none
  IOl_IsNorth = .false.
end subroutine IO_SetSouth

!----------------------------------------------------------------------

subroutine IO_SetMLTs(MLTsIn, iError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(in) :: MLTsIn

  integer :: i,j

  iError = 0
  if (allocated(IOr2_NeedMLTs)) deallocate(IOr2_NeedMLTs)
  allocate(IOr2_NeedMLTs(IOi_NeednMLTs,IOi_NeednLats), stat=iError)
  if (iError /= 0) then
     iError = ecAllocationError_
     return
  else
     do i=1,IOi_NeednMLTs
        do j=1,IOi_NeednLats
           IOr2_NeedMLTs(i,j) = mod((MLTsIn(i,j)+24.0),24.0)
        enddo
     enddo
  endif

end subroutine IO_SetMLTs

!----------------------------------------------------------------------

subroutine IO_SetLats(LatsIn, iError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(in) :: LatsIn

  iError = 0
  if (allocated(IOr2_NeedLats)) deallocate(IOr2_NeedLats)
  allocate(IOr2_NeedLats(IOi_NeednMLTs,IOi_NeednLats), stat=iError)
  if (iError /= 0) then
     iError = ecAllocationError_
     return
  else
     IOr2_NeedLats(1:IOi_NeednMLTs,1:IOi_NeednLats) = &
          LatsIn(1:IOi_NeednMLTs,1:IOi_NeednLats)
  endif

end subroutine IO_SetLats

!----------------------------------------------------------------------

subroutine IO_SetGrid(MLTsIn, LatsIn, iError)
  
  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(in) :: MLTsIn, LatsIn
  real, dimension(2) :: IOr1_Location
  real, dimension(5) :: EIEr1_Location

  integer :: i,j

  iError = 0

  call IO_SetMLTs(MLTsIn, iError)
  if (iError /= 0) return

  call IO_SetLats(LatsIn, iError)
  if (iError /= 0) return

  if (UseGridBasedEIE) then

     if (allocated(IOi3_InterpolationIndices)) &
          deallocate(IOi3_InterpolationIndices)
     allocate(IOi3_InterpolationIndices(IOi_NeednMLTs,IOi_NeednLats,3), &
          stat=iError)
     if (iError /= 0) then
        iError = ecAllocationError_
        return
     endif

     if (allocated(IOr3_InterpolationRatios)) &
          deallocate(IOr3_InterpolationRatios)
     allocate(IOr3_InterpolationRatios(IOi_NeednMLTs,IOi_NeednLats,2), &
          stat=iError)
     if (iError /= 0) then
        iError = ecAllocationError_
        return
     endif

     do i=1,IOi_NeednLats
        do j=1,IOi_NeednMLTs

           IOr1_Location(1) = mod((IOr2_NeedMLTs(j,i) + 24.0),24.0)
           IOr1_Location(2) = IOr2_NeedLats(j,i)

           call EIE_FindPoint(IOr1_Location, EIEr1_Location, iError)

           if (iError == 0) then
              IOi3_InterpolationIndices(j,i,1:3) = EIEr1_Location(1:3)
              IOr3_InterpolationRatios(j,i,1:2)  = EIEr1_Location(4:5)
           else
              IOi3_InterpolationIndices(j,i,1) = -1
           endif
           
        enddo
     enddo

  endif

end subroutine IO_SetGrid

!----------------------------------------------------------------------

subroutine IO_GetPotential(PotentialOut, iError)
  
  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs) :: EIE_Value
  real, dimension(IOi_NeednMLTs,IOi_NeednLats)               :: ValueOut
  real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: PotentialOut
  real :: Filler = 0.0

  iError = 0

  if (UseGridBasedEIE) then

     EIE_Value = EIEr3_HavePotential

     call IO_GetValue(EIE_Value, ValueOut, Filler, iError)

     if (iError /= 0) then
        write(*,*) "Error in routine IO_GetPotential:"
        write(*,*) cErrorCodes(iError)
        stop
     else
        PotentialOut = ValueOut
     endif

  else

     call IO_GetNonGridBasedPotential(ValueOut, iError)

     if (iError == 0) then
        PotentialOut = ValueOut
     else
        PotentialOut = -1.0e32
     endif

  endif

end subroutine IO_GetPotential

subroutine IO_GetNonGridBasedPotential(PotentialOut, iError)

  use EIE_ModWeimer, only: WeiEpot96, get_tilt
  use w05sc, only: setmodel, epotval

  use ModEIE_Interface
  use ModErrors
  use ModTimeConvert, ONLY: time_real_to_int

  implicit none

  integer, intent(out)                                      :: iError
  real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out) :: PotentialOut

  integer, dimension(7) :: itime
  integer :: iYear, iMonth, iDay, iHour, iMinute, iSecond

  integer :: iChange, iHemisphere

  integer :: iMLT, iLat

  real    :: Potential, ETheta, EPhi, MLT, Lat, tilt, hr

  iError = 0

  PotentialOut = 0.0

  if (index(EIE_NameOfEFieldModel,'zero') > 0) then
     return
  endif

  if (IOd_NeedTime < 0.0) then
     iError = ecTimeNotSet_
     return
  endif

  if (index(EIE_NameOfEFieldModel,'hpi') <= 0) then
     if (IOr_NeedIMFBz < -1000.0) then
        iError = ecIMFBzNotSet_
        return
     endif
     if (IOr_NeedIMFBy < -1000.0) then
        iError = ecIMFByNotSet_
        return
     endif
  endif

  if (index(EIE_NameOfEFieldModel,'weimer05') > 0) then
     if (IOr_NeedSWV < -1000.0) then
        iError = ecSWVNotSet_
        return
     endif
     if (IOr_NeedSWN < -1000.0) then
        iError = ecSWNNotSet_
        return
     endif
     iChange = 1
  endif

  if (index(EIE_NameOfEFieldModel,'weimer01') > 0) then
     if (IOr_NeedSWV < -1000.0) then
        iError = ecSWVNotSet_
        return
     endif
     iChange = 1
  endif

  if (index(EIE_NameOfEFieldModel,'weimer96') > 0) then
     if (IOr_NeedSWV < -1000.0) then
        iError = ecSWVNotSet_
        return
     endif
     iChange = 1
  endif

  if (index(EIE_NameOfEFieldModel,'millstone') > 0) then
     if (IOr_NeedHPI < -1000.0) then
        iError = ecHPINotSet_
        return
     endif
  endif

  if (index(EIE_NameOfEFieldModel,'hmr89') > 0) then
     if (IOr_NeedKP < -1000.0) then
        iError = ecKPNotSet_
        return
     endif
     iChange = 1
  endif

  call time_real_to_int(IOd_NeedTime, itime)

  iYear   = itime(1)
  iMonth  = itime(2)
  iDay    = itime(3)
  iHour   = itime(4)
  iMinute = itime(5)
  iSecond = itime(6)

  Lat = 0.0

  do iMLT = 1, IOi_NeednMLTs
     do iLat = 1, IOi_NeednLats

        if (iLat == 1) then 
           iChange = 1 
        else
           if (IOr2_NeedLats(iMLT,iLat)*IOr2_NeedLats(iMLT,iLat-1) <= 0.0) &
                iChange = 1
        endif

        MLT = IOr2_NeedMLTs(iMLT,iLat)
        Lat = IOr2_NeedLats(iMLT,iLat)

        if (IOl_IsNorth) then
           if (Lat > 0.0) then
              iHemisphere = 1
           else
              Lat = abs(Lat)
              iHemisphere = -1
           endif
        else
           iHemisphere = -1
           if (Lat < 0.0) Lat = abs(Lat)
        endif

        if (index(EIE_NameOfEFieldModel,'weimer05') > 0) then
           if (abs(lat) >= 45.0) then 
              hr = float(iHour) + float(iMinute)/60.0
              tilt = get_tilt (iYear,iMonth,iDay,hr)
              call setmodel(IOr_NeedIMFBy,IOr_NeedIMFBz,tilt, &
                   IOr_NeedSWV,IOr_NeedSWN,'epot')
              call epotval(lat,mlt,0.0,Potential)
              Potential = Potential * 1000.0
              iChange = 0
           else
              Potential = 0.0
              ETheta = 0.0
              EPhi = 0.0
           endif
	endif

        if (index(EIE_NameOfEFieldModel,'weimer96') > 0) then
           if (abs(lat) >= 45.0) then 
              !call timing_start("weimer")
              call WEIEPOT96(iYear, iMonth, iDay, iHour, iMinute, IOr_NeedSWV, &
                   IOr_NeedIMFBy, IOr_NeedIMFBz, iHemisphere, iChange, &
                   Lat, MLT, ETheta, EPhi, Potential)
              !call timing_stop("weimer")
              Potential = Potential * 1000.0
              iChange = 0
           else
              Potential = 0.0
              ETheta = 0.0
              EPhi = 0.0
           endif
        endif

        if (index(EIE_NameOfEFieldModel,'millstone_hpi') > 0) then
           call MHEMODL(Lat,mlt,IOr_NeedHPI,IOr_NeedIMFBy, &
                IOr_NeedIMFBz,1,ETheta, EPhi, Potential)
           Potential = Potential * 1000.0
        endif

        if (index(EIE_NameOfEFieldModel,'millstone_imf') > 0) then
           call MHEMODL(lat,mlt,IOr_NeedHPI,IOr_NeedIMFBy, &
                IOr_NeedIMFBz,2,ETheta, EPhi, Potential)
           Potential = Potential * 1000.0
        endif

        if (index(EIE_NameOfEFieldModel,'hmr89') > 0) then
           call HMREPOT(lat,mlt,IOr_NeedIMFBy,IOr_NeedIMFBz, &
                IOr_NeedKp,iChange,ETheta, EPhi, Potential)
           Potential = Potential * 1000.0
           iChange = 0
        endif

        if (index(EIE_NameOfEFieldModel,'izmem') > 0) then
           call IZEPOT(iMonth,lat,mlt, iHemisphere * IOr_NeedIMFBy, &
                IOr_NeedIMFBz,ETheta, EPhi, Potential)
           Potential = Potential * 1000.0
        endif

        PotentialOut(iMLT,iLat) = Potential

     enddo
  enddo

end subroutine IO_GetNonGridBasedPotential


!----------------------------------------------------------------------

subroutine IO_GetAveE(AveEOut, iError)
  
  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs) :: EIE_Value
  real, dimension(IOi_NeednMLTs,IOi_NeednLats)               :: ValueOut
  real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: AveEOut

  iError = 0

  if (UseGridBasedEIE) then

     EIE_Value = EIEr3_HaveAveE
     call IO_GetValue(EIE_Value, ValueOut, 0.1, iError)

     if (iError /= 0) then
        write(*,*) "Error in routine IO_GetAveE:"
        write(*,*) cErrorCodes(iError)
     else
        AveEOut = ValueOut
     endif

  else

     call IO_GetNonGridBasedAveE(ValueOut, iError)

     if (iError == 0) then
        AveEOut = ValueOut
     else
        AveEOut = -1.0e32
     endif

  endif

end subroutine IO_GetAveE

subroutine IO_GetNonGridBasedAveE(AveEOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(IOi_NeednMLTs,IOi_NeednLats) :: AveEOut

  integer :: iMLT, iLat
  real :: Lat, MLT, hpi
  real :: ped, hal, avkev, eflx

  iError = 0
  if (IOr_NeedHPINorm < -1000.0) then
     iError = ecHPINotSet_
     return
  endif

  do iMLT = 1, IOi_NeednMLTs
     do iLat = 1, IOi_NeednLats

        MLT = IOr2_NeedMLTs(iMLT,iLat)
        Lat = IOr2_NeedLats(iMLT,iLat)
        hpi = IOr_NeedHPINorm

        call get_auroral_conductance(Lat, MLT, hpi, &
             ped, hal, avkev, eflx)

        AveEOut(iMLT, iLat) = avkev

     enddo
  enddo

end subroutine IO_GetNonGridBasedAveE

subroutine IO_GetNonGridBasedEFlux(EFluxOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(IOi_NeednMLTs,IOi_NeednLats) :: EFluxOut

  integer :: iMLT, iLat
  real :: Lat, MLT, hpi
  real :: ped, hal, avkev, eflx

  iError = 0

  if (IOr_NeedHPINorm < -1000.0) then
     iError = ecHPINotSet_
     return
  endif

  do iMLT = 1, IOi_NeednMLTs
     do iLat = 1, IOi_NeednLats

        MLT = IOr2_NeedMLTs(iMLT,iLat)
        Lat = IOr2_NeedLats(iMLT,iLat)

        call get_auroral_conductance(Lat, MLT, IOr_NeedHPINorm, &
             ped, hal, avkev, eflx)

        ! Need to convert from erg/cm2/s to W/m2
        EFluxOut(iMLT, iLat) = eflx !* 1.0e-7 * 100.0 * 100.0

     enddo
  enddo

end subroutine IO_GetNonGridBasedEFlux

!----------------------------------------------------------------------

subroutine IO_GetEFlux(EFluxOut, iError)
  
  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs) :: EIE_Value
  real, dimension(IOi_NeednMLTs,IOi_NeednLats)               :: ValueOut
  real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: EFluxOut

  iError = 0

  if (UseGridBasedEIE) then

     EIE_Value = EIEr3_HaveEFlux
     call IO_GetValue(EIE_Value, ValueOut, 1.0e-10, iError)

     if (iError /= 0) then
        write(*,*) "Error in routine IO_GetEFlux:"
        write(*,*) cErrorCodes(iError)
        stop
     else
        EFluxOut = ValueOut
     endif

  else

     call IO_GetNonGridBasedEFlux(ValueOut, iError)

     if (iError == 0) then
        EFluxOut = ValueOut
     else
        EFluxOut = -1.0e32
     endif

  endif

end subroutine IO_GetEFlux

!----------------------------------------------------------------------

subroutine IO_GetValue(ValueIn, ValueOut, Filler, iError)
  
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, intent(in)     :: Filler
  real, dimension(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), intent(in) :: &
       ValueIn
  real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out) :: ValueOut

  integer :: iMLT, iLat, iB, iM, iL
  real    :: dM, dL

  iError = 0

  do iMLT = 1, IOi_NeednMLTs
     do iLat = 1, IOi_NeednLats

        iB = IOi3_InterpolationIndices(iMLT,iLat,1)
        iM = IOi3_InterpolationIndices(iMLT,iLat,2)
        iL = IOi3_InterpolationIndices(iMLT,iLat,3)

        if (iB < 0 .or. iM < 0 .or. iL < 0) then
           ValueOut(iMLT, iLat) =  Filler
        else

           dM = IOr3_InterpolationRatios(iMLT,iLat,1)
           dL = IOr3_InterpolationRatios(iMLT,iLat,2)

           ValueOut(iMLT, iLat) =  &
                (1.0 - dM) * (1.0 - dL) * ValueIn(iM  , IL  , iB) + &
                (1.0 - dM) * (      dL) * ValueIn(iM  , IL+1, iB) + &
                (      dM) * (      dL) * ValueIn(iM+1, IL+1, iB) + &
                (      dM) * (1.0 - dL) * ValueIn(iM+1, IL  , iB)

        endif

     enddo
  enddo

end subroutine IO_GetValue
  
