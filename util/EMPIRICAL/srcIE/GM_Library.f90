!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!----------------------------------------------------------------------

subroutine GM_SetnMLTs(nMLTsIn)
  use ModEIE_Interface
  implicit none
  integer, intent(in) :: nMLTsIn
  GMi_NeednMLTs = nMLTsIn
end subroutine GM_SetnMLTs

!----------------------------------------------------------------------

subroutine GM_SetnLats(nLatsIn)
  use ModEIE_Interface
  implicit none
  integer, intent(in) :: nLatsIn
  GMi_NeednLats = nLatsIn
end subroutine GM_SetnLats

!----------------------------------------------------------------------

subroutine GM_SetMLTs(MLTsIn, iError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(GMi_NeednMLTs,GMi_NeednLats), intent(in) :: MLTsIn

  integer :: i,j

  iError = 0
  if (allocated(GMr2_NeedMLTs)) deallocate(GMr2_NeedMLTs)
  allocate(GMr2_NeedMLTs(GMi_NeednMLTs,GMi_NeednLats), stat=iError)
  if (iError /= 0) then
     iError = ecAllocationError_
     return
  else
     do i=1,GMi_NeednMLTs
        do j=1,GMi_NeednLats
           GMr2_NeedMLTs(i,j) = mod((MLTsIn(i,j)+24.0),24.0)
        enddo
     enddo
  endif

end subroutine GM_SetMLTs

!----------------------------------------------------------------------

subroutine GM_SetLats(LatsIn, iError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(GMi_NeednMLTs,GMi_NeednLats), intent(in) :: LatsIn

  iError = 0
  if (allocated(GMr2_NeedLats)) deallocate(GMr2_NeedLats)
  allocate(GMr2_NeedLats(GMi_NeednMLTs,GMi_NeednLats), stat=iError)
  if (iError /= 0) then
     iError = ecAllocationError_
     return
  else
     GMr2_NeedLats(1:GMi_NeednMLTs,1:GMi_NeednLats) = &
          LatsIn(1:GMi_NeednMLTs,1:GMi_NeednLats)
  endif

end subroutine GM_SetLats

!----------------------------------------------------------------------

subroutine GM_SetGrid(MLTsIn, LatsIn, iError)
  
  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(GMi_NeednMLTs,GMi_NeednLats), intent(in) :: MLTsIn, LatsIn
  real, dimension(2) :: GMr1_Location
  real, dimension(5) :: EIEr1_Location

  integer :: i,j

  iError = 0

  call GM_SetMLTs(MLTsIn, iError)
  if (iError /= 0) return

  call GM_SetLats(LatsIn, iError)
  if (iError /= 0) return

  allocate(GMi3_InterpolationIndices(GMi_NeednMLTs,GMi_NeednLats,3), &
       stat=iError)
  if (iError /= 0) then
     iError = ecAllocationError_
     return
  endif

  allocate(GMr3_InterpolationRatios(GMi_NeednMLTs,GMi_NeednLats,2), &
       stat=iError)
  if (iError /= 0) then
     iError = ecAllocationError_
     return
  endif

  do i=1,GMi_NeednLats
     do j=1,GMi_NeednMLTs

        GMr1_Location(1) = GMr2_NeedMLTs(j,i)
        GMr1_Location(2) = GMr2_NeedLats(j,i)

        call EIE_FindPoint(GMr1_Location, EIEr1_Location, iError)

        if (iError == 0) then
           GMi3_InterpolationIndices(j,i,1:3) = EIEr1_Location(1:3)
           GMr3_InterpolationRatios(j,i,1:2)  = EIEr1_Location(4:5)
        else
           GMi3_InterpolationIndices(j,i,1) = -1
        endif

     enddo
  enddo

end subroutine GM_SetGrid

!----------------------------------------------------------------------

subroutine GM_GetPotential(PotentialOut, iError)
  
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs) :: EIE_Value
  real, dimension(GMi_NeednMLTs,GMi_NeednLats)               :: ValueOut
  real, dimension(GMi_NeednMLTs,GMi_NeednLats), intent(out)  :: PotentialOut

  iError = 0

  EIE_Value = EIEr3_HavePotential
  call GM_GetValue(EIE_Value, ValueOut, iError)
  PotentialOut = ValueOut

end subroutine GM_GetPotential

!----------------------------------------------------------------------

subroutine GM_GetAveE(AveEOut, iError)
  
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs) :: EIE_Value
  real, dimension(GMi_NeednMLTs,GMi_NeednLats)               :: ValueOut
  real, dimension(GMi_NeednMLTs,GMi_NeednLats), intent(out)  :: AveEOut

  iError = 0

  EIE_Value = EIEr3_HaveAveE
  call GM_GetValue(EIE_Value, ValueOut, iError)
  AveEOut = ValueOut

end subroutine GM_GetAveE

!----------------------------------------------------------------------

subroutine GM_GetEFlux(EFluxOut, iError)
  
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs) :: EIE_Value
  real, dimension(GMi_NeednMLTs,GMi_NeednLats)               :: ValueOut
  real, dimension(GMi_NeednMLTs,GMi_NeednLats), intent(out)  :: EFluxOut

  iError = 0

  EIE_Value = EIEr3_HaveEFlux
  call GM_GetValue(EIE_Value, ValueOut, iError)
  EFluxOut = ValueOut

end subroutine GM_GetEFlux

!----------------------------------------------------------------------

subroutine GM_GetValue(ValueIn, ValueOut, iError)
  
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), intent(in) :: &
       ValueIn
  real, dimension(GMi_NeednMLTs,GMi_NeednLats), intent(out) :: ValueOut

  integer :: iMLT, iLat, iB, iM, iL
  real    :: dM, dL

  iError = 0

  do iMLT = 1, GMi_NeednMLTs
     do iLat = 1, GMi_NeednLats

        iB = GMi3_InterpolationIndices(iMLT,iLat,1)
        iM = GMi3_InterpolationIndices(iMLT,iLat,2)
        iL = GMi3_InterpolationIndices(iMLT,iLat,3)

        dM = GMr3_InterpolationRatios(iMLT,iLat,1)
        dL = GMr3_InterpolationRatios(iMLT,iLat,2)

        ValueOut(iMLT, iLat) =  &
             (1.0 - dM) * (1.0 - dL) * ValueIn(iM  , IL  , iB) + &
             (      dM) * (1.0 - dL) * ValueIn(iM+1, IL  , iB) + &
             (1.0 - dM) * (      dL) * ValueIn(iM  , IL+1, iB)

     enddo
  enddo

end subroutine GM_GetValue
  
