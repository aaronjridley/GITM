!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! ---------------------------------------------------

subroutine indices_set_IMF_Bz_single(BzIn)

  use ModIndices
  implicit none
  real, intent(in) :: BzIn

  nIndices_V(imf_bz_) = 1
  Indices_TV(1,imf_bz_) = BzIn

end subroutine indices_set_IMF_Bz_single

! ---------------------------------------------------

subroutine indices_set_IMF_By_single(ByIn)

  use ModIndices
  implicit none
  real, intent(in) :: ByIn

  nIndices_V(imf_by_) = 1
  Indices_TV(1,imf_by_) = ByIn

end subroutine indices_set_IMF_By_single

! ---------------------------------------------------

subroutine indices_set_IMF_Bx_single(BxIn)

  use ModIndices
  implicit none
  real, intent(in) :: BxIn

  nIndices_V(imf_bx_) = 1
  Indices_TV(1,imf_bx_) = BxIn

end subroutine indices_set_IMF_Bx_single

! ---------------------------------------------------

subroutine indices_set_SW_v_single(VIn)

  use ModIndices
  implicit none
  real, intent(in) :: VIn

  nIndices_V(sw_v_) = 1
  Indices_TV(1,sw_v_) = VIn

end subroutine indices_set_SW_v_single

! ---------------------------------------------------

subroutine indices_set_hpi_single(HpiIn)

  use ModIndices
  implicit none
  real, intent(in) :: HpiIn

  nIndices_V(Hpi_) = 1
  Indices_TV(1,Hpi_) = HpiIn

end subroutine indices_set_hpi_single

! ---------------------------------------------------

subroutine indices_set_f107_single(F107In)

  use ModIndices
  implicit none
  real, intent(in) :: F107In

  nIndices_V(f107_) = 1
  Indices_TV(1,f107_) = F107In

end subroutine indices_set_f107_single

! ---------------------------------------------------

subroutine indices_set_f107a_single(F107aIn)

  use ModIndices
  implicit none
  real, intent(in) :: F107aIn

  nIndices_V(f107a_) = 1
  Indices_TV(1,f107a_) = F107aIn

end subroutine indices_set_f107a_single

! ---------------------------------------------------

subroutine indices_set_kp_single(kpIn)

  use ModIndices
  implicit none
  real, intent(in) :: kpIn

  nIndices_V(kp_) = 1
  Indices_TV(1,kp_) = kpIn

end subroutine indices_set_kp_single

! ---------------------------------------------------

subroutine set_time(TimeIn, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  integer, intent(out) :: iOutputError
  integer :: iIndex, iError
  real    :: IndexValue

  iOutputError = 0

  !\
  ! If values are already set, then do nothing
  !/

  if (SavedTime == TimeIn) return

  SavedTime = TimeIn

  SavedIndices_V = -1.0e32
  SavedErrors_V  = 0

  do iIndex=1,nIndices
     call get_index(iIndex, SavedTime, 0, IndexValue, iError)
     SavedIndices_V(iIndex) = IndexValue
     SavedErrors_V(iIndex)  = iError
  enddo

end subroutine set_time

! If Interpolate = 0, it will interpolate
! If Interpolare = 1, it will take the past value
! If Interpolate = 2, it will take the next value

subroutine get_index(label_, TimeIn, Interpolate, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  integer, intent(in)  :: label_
  real (Real8_)  :: TimeIn
  integer, intent(in)  :: Interpolate
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  integer :: iMin, iMax, iCenter
  real    :: DtNorm
  logical :: IsFound

  iOutputError = 0
  value = -1.0e32

  if (label_ < 1 .or. label_ > nIndices) then
     iOutputError = 1
     return
  endif

  if (nIndices_V(label_) == 0) then
     iOutputError = 2
     return
  endif
     
  if (nIndices_V(label_) == 1) then
     value = Indices_TV(1,label_)
     return
  endif

  if (TimeIn <= IndexTimes_TV(1,label_)) then
     value = Indices_TV(1,label_)
     return
  endif

  if (TimeIn >= IndexTimes_TV(nIndices_V(label_),label_)) then
     value = Indices_TV(nIndices_V(label_),label_)
     return
  endif

  iMin = 1
  iMax = nIndices_V(label_)
  IsFound = .false.

  do while (.not.IsFound)

     iCenter = (iMin + iMax)/2

     if (iCenter >= iMax .or. iCenter <= iMin) then
        IsFound = .true.
     else

        if (TimeIn == IndexTimes_TV(iCenter, label_)) then
           iMin = iCenter
           iMax = iCenter
        else

           if (TimeIn < IndexTimes_TV(iCenter, label_)) then
              iMax = iCenter
           else
              iMin = iCenter
           endif

        endif

     endif

  enddo

  if (Interpolate == 1) then 
     value = Indices_TV(iMin, label_)
  elseif (Interpolate == 2) then 
     value = Indices_TV(iMin+1, label_)
  else

     if (iMin == iMax) then
        value = Indices_TV(iCenter, label_)
     else
        DtNorm = 1.0 - (IndexTimes_TV(iMax, label_) - TimeIn) / &
             (IndexTimes_TV(iMax, label_) - IndexTimes_TV(iMin, label_)+1.0e-6)
        value  =  DtNorm * Indices_TV(iMax, label_) + &
             (1.0 - DtNorm) * Indices_TV(iMin, label_)
     endif

  endif

end subroutine get_index

!------------------------------------------------------------------------------
! IMF
!------------------------------------------------------------------------------

subroutine get_IMF_Bx_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  call get_index(imf_bx_, TimeIn, 0, value, iOutputError)

end subroutine get_IMF_Bx_wtime

subroutine get_IMF_Bx_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(imf_bx_)
     iOutputError = SavedErrors_V(imf_bx_)
  endif
  
end subroutine get_IMF_Bx_wotime

subroutine get_IMF_By_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  call get_index(imf_by_, TimeIn, 0, value, iOutputError)

end subroutine get_IMF_By_wtime

subroutine get_IMF_By_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(imf_by_)
     iOutputError = SavedErrors_V(imf_by_)
  endif
  
end subroutine get_IMF_By_wotime

subroutine get_IMF_Bz_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  call get_index(imf_bz_, TimeIn, 0, value, iOutputError)

end subroutine get_IMF_Bz_wtime

subroutine get_IMF_Bz_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(imf_bz_)
     iOutputError = SavedErrors_V(imf_bz_)
  endif
  
end subroutine get_IMF_Bz_wotime

subroutine get_IMF_B_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  real :: value_x, value_y, value_z

  iOutputError = 0

  call get_index(imf_bx_, TimeIn, 0, value_x, iOutputError)
  call get_index(imf_by_, TimeIn, 0, value_y, iOutputError)
  call get_index(imf_bz_, TimeIn, 0, value_z, iOutputError)

  if (iOutputError == 0) then
     value = sqrt(value_x**2 + value_y**2 + value_z**2)
  else
     value = -1.0e32
  endif

end subroutine get_IMF_B_wtime

subroutine get_IMF_B_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  real :: value_x, value_y, value_z

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value_x = SavedIndices_V(imf_bx_)
     value_y = SavedIndices_V(imf_by_)
     value_z = SavedIndices_V(imf_bz_)
     iOutputError = SavedErrors_V(imf_bx_)
  endif

  if (iOutputError == 0) then
     value = sqrt(value_x**2 + value_y**2 + value_z**2)
  else
     value = -1.0e32
  endif

end subroutine get_IMF_B_wotime

!------------------------------------------------------------------------------
! Solar Wind Velocity
!------------------------------------------------------------------------------

!\
! X Component
!/

subroutine get_SW_Vx_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  call get_index(sw_vx_, TimeIn, 0, value, iOutputError)

end subroutine get_SW_Vx_wtime

subroutine get_SW_Vx_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(sw_vx_)
     iOutputError = SavedErrors_V(sw_vx_)
  endif
  
end subroutine get_SW_Vx_wotime

!\
! Y Component
!/

subroutine get_SW_Vy_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  call get_index(sw_vy_, TimeIn, 0, value, iOutputError)

end subroutine get_SW_Vy_wtime

subroutine get_SW_Vy_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(sw_vy_)
     iOutputError = SavedErrors_V(sw_vy_)
  endif
  
end subroutine get_SW_Vy_wotime

!\
! Z Component
!/

subroutine get_SW_Vz_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  call get_index(sw_vz_, TimeIn, 0, value, iOutputError)

end subroutine get_SW_Vz_wtime

subroutine get_SW_Vz_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(sw_vz_)
     iOutputError = SavedErrors_V(sw_vz_)
  endif
  
end subroutine get_SW_Vz_wotime

!\
! Full V
!/

subroutine get_SW_V_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_)  :: TimeIn
  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  real :: value_x, value_y, value_z

  iOutputError = 0

  !\
  ! this is a little odd, since you can have V as an input, but it
  ! is probably better to see if it is a combo of Vx, Vy, and Vz first
  !/

  call get_index(sw_vx_, TimeIn, 0, value_x, iOutputError)

  if (iOutputError > 0) then 

     !\
     ! This means that we don't have the 3 components, but we may
     ! have the full velocity anyways.
     !/

     call get_index(sw_v_, TimeIn, 0, value_x, iOutputError)
     value = value_x

  else
     call get_index(sw_vy_, TimeIn, 0, value_y, iOutputError)
     call get_index(sw_vz_, TimeIn, 0, value_z, iOutputError)
     value = sqrt(value_x**2 + value_y**2 + value_z**2)
  endif

end subroutine get_SW_V_wtime

subroutine get_SW_V_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  real :: value_x, value_y, value_z

  iOutputError = 0

  !\
  ! this is a little odd, since you can have V as an input, but it
  ! is probably better to see if it is a combo of Vx, Vy, and Vz first
  !/

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else

     value_x = SavedIndices_V(sw_vx_)
     iOutputError = SavedErrors_V(sw_vx_)

     if (iOutputError > 0) then 

        !\
        ! This means that we don't have the 3 components, but we may
        ! have the full velocity anyways.
        !/

        value = SavedIndices_V(sw_v_)
        iOutputError = SavedErrors_V(sw_v_)

     else
        value_y = SavedIndices_V(sw_vy_)
        value_z = SavedIndices_V(sw_vz_)
        value = sqrt(value_x**2 + value_y**2 + value_z**2)
     endif

  endif

end subroutine get_SW_V_wotime

!------------------------------------------------------------------------------
! Solar Wind Density
!------------------------------------------------------------------------------

subroutine get_SW_N_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in)  :: TimeIn
  real, intent(out)                :: value
  integer, intent(out)             :: iOutputError

  call get_index(sw_n_, TimeIn, 0, value, iOutputError)

end subroutine get_SW_N_wtime

subroutine get_SW_N_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(sw_n_)
     iOutputError = SavedErrors_V(sw_n_)
  endif
  
end subroutine get_SW_N_wotime

!------------------------------------------------------------------------------
! F10.7
!------------------------------------------------------------------------------

subroutine get_f107_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  iOutputError = 0

  call get_index(f107_, TimeIn, 0, value, iOutputError)

end subroutine get_f107_wtime

subroutine get_f107_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(f107_)
     iOutputError = SavedErrors_V(f107_)
  endif
  
end subroutine get_f107_wotime

!------------------------------------------------------------------------------
! F10.7a
!------------------------------------------------------------------------------

subroutine get_f107a_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(f107a_, TimeIn, 0, value, iOutputError)

end subroutine get_f107a_wtime

subroutine get_f107a_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(f107a_)
     iOutputError = SavedErrors_V(f107a_)
  endif
  
end subroutine get_f107a_wotime

!------------------------------------------------------------------------------
! Hemispheric Power Index
!------------------------------------------------------------------------------

subroutine get_hpi_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(hpi_, TimeIn, 0, value, iOutputError)

end subroutine get_hpi_wtime

subroutine get_hpi_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(hpi_)
     iOutputError = SavedErrors_V(hpi_)
  endif
  
end subroutine get_hpi_wotime

subroutine get_hpi_calc_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(hpi_calc_, TimeIn, 0, value, iOutputError)

end subroutine get_hpi_calc_wtime

subroutine get_hpi_calc_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(hpi_calc_)
     iOutputError = SavedErrors_V(hpi_calc_)
  endif
  
end subroutine get_hpi_calc_wotime

subroutine get_hpi_norm_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(hpi_norm_, TimeIn, 0, value, iOutputError)

  if (iOutputError > 0) then
     call get_index(hpi_, TimeIn, 0, value, iOutputError)
     if (iOutputError == 0) then
        value = 2.09 * ALOG(value) * 1.0475
     endif
  endif

end subroutine get_hpi_norm_wtime

subroutine get_hpi_norm_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(hpi_norm_)
     iOutputError = SavedErrors_V(hpi_norm_)
  endif
  
end subroutine get_hpi_norm_wotime

!------------------------------------------------------------------------------
! kp
!------------------------------------------------------------------------------

subroutine get_kp_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(kp_, TimeIn, 0, value, iOutputError)

end subroutine get_kp_wtime

subroutine get_kp_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(kp_)
     iOutputError = SavedErrors_V(kp_)
  endif
  
end subroutine get_kp_wotime

!------------------------------------------------------------------------------
! ap
!------------------------------------------------------------------------------

subroutine get_ap_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(ap_, TimeIn, 0, value, iOutputError)

end subroutine get_ap_wtime

subroutine get_ap_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(ap_)
     iOutputError = SavedErrors_V(ap_)
  endif
  
end subroutine get_ap_wotime

!------------------------------------------------------------------------------
! Auroral Indices
!------------------------------------------------------------------------------

subroutine get_au_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(au_, TimeIn, 0, value, iOutputError)

end subroutine get_au_wtime

subroutine get_au_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(au_)
     iOutputError = SavedErrors_V(au_)
  endif
  
end subroutine get_au_wotime

subroutine get_al_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(al_, TimeIn, 0, value, iOutputError)

end subroutine get_al_wtime

subroutine get_al_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(al_)
     iOutputError = SavedErrors_V(al_)
  endif
  
end subroutine get_al_wotime

subroutine get_ae_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(ae_, TimeIn, 0, value, iOutputError)

end subroutine get_ae_wtime

subroutine get_ae_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(ae_)
     iOutputError = SavedErrors_V(ae_)
  endif
  
end subroutine get_ae_wotime

subroutine get_onsetut_wtime(TimeIn, iNext, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in)       :: TimeIn
  integer, intent(in)             :: iNext
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  iOutputError = 0

  call get_index(onsetut_, TimeIn, iNext, value, iOutputError)

end subroutine get_onsetut_wtime

subroutine get_onsetut_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(onsetut_)
     iOutputError = SavedErrors_V(onsetut_)
  endif
  
end subroutine get_onsetut_wotime

subroutine get_nAE_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(nAE_, TimeIn, 0, value, iOutputError)

end subroutine get_nAE_wtime

subroutine get_nAE_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(nAE_)
     iOutputError = SavedErrors_V(nAE_)
  endif
  
end subroutine get_nAE_wotime

!------------------------------------------------------------------------------
! Dst
!------------------------------------------------------------------------------

subroutine get_dst_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(dst_, TimeIn, 0, value, iOutputError)

end subroutine get_dst_wtime

subroutine get_dst_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(dst_)
     iOutputError = SavedErrors_V(dst_)
  endif
  
end subroutine get_dst_wotime

subroutine get_nDst_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(nDst_, TimeIn, 0, value, iOutputError)

end subroutine get_nDst_wtime

subroutine get_nDst_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(nDst_)
     iOutputError = SavedErrors_V(nDst_)
  endif
  
end subroutine get_nDst_wotime

!------------------------------------------------------------------------------
! Joule Heating
!------------------------------------------------------------------------------

subroutine get_jh_calc_wtime(TimeIn, value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real (Real8_), intent(in) :: TimeIn
  real, intent(out)               :: value
  integer, intent(out)            :: iOutputError

  call get_index(jh_calc_, TimeIn, 0, value, iOutputError)

end subroutine get_jh_calc_wtime

subroutine get_jh_calc_wotime(value, iOutputError)

  use ModKind
  use ModIndices

  implicit none

  real, intent(out)    :: value
  integer, intent(out) :: iOutputError

  iOutputError = 0

  if (SavedTime < 0.0) then
     value = -1.0e32
     iOutputError = 3
     return
  else
     value = SavedIndices_V(jh_calc_)
     iOutputError = SavedErrors_V(jh_calc_)
  endif
  
end subroutine get_jh_calc_wotime

