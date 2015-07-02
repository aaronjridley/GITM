!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! ---------------------------------------------------

subroutine IO_set_IMF_Bz_single(BzIn)

  use ModIndices
  implicit none
  real, intent(in) :: BzIn

  call init_mod_indices

  nIndices_V(imf_bz_) = 1
  Indices_TV(1,imf_bz_) = BzIn

end subroutine IO_set_IMF_Bz_single

! ---------------------------------------------------

subroutine IO_set_IMF_By_single(ByIn)

  use ModIndices
  implicit none
  real, intent(in) :: ByIn

  call init_mod_indices

  nIndices_V(imf_by_) = 1
  Indices_TV(1,imf_by_) = ByIn

end subroutine IO_set_IMF_By_single

! ---------------------------------------------------

subroutine IO_set_IMF_Bx_single(BxIn)

  use ModIndices
  implicit none
  real, intent(in) :: BxIn

  call init_mod_indices

  nIndices_V(imf_bx_) = 1
  Indices_TV(1,imf_bx_) = BxIn

end subroutine IO_set_IMF_Bx_single

! ---------------------------------------------------

subroutine IO_set_SW_v_single(VIn)

  use ModIndices
  implicit none
  real, intent(in) :: VIn

  call init_mod_indices

  nIndices_V(sw_v_) = 1
  Indices_TV(1,sw_v_) = VIn

end subroutine IO_set_SW_v_single

! ---------------------------------------------------

subroutine IO_set_SW_n_single(NIn)

  use ModIndices
  implicit none
  real, intent(in) :: NIn

  call init_mod_indices

  nIndices_V(sw_n_) = 1
  Indices_TV(1,sw_n_) = NIn

end subroutine IO_set_SW_n_single

! ---------------------------------------------------

subroutine IO_set_hpi_single(HpiIn)

  use ModIndices
  implicit none
  real, intent(in) :: HpiIn

  call init_mod_indices

  nIndices_V(Hpi_) = 1
  Indices_TV(1,Hpi_) = HpiIn

end subroutine IO_set_hpi_single

! ---------------------------------------------------

subroutine IO_set_f107_single(F107In)

  use ModIndices
  implicit none
  real, intent(in) :: F107In

  call init_mod_indices

  nIndices_V(f107_) = 1
  Indices_TV(1,f107_) = F107In

end subroutine IO_set_f107_single

! ---------------------------------------------------

subroutine IO_set_f107a_single(F107aIn)

  use ModIndices
  implicit none
  real, intent(in) :: F107aIn

  call init_mod_indices

  nIndices_V(f107a_) = 1
  Indices_TV(1,f107a_) = F107aIn

end subroutine IO_set_f107a_single

! ---------------------------------------------------

subroutine IO_set_kp_single(kpIn)

  use ModIndices
  implicit none
  real, intent(in) :: kpIn

  call init_mod_indices

  nIndices_V(kp_) = 1
  Indices_TV(1,kp_) = kpIn

end subroutine IO_set_kp_single


! ---------------------------------------------------

subroutine IO_set_ap_single(apIn)

  use ModIndices
  implicit none
  real, intent(in) :: apIn

  call init_mod_indices

  nIndices_V(ap_) = 1
  Indices_TV(1,ap_) = apIn

end subroutine IO_set_ap_single


