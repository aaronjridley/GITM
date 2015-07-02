!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

Module ModIndicesInterfaces

  !\
  ! IMF Interfaces
  !/

  interface get_IMF_Bx
     subroutine get_IMF_Bx_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_IMF_Bx_wtime
     subroutine get_IMF_Bx_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_IMF_Bx_wotime
  end interface

  interface get_IMF_By
     subroutine get_IMF_By_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_IMF_By_wtime
     subroutine get_IMF_By_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_IMF_By_wotime
  end interface

  interface get_IMF_Bz
     subroutine get_IMF_Bz_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_IMF_Bz_wtime
     subroutine get_IMF_Bz_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_IMF_Bz_wotime
  end interface

  interface get_IMF_B
     subroutine get_IMF_B_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_IMF_B_wtime
     subroutine get_IMF_B_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_IMF_B_wotime
  end interface

  !\
  ! Solar Wind Interfaces
  !/

  interface get_SW_Vx
     subroutine get_SW_Vx_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_Vx_wtime
     subroutine get_SW_Vx_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_Vx_wotime
  end interface

  interface get_SW_Vy
     subroutine get_SW_Vy_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_Vy_wtime
     subroutine get_SW_Vy_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_Vy_wotime
  end interface

  interface get_SW_Vz
     subroutine get_SW_Vz_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_Vz_wtime
     subroutine get_SW_Vz_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_Vz_wotime
  end interface

  interface get_SW_V
     subroutine get_SW_V_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_V_wtime
     subroutine get_SW_V_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_V_wotime
  end interface

  !\
  ! Solar Wind N
  !/

  interface get_SW_N
     subroutine get_SW_N_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_N_wtime
     subroutine get_SW_N_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_SW_N_wotime
  end interface

  !\
  ! F10.7
  !/

  interface get_f107
     subroutine get_f107_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_f107_wtime
     subroutine get_f107_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_f107_wotime
  end interface

  !\
  ! F10.7a
  !/

  interface get_f107a
     subroutine get_f107a_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_f107a_wtime
     subroutine get_f107a_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_f107a_wotime
  end interface

  !\
  ! Hemispheric Power Index
  !/

  interface get_hpi
     subroutine get_hpi_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_hpi_wtime
     subroutine get_hpi_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_hpi_wotime
  end interface

  interface get_hpi_calc
     subroutine get_hpi_calc_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_hpi_calc_wtime
     subroutine get_hpi_calc_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_hpi_calc_wotime
  end interface

  interface get_hpi_norm
     subroutine get_hpi_norm_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_hpi_norm_wtime
     subroutine get_hpi_norm_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_hpi_norm_wotime
  end interface

  !\
  ! kp
  !/

  interface get_kp
     subroutine get_kp_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_kp_wtime
     subroutine get_kp_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_kp_wotime
  end interface

  !\
  ! ap
  !/

  interface get_ap
     subroutine get_ap_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_ap_wtime
     subroutine get_ap_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_ap_wotime
  end interface

  !\
  ! Auroral Indices
  !/

  interface get_au
     subroutine get_au_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_au_wtime
     subroutine get_au_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_au_wotime
  end interface

  interface get_al
     subroutine get_al_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_al_wtime
     subroutine get_al_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_al_wotime
  end interface

  interface get_ae
     subroutine get_ae_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_ae_wtime
     subroutine get_ae_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_ae_wotime
  end interface

  interface get_nAE
     subroutine get_nAE_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_nAE_wtime
     subroutine get_nAE_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_nAE_wotime
  end interface

  interface get_onsetut
     subroutine get_onsetut_wtime(TimeIn, iNext, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       integer, intent(in)  :: iNext
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_onsetut_wtime
     subroutine get_onsetut_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_onsetut_wotime
  end interface

  !\
  ! Dst Interfaces
  !/

  interface get_Dst
     subroutine get_Dst_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_Dst_wtime
     subroutine get_Dst_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_Dst_wotime
  end interface

  interface get_nDst
     subroutine get_nDst_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_nDst_wtime
     subroutine get_nDst_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_nDst_wotime
  end interface

  !\
  ! Joule Heating Interface
  !/

  interface get_jh_calc
     subroutine get_jh_calc_wtime(TimeIn, value, iOutputError)
       use ModKind
       real (Real8_)  :: TimeIn
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_jh_calc_wtime
     subroutine get_jh_calc_wotime(value, iOutputError)
       real, intent(out)    :: value
       integer, intent(out) :: iOutputError
     end subroutine get_jh_calc_wotime
  end interface

end Module ModIndicesInterfaces
