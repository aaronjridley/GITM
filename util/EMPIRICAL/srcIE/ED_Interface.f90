!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------

subroutine ED_Init(error)

  implicit none

  INCLUDE 'ED_R_elec_ed_lup_subs.inc'

  ! Variables passed:

  integer, intent(out) :: error

  ! External Functions used:

  integer, external :: r_load_edep_lookup

  ! The subroutine:

  call init_all

  error = r_load_edep_lookup ()

end subroutine ED_Init

!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------

subroutine ED_Get_Grid(grid, PressureCoordinates, iError)

  implicit none

  INCLUDE 'ED_R_elec_ed_lup_subs.inc'

  ! Variables passed:
  logical, intent(in) :: PressureCoordinates
  real, dimension(altitudes), intent(out) :: grid(ALTITUDES)
  integer, intent(out) :: iError

  ! External Functions used:
  REAL, external  :: R_EDEP_ALT_VALUE

  ! Local variables:
  integer :: i

  iError = 0

  if (.not.PressureCoordinates) then

     do i=1,altitudes

        ! Get the values in meters (the "2")

        grid(i) = R_EDEP_ALT_VALUE (I,2)

     enddo

  else

     do i=1,altitudes

        ! Get the values in millibars (the "4") and convert to Pascals

        grid(i) = R_EDEP_ALT_VALUE (I,4) * 100.0

     enddo

  endif

end subroutine ED_Get_Grid

!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------

subroutine ED_Get_Grid_Size(npts)

  implicit none

  INCLUDE 'ED_R_elec_ed_lup_subs.inc'

  integer, intent(out) :: npts

  npts = altitudes

end subroutine ED_Get_Grid_Size

!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------

subroutine ED_Get_Number_of_Energies(npts)

  implicit none

  INCLUDE 'ED_R_elec_ed_lup_subs.inc'

  integer, intent(out) :: npts

  npts = spectra_levels

end subroutine ED_Get_Number_of_Energies


!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------

subroutine ED_Get_Energies(energies)

  implicit none

  INCLUDE 'ED_R_elec_ed_lup_subs.inc'

  REAL    EMAT_SAVE(SPECTRA_LEVELS)
  REAL    DEMAT_SAVE(SPECTRA_LEVELS)
  COMMON /R_DE/ EMAT_SAVE,DEMAT_SAVE

  real, dimension(1:spectra_levels), intent(out) :: energies
  real :: ED_Max_Energy, ED_Min_Energy, de
  integer :: i

  ! Now we want to figure out what energies we are going to use:

  ED_Max_Energy = 5.0e5
  ED_Min_Energy = 10.0

  de = (alog10(ED_Max_Energy) - alog10(ED_Min_Energy))/(spectra_levels-1)

  do i=1,spectra_levels
     energies(spectra_levels - (i-1)) = ED_Min_Energy + 10.0**((i-1)*de)
  enddo

end subroutine ED_Get_Energies
