module ModSpice
  use ModKind, ONLY: Real8_
  ! Empty version of ModSpice_orig.f90 in case there is no SPICE library

  implicit none
  save
  private

  public:: spice_init            ! read in SPICE "kernels", define start time
  public:: spice_rot_matrix      ! return 3x3 rotation matrix
  public:: spice_rot_vel_matrix  ! return 6x6 matrix for position and 
  !                              !      velocity transform
  public:: spice_get_distance    ! return the distance between two bodies

  ! Number of seconds between SWMF and SPICE base times:
  real, parameter,  public::   DtSpiceSwmf = -1104494336.0 

contains
  !============================================================================
  subroutine spice_init(tStart, NameDirIn)

    real(Real8_),     intent(in):: tStart
    character(len=*), intent(in), optional:: NameDirIn
    character(len=100):: NameDir

    character(len=*), parameter:: NameSub = 'spice_init'
    !--------------------------------------------------------------------------

  end subroutine spice_init
  !============================================================================
  subroutine spice_rot_matrix(tSimulation, NameCoord1, NameCoord2, Rot_DD)

    real,             intent(in) :: tSimulation
    character(len=*), intent(in) :: NameCoord1
    character(len=*), intent(in) :: NameCoord2
    real,             intent(out):: Rot_DD(3,3)

    real(Real8_):: tSpice

    character(len=*), parameter:: NameSub = 'spice_rot_matrix'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': SPICE is not swithed on!')

  end subroutine spice_rot_matrix
  !============================================================================
  subroutine spice_rot_vel_matrix(tSimulation, NameCoord1, NameCoord2, Rot_II)

    real,             intent(in) :: tSimulation
    character(len=*), intent(in) :: NameCoord1
    character(len=*), intent(in) :: NameCoord2
    real,             intent(out):: Rot_II(6,6)

    character(len=*), parameter:: NameSub = 'spice_rot_vel_matrix'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': SPICE is not swithed on!')

  end subroutine spice_rot_vel_matrix
  !============================================================================
  subroutine spice_get_distance(tSimulation, NameBody1, NameBody2, Distance)

    real,             intent(in) :: tSimulation
    character(len=*), intent(in) :: NameBody1
    character(len=*), intent(in) :: NameBody2
    real,             intent(out):: Distance

    character(len=*), parameter:: NameSub = 'spice_get_distance'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': SPICE is not swithed on!')

  end subroutine spice_get_distance
end module ModSpice
