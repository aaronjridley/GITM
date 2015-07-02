!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModCms
  use EEE_ModCommonVariables
  implicit none
  save

  private

  public :: set_parameters_cms
  public :: get_cms

  integer :: nLevelCms = 4

  real :: RescaleB = 1.0

  character(len=100):: NameCmsFile

  ! CMS coordinates are respectively: longitude, latitude, radius
  type CmsType
     integer :: nCms_D(3)
     real, allocatable :: CoordCms1_I(:), CoordCms2_I(:), CoordCms3_I(:)
     real, allocatable :: VarCms_VIII(:,:,:,:)
     real :: Min_D(3), Max_D(3)
  end type CmsType

  type(CmsType), allocatable :: CmsData(:)

  logical :: IsFirst = .true.

contains

  !============================================================================

  subroutine set_parameters_cms(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_cms'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#CMS")
       call read_var('NameCmsFile',NameCmsFile)
       call read_var('nLevelCms', nLevelCms)
       call read_var('RescaleB', RescaleB) ! rescale CMS difference field
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_cms

  !============================================================================

  subroutine get_cms(x_D, B_D)

    use ModConst,    ONLY: rSun
    use ModInterpolate, ONLY: trilinear
    use ModCoordTransform, ONLY: xyz_to_sph
    use ModNumConst, ONLY: cDegToRad, cTwoPi, cHalfPi
    use ModPlotFile, ONLY: read_plot_file

    real, intent(in) :: x_D(3)
    real, intent(out) :: B_D(3)

    integer :: iLevel, n_D(3)

    real :: xSph_D(3), xCms_D(3)
    real :: SinTheta, CosTheta, SinPhi, CosPhi
    real :: Bsph_D(3), Bcms_D(3)
    real :: Min_D(3), Max_D(3)

    character(len=100):: NameCmsFileOut
    !--------------------------------------------------------------------------
    B_D = 0.0

    if(IsFirst)then
       allocate(CmsData(nLevelCms))

       ! CMS domain is split in nLevelCms overlapping domains
       do iLevel = 1, nLevelCms
          write(NameCmsFileOut,'(a,i1.1,a)') &
               trim(NameCmsFile)//'_',iLevel-1,'.out'

          call read_plot_file(NameCmsFileOut, nOut_D = CmsData(iLevel)%nCms_D)

          n_D = CmsData(iLevel)%nCms_D

          allocate(CmsData(iLevel)%CoordCms1_I(n_D(1)))
          allocate(CmsData(iLevel)%CoordCms2_I(n_D(2)))
          allocate(CmsData(iLevel)%CoordCms3_I(n_D(3)))
          allocate(CmsData(iLevel)%VarCms_VIII(3,n_D(1),n_D(2),n_D(3)))

          call read_plot_file(NameCmsFileOut, &
               Coord1Out_I = CmsData(iLevel)%CoordCms1_I, &
               Coord2Out_I = CmsData(iLevel)%CoordCms2_I, &
               Coord3Out_I = CmsData(iLevel)%CoordCms3_I, &
               VarOut_VIII = CmsData(iLevel)%VarCms_VIII)

          ! Coordinates are (longitude, latitude, radius)
          ! angular coordinates are in degrees, radius in solar radius
          CmsData(iLevel)%CoordCms1_I = CmsData(iLevel)%CoordCms1_I*cDegToRad
          CmsData(iLevel)%CoordCms2_I = CmsData(iLevel)%CoordCms2_I*cDegToRad
          CmsData(iLevel)%CoordCms3_I = CmsData(iLevel)%CoordCms3_I &
               *rSun*Si2No_V(UnitX_)
          ! Magnetic field is in Gauss
          CmsData(iLevel)%VarCms_VIII = CmsData(iLevel)%VarCms_VIII &
               *1e-4*Si2No_V(UnitB_)

          ! minimum and maximum coordinates for each domain
          CmsData(iLevel)%Min_D = (/ minval(CmsData(iLevel)%CoordCms1_I), &
               minval(CmsData(iLevel)%CoordCms2_I), &
               minval(CmsData(iLevel)%CoordCms3_I) /)
          CmsData(iLevel)%Max_D = (/ maxval(CmsData(iLevel)%CoordCms1_I), &
               maxval(CmsData(iLevel)%CoordCms2_I), &
               maxval(CmsData(iLevel)%CoordCms3_I) /)

       end do

       IsFirst = .false.
    end if

    call xyz_to_sph(x_D, xSph_D)

    xCms_D(1) = xSph_D(3)
    xCms_D(2) = cHalfPi - xSph_D(2)
    xCms_D(3) = xSph_D(1)

    do iLevel = nLevelCms, 1, -1

       n_D = CmsData(iLevel)%nCms_D
       Min_D = CmsData(iLevel)%Min_D
       Max_D = CmsData(iLevel)%Max_D

       ! check if x_D is in CMS domain level. If not, then cycle to next level
       if(xCms_D(3) < Min_D(3) .or. xCms_D(3) > Max_D(3)) CYCLE
       if(xCms_D(2) < Min_D(2) .or. xCms_D(2) > Max_D(2)) CYCLE

       if(xCms_D(1) >= Min_D(1) .and. xCms_D(1) <= Max_D(1))then
       elseif(xCms_D(1)-cTwoPi>=Min_D(1).and.xCms_D(1)-cTwoPi<=Max_D(1))then
          xCms_D(1) = xCms_D(1) - cTwoPi
       elseif(xCms_D(1)+cTwoPi>=Min_D(1).and.xCms_D(1)+cTwoPi<=Max_D(1))then
          xCms_D(1) = xCms_D(1) + cTwoPi
       else
          CYCLE
       end if

       ! Obtain magnetic at x_D via trilinear interpolation of
       ! CMS domin level
       Bcms_D = trilinear(CmsData(iLevel)%VarCms_VIII, 3, &
            1, n_D(1), 1, n_D(2), 1, n_D(3), xCms_D, &
            CmsData(iLevel)%CoordCms1_I, &
            CmsData(iLevel)%CoordCms2_I, &
            CmsData(iLevel)%CoordCms3_I)

       ! CMS magnetic field is (Blongitude, Blatitude, Bradius)
       ! convert to (Bradius, Bcolatitude, Blongitude)
       Bsph_D(1) = Bcms_D(3)
       Bsph_D(2) =-Bcms_D(2)
       Bsph_D(3) = Bcms_D(1)

       SinTheta    = sin(xSph_D(2))
       CosTheta    = cos(xSph_D(2))
       SinPhi      = sin(xSph_D(3))
       CosPhi      = cos(xSph_D(3))

       ! convert to cartesian magnetic field component
       B_D(1) = Bsph_D(1)*SinTheta*CosPhi &
            + Bsph_D(2)*CosTheta*CosPhi - Bsph_D(3)*SinPhi

       B_D(2) =  Bsph_D(1)*SinTheta*SinPhi &
            + Bsph_D(2)*CosTheta*SinPhi + Bsph_D(3)*CosPhi

       B_D(3) = Bsph_D(1)*CosTheta - Bsph_D(2)*SinTheta

       B_D = B_D*No2Si_V(UnitB_)*RescaleB

       EXIT
    end do

  end subroutine get_cms

end module EEE_ModCms
