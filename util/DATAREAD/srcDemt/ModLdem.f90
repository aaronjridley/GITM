!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module ModLdem

  ! Reads the 3D LDEM (local differential emission measure) moments.
  ! Currently, only the lowest measured layer is used (usely at 1.035 Rsun)
  ! The rejected data points (negative data) are filled in by interpolation
  ! using the Delaunay triangulation. The resulting boundary electron density
  ! and mean electron temperature can be used as boundary condition for the
  ! solar wind.

  implicit none
  save

  private !except

  public :: read_ldem
  public :: get_ldem_moments

  logical, public :: UseLdem = .false.
  character(len=100), public :: NameLdemFile
  integer, public :: iRadiusLdem = 1

  real, allocatable :: ThetaLdem_I(:), PhiLdem_I(:), VarLdem_VII(:,:,:)
  integer :: nThetaLdem, nPhiLdem

contains

  subroutine read_ldem(NamePlotDir)

    use ModNumConst,    ONLY: cDegToRad, cRadToDeg, cHalfPi, cTwoPi
    !use ModPhysics,     ONLY: AverageIonCharge, Si2No_V, No2Si_V, UnitN_, &
    !     UnitTemperature_
    use ModPlotFile,    ONLY: read_plot_file, save_plot_file
    use ModMpi     
    use ModTriangulate, ONLY: calc_triangulation, find_triangle

    character(len=*),intent(in) :: NamePlotDir
      
    ! Original 3D EUVI LDEM data
    integer :: nCell_D(3), nVar
    real, allocatable :: Coord_DIII(:,:,:,:), Var_VIII(:,:,:,:)

    ! 2D boundary map
    integer :: j, k, iVar
    real, allocatable :: Coord_DII(:,:,:), Var_VII(:,:,:)

    ! Triangulation of boundary map to patch the LDEM map for used radius
    integer :: iCell, nCell, iNode1, iNode2, iNode3
    integer, save :: nTriangle
    real :: Weight1, Weight2, Weight3
    integer, allocatable :: iNodeTriangle_II(:,:)
    real, allocatable :: Coord_DI(:,:), Var_VI(:,:)

    integer, parameter :: Theta_ = 2, Phi_ = 3
  
    integer  :: iProc, iError
    character (len=*), parameter :: NameSub = 'read_ldem'
    !--------------------------------------------------------------------------

    nCell_D = 1
    call read_plot_file(NameLdemFile, nVarOut = nVar, nOut_D = nCell_D)

    if(iRadiusLdem < 1 .or. iRadiusLdem > nCell_D(1))then
       write(*,*) NameSub//' ERROR: iRadiusLdem is out of range [1,', &
            nCell_D(1),']'
       call CON_stop(' Correct PARAM.in file')
    end if
       
    allocate(Coord_DIII(3,nCell_D(1),nCell_D(2),nCell_D(3)), &
         Var_VIII(nVar,nCell_D(1),nCell_D(2),nCell_D(3)))

    call read_plot_file(NameLdemFile, &
         CoordOut_DIII = Coord_DIII, VarOut_VIII = Var_VIII)

    allocate(Coord_DII(Theta_:Phi_,nCell_D(2),nCell_D(3)), &
         Var_VII(2,nCell_D(2),nCell_D(3)))

    ! For now: the smallest radius is assumed to be the solar boundary
    Coord_DII(:,:,:) = Coord_DIII(Theta_:Phi_,iRadiusLdem,:,:)

    ! Keep only the two lowest moments
    Var_VII(:,:,:) = Var_VIII(1:2,iRadiusLdem,:,:)

    ! 3D data is no longer needed
    deallocate(Coord_DIII, Var_VIII)

    ! How many elements have only positive data ?
    nCell = count(Var_VII(1,:,:)>0.0 .and. Var_VII(2,:,:)>0.0)

    ! Create arrays for Delaunay triangulation
    allocate(Coord_DI(Theta_:Phi_,nCell), Var_VI(2,nCell), &
         iNodeTriangle_II(3,2*nCell))

    iCell = 0
    do k = 1, nCell_D(3); do j = 1, nCell_D(2)
       if(any(Var_VII(1:2,j,k) <= 0.0)) CYCLE

       iCell = iCell + 1
       Coord_DI(:,iCell) = Coord_DII(:,j,k)
       Var_VI(:,iCell) = Var_VII(:,j,k)
    end do; end do

    ! Delaunay triangulation
    call calc_triangulation(nCell, Coord_DI, iNodeTriangle_II, nTriangle)

    ! Overwrite the rejected (negative) data points
    do k = 1, nCell_D(3); do j = 1, nCell_D(2)
       if(all(Var_VII(1:2,j,k) > 0.0)) CYCLE

       call find_triangle(nCell, nTriangle, Coord_DII(:,j,k), &
            Coord_DI(:,:), iNodeTriangle_II(:,1:nTriangle), &
            iNode1, iNode2, iNode3, Weight1, Weight2, Weight3)

       Var_VII(:,j,k) = &
            Weight1*Var_VI(:,iNode1) + &
            Weight2*Var_VI(:,iNode2) + &
            Weight3*Var_VI(:,iNode3)
    end do; end do

    ! We are done with the Delaunay triangulation
    deallocate(Coord_DI, Var_VI, iNodeTriangle_II)

    ! Store coordinate arrays for bilinear interpolation
    ! extend the domain size and make dimensionless
    nThetaLdem = nCell_D(2)
    nPhiLdem = nCell_D(3)
    allocate(ThetaLdem_I(0:nThetaLdem+1), PhiLdem_I(0:nPhiLdem+1))
    ThetaLdem_I(1:nThetaLdem) = Coord_DII(Theta_,:,1)*cDegToRad
    ThetaLdem_I(0) = -cHalfPi
    ThetaLdem_I(nThetaLdem+1) = cHalfPi
    PhiLdem_I(1:nPhiLdem) = Coord_DII(Phi_,1,:)*cDegToRad
    PhiLdem_I(0) = PhiLdem_I(nPhiLdem) - cTwoPi
    PhiLdem_I(nPhiLdem+1) = PhiLdem_I(1) + cTwoPi

    allocate(VarLdem_VII(2,0:nThetaLdem+1,0:nPhiLdem+1))

    VarLdem_VII(:,1:nThetaLdem,1:nPhiLdem) = Var_VII(:,:,:)
    ! periodicity
    VarLdem_VII(:,1:nThetaLdem,0) = Var_VII(:,:,nCell_D(3))
    VarLdem_VII(:,1:nThetaLdem,nPhiLdem+1) = Var_VII(:,:,1)
    ! average to poles
    do iVar = 1, 2
       VarLdem_VII(iVar,0,:) = sum(Var_VII(iVar,1,:))/nCell_D(3)
       VarLdem_VII(iVar,nThetaLdem+1,:) = &
            sum(Var_VII(iVar,nThetaLdem,:))/nCell_D(3)
    end do

    ! Convert moments (electron density and temperature) to SI
    ! density:     1/cm^3 -> 1/m^3
    ! temperature: MK     -> K
    VarLdem_VII(:,:,:) = VarLdem_VII(:,:,:)*1e6
  
    deallocate(Coord_DII, Var_VII)

    call MPI_COMM_RANK(MPI_COMM_WORLD,iProc,iError)

    if(iProc == 0)then
       ! Flip coordinates for plotting
       ! (The original file has a left handed coordinate system)
       allocate(Coord_DII(2,0:nPhiLdem+1,0:nThetaLdem+1), &
            Var_VII(2,0:nPhiLdem+1,0:nThetaLdem+1))
       do k = 0, nThetaLdem+1; do j = 0, nPhiLdem+1
          Coord_DII(1,j,k) = PhiLdem_I(nPhiLdem+1-j)*cRadToDeg
          Coord_DII(2,j,k) = ThetaLdem_I(nThetaLdem+1-k)*cRadToDeg
          ! convert to output file units
          Var_VII(1,j,k) = VarLdem_VII(1,nThetaLdem+1-k,nPhiLdem+1-j)*1e-6
          Var_VII(2,j,k) = VarLdem_VII(2,nThetaLdem+1-k,nPhiLdem+1-j)*1e-6
       end do; end do

       ! Show the patched LDEM moments
       call save_plot_file(trim(NamePlotDir)//'LDEM_patched.outs',&
            nDimIn = 2, CoordIn_DII = Coord_DII, VarIn_VII = Var_VII, &
            NameVarIn = 'Longitude Latitude Ne Tm', StringHeaderIn = &
            'Longitude [Deg], Latitude [Deg], Ne [1/cm^3], Te [MK]')

       deallocate(Coord_DII, Var_VII)
    end if

  end subroutine read_ldem

  !============================================================================

  subroutine get_ldem_moments(x_D, Ne, Te)

    Use ModInterpolate, ONLY: bilinear
    use ModNumConst,    ONLY: cHalfPi, cTwoPi

    real, intent(in)  :: x_D(3)
    real, intent(out) :: Ne, Te

    real :: r, Theta, Phi, VarLdem_V(2)
    !--------------------------------------------------------------------------
    r = sqrt(sum(x_D**2))
    Theta = cHalfPi - acos(x_D(3)/r)
    Phi = modulo(atan2(x_D(2),x_D(1)), cTwoPi)

    VarLdem_V = bilinear(VarLdem_VII, 2, 0, nThetaLdem+1, 0, nPhiLdem+1, &
         (/Theta,Phi/), x_I=ThetaLdem_I, y_I=PhiLdem_I, DoExtrapolate=.false.)

    Ne = VarLdem_V(1)
    Te = VarLdem_V(2)

  end subroutine get_ldem_moments

end module ModLdem
