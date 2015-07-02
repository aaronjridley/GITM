!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPotentialField
  
  use ModMpi

  implicit none

  ! input parameter
  logical:: DoReadMagnetogram = .true.

  ! grid and domain parameters
  integer:: nR = 150, nThetaAll = 180, nPhiAll = 360
  real   :: rMin = 1.0, rMax = 2.5
  logical:: UseLogRadius = .false.  ! logarithmic or linear in radius

  ! domain decomposition
  integer:: nProcTheta = 2, nProcPhi = 2

  ! solver parameters
  character(len=20):: NameSolver         = 'BICGSTAB' ! or 'GMRES' or 'AMG'
  character(len=20):: NamePreconditioner = 'ILU' ! or 'AMG' or 'MG'
  logical          :: UseHypre           = .false.
  logical          :: UsePreconditioner  = .true.
  real             :: Tolerance          = 1e-10

  ! magnetogram parameters
  logical           :: IsNewMagnetogramStyle = .false.
  character(len=100):: NameFileIn = 'fitsfile.dat'  ! filename
  logical           :: UseCosTheta  = .true. 
  real              :: BrMax = 3500.0               ! Saturation level of MDI

  ! output paramters
  logical           :: DoSaveField   = .true.
  character(len=100):: NameFileField = 'potentialfield'
  character(len=5)  :: TypeFileField = 'real8'

  logical           :: DoSavePotential   = .true.
  character(len=100):: NameFilePotential = 'potentialtest'
  character(len=5)  :: TypeFilePotential = 'real8'
  
  logical           :: DoSaveTecplot   = .false.
  character(len=100):: NameFileTecplot = 'potentialfield.dat'

  ! testing parameters
  logical :: UseTiming = .true.
  real    :: TimeStart, TimeEnd

  logical :: DoTestMe  = .false.
  integer :: iRTest = 1, iPhiTest = 1, iThetaTest = 2

  ! local variables --------------------
  character(len=100):: NameFile

  logical :: UseBr = .true.

  real, dimension(:), allocatable :: &
       Radius_I, Theta_I, Phi_I, SinTheta_I, &
       dRadius_I, dPhi_I, dCosTheta_I, &
       RadiusNode_I, ThetaNode_I, PhiNode_I, SinThetaNode_I, &
       dRadiusNode_I, dTheta_I, dThetaNode_I, dPhiNode_I, dCosThetaNode_I

  real, allocatable:: Br_II(:,:), Potential_C(:,:,:), Rhs_C(:,:,:), &
       B0_DF(:,:,:,:), DivB_C(:,:,:), PlotVar_VC(:,:,:,:), BrLocal_II(:,:)

  ! Variables for hepta preconditioner
  real, parameter:: PrecondParam = 1.0 ! see ModLinearSolver

  ! Seven diagonals for the preconditioner
  real, dimension(:), allocatable :: &
       d_I, e_I, e1_I, e2_I, f_I, f1_I, f2_I

  integer, parameter :: iComm = MPI_COMM_WORLD
  integer :: iProc, nProc, iProcTheta, iProcPhi
  integer :: iTheta0, iPhi0
  integer :: nTheta, nPhi
  real,  allocatable :: tmpXPhi0_II(:,:),tmpXPhipi_II(:,:)
  integer :: nThetaLgr,nThetaSml,nPhiLgr,nPhiSml
  integer :: nProcThetaLgr,nProcThetaSml,nProcPhiLgr,nProcPhiSml

  real, allocatable :: &
       sendBC010_II(:,:), sendBC180_II(:,:), sendBC12_II(:,:), &
       sendBC21_II(:,:),  sendBC34_II(:,:), sendBC43_II(:,:), &
       sendBC020_II(:,:), sendBC360_II(:,:), &
       recvBCLgr010_II(:,:), recvBCSml010_II(:,:), &
       recvBCLgr180_II(:,:), recvBCSml180_II(:,:), &
       recvBC12_II(:,:), recvBC21_II(:,:), &
       recvBC34_II(:,:), recvBC43_II(:,:), &
       recvBC020_II(:,:), recvBC360_II(:,:)

contains

  !===========================================================================
  subroutine read_fdips_param

    use ModReadParam

    character(len=lStringLine) :: NameCommand
    character(len=10):: TypeOutput
    integer:: i, iProcTest

    character(len=*), parameter:: NameSub = 'read_fdips_param'
    !-----------------------------------------------------------------------
    call read_file('FDIPS.in')
    call read_init
    if(iProc==0) call read_echo_set(.true.)

    ! Default decomposition
    nProcTheta = nProc
    nProcPhi   = 1

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#DOMAIN")
          call read_var('rMin', rMin)
          call read_var('rMax', rMax)
          call read_var('UseLogRadius', UseLogRadius)
       case("#GRID")
          call read_var('nR    ', nR)
          call read_var('nThetaAll', nThetaAll)
          call read_var('nPhiAll  ', nPhiAll)
       case("#PARALLEL")
          call read_var('nProcTheta', nProcTheta)
          call read_var('nProcPhi'  , nProcPhi)
       case("#MAGNETOGRAMFILE")
          IsNewMagnetogramStyle = .true.
          call read_var('NameFileIn' , NameFileIn)
          call read_var('BrMax'      , BrMax)
       case("#MAGNETOGRAM")
          call read_var('NameFileIn' , NameFileIn)
          call read_var('UseCosTheta', UseCosTheta)
          call read_var('BrMax'      , BrMax)
       case("#TIMING")
          call read_var('UseTiming', UseTiming)
       case("#TEST")
          call read_var('iProcTest', iProcTest)
          DoTestMe = iProc == iProcTest
       case("#TESTIJK")
          call read_var('iRTest'    , iRTest)
          call read_var('iPhiTest'  , iPhiTest)
          call read_var('iThetaTest', iThetaTest)
       case("#SOLVER")
          call read_var('NameSolver',         NameSolver, &
               IsUpperCase=.true.)
          call read_var('NamePreconditioner', NamePreconditioner, &
               IsUpperCase=.true.)
          call read_var('Tolerance',          Tolerance)
          UseHypre = index(NameSolver,'MG') > 0 .or. &
               index(NamePreconditioner,'MG') > 0
          UsePreconditioner = NameSolver == 'BICGSTAB' .and. &
               NamePreconditioner /= 'NONE'
       case("#HYPRE")
          call read_hypre_param
       case("#OUTPUT")
          call read_var('TypeOutput', TypeOutput, IsLowerCase=.true.)
          select case(TypeOutput)
          case('field')
             DoSaveField = .true.
             call read_var('NameFileField', NameFileField)
             call read_var('TypeFileField', TypeFileField)
             ! remove .out extension if present
             i = index(NameFileField,'.out')
             if(i>0) NameFileField = NameFileField(1:i-1)
          case('potential')
             DoSavePotential = .true.
             call read_var('NameFilePotential', NameFilePotential)
             call read_var('TypeFilePotential', TypeFilePotential)
             ! remove .out extension if present
             i = index(NameFilePotential,'.out')
             if(i>0) NameFilePotential = NameFilePotential(1:i-1)
          case('tecplot')
             if(nProc > 1)call CON_stop(NameSub// &
                  ': TypeOutput=tecplot works for serial runs only')
             DoSaveTecplot = .true.
             call read_var('NameFileTecplot', NameFileTecplot)
          case default
             call CON_stop(NameSub//': unknown TypeOutput='//trim(TypeOutput))
          end select
       case default
          call CON_stop(NameSub//': unknown command='//trim(NameCommand))
       end select
    end do

    if ( nProcTheta*nProcPhi /= nProc .and. iProc==0) then
       write(*,*)NameSub,': nProcTheta, nProcPhi, nProc=', &
            nProcTheta, nProcPhi, nProc
       call CON_stop(NameSub//': nProc should be nProcTheta*nProcPhi')
    end if

    ! Do timing on proc 0 only, if at all
    if(iProc > 0) UseTiming = .false.

  end subroutine read_fdips_param
  !===========================================================================
  subroutine read_magnetogram

    use ModIoUnit, ONLY: UnitTmp_
    use ModPlotFile, ONLY: read_plot_file

    ! Read the raw magnetogram file into a 2d array

    integer:: iError
    integer:: nCarringtonRotation
    integer:: nTheta0, nPhi0, nThetaRatio, nPhiRatio
    integer:: iTheta, iPhi, iTheta0, iTheta1, jPhi0, jPhi1, jPhi, kPhi
    real :: BrAverage, Weight
    character (len=100) :: String

    real, allocatable:: Br0_II(:,:), Var_II(:,:), Phi0_I(:), Theta0_I(:)
    real:: Param_I(1)

    character(len=*), parameter:: NameSub = 'read_magnetogram'
    !------------------------------------------------------------------------
    if(IsNewMagnetogramStyle)then

       call read_plot_file(NameFileIn, n1Out = nPhi0, n2Out = nTheta0, &
            ParamOut_I=Param_I, iErrorOut=iError)

       if(iError /= 0) call CON_stop(NameSub// &
            ': could not read header from file'//trim(NameFileIn))

       write(*,*)'nTheta0, nPhi0, LongitudeShift: ', nTheta0, nPhi0, Param_I

       allocate(Phi0_I(nPhi0), Theta0_I(nTheta0), Var_II(nPhi0,nTheta0), &
            Br0_II(nTheta0,nPhi0))
       
       call read_plot_file(NameFileIn, &
            Coord1Out_I=Phi0_I, Coord2Out_I=Theta0_I, VarOut_II = Var_II, &
            iErrorOut=iError)

       if(iError /= 0) call CON_stop(NameSub// &
            ': could not read date from file'//trim(NameFileIn))

       ! Check if the theta coordinate is uniform or not
       UseCosTheta = abs(Theta0_I(3) - 2*Theta0_I(2) + Theta0_I(1)) > 1e-6

       ! Br0 is defined with the opposite index order
       Br0_II = transpose(Var_II)

       deallocate(Var_II)

    else
       open(UnitTmp_, file=NameFileIn, status='old', iostat=iError)
       if(iError /= 0) call CON_stop(NameSub// &
            ': could not open input file'//trim(NameFileIn))
       do 
          read(UnitTmp_,'(a)', iostat = iError ) String
          if(index(String,'#CR')>0)then
             read(UnitTmp_,*) nCarringtonRotation
          endif
          if(index(String,'#ARRAYSIZE')>0)then
             read(UnitTmp_,*) nPhi0
             read(UnitTmp_,*) nTheta0
          endif
          if(index(String,'#START')>0) EXIT
       end do

       write(*,*)'nCarringtonRotation, nTheta0, nPhi0: ',&
            nCarringtonRotation, nTheta0, nPhi0

       allocate(Br0_II(nTheta0,nPhi0))

       ! input file is in longitude, latitude
       do iTheta = nTheta0, 1, -1
          do iPhi = 1, nPhi0
             read(UnitTmp_,*) Br0_II(iTheta,iPhi)
          end do
       end do
       close(UnitTmp_)
    end if

    ! Fix too large values of Br
    where (abs(Br0_II) > BrMax) Br0_II = sign(BrMax, Br0_II)

    if(nTheta0 > nThetaAll .and. nThetaAll > 1)then

       if(modulo(nTheta0, nThetaAll) /= 0)then
          write(*,*) NameSub,' nTheta in file    =', nTheta0
          write(*,*) NameSub,' nTheta in FDIPS.in=', nThetaAll
          call CON_stop(NameSub//': not an integer coarsening ratio')
       end if
       ! Set integer coarsening ratio
       nThetaRatio = nTheta0 / nThetaAll
       nThetaAll   = nTheta0 / nThetaRatio
    else
       nThetaRatio = 1
       nThetaAll   = nTheta0
    end if

    if(nPhi0 > nPhiAll .and. nPhiAll > 1)then
       if(modulo(nPhi0, nPhiAll) /= 0)then
          write(*,*) NameSub,' nPhi in file    =', nPhi0
          write(*,*) NameSub,' nPhi in FDIPS.in=', nPhiAll
          call CON_stop(NameSub//': not an integer coarsening ratio')
       end if
       nPhiRatio = nPhi0 / nPhiAll
       nPhiAll   = nPhi0 / nPhiRatio
    else
       nPhiRatio = 1
       nPhiAll   = nPhi0
    end if

    allocate(Br_II(nThetaAll,nPhiAll))
    Br_II = 0.0

    do iPhi = 1, nPhiAll
       jPhi0 = nPhiRatio*(iPhi-1) - nPhiRatio/2 + 1
       jPhi1 = nPhiRatio*(iPhi-1) + nPhiRatio/2 + 1

       do jPhi = jPhi0, jPhi1

          if( modulo(nPhiRatio,2) == 0 .and. &
               (jPhi == jPhi0 .or. jPhi == jPhi1) )then
             ! For even coarsening ratio use 0.5 weight at the two ends
             Weight = 0.5
          else
             Weight = 1.0
          end if

          ! Apply periodicity
          kPhi = modulo(jPhi-1,nPhi0) + 1

          do iTheta = 1, nThetaAll
             iTheta0 = nThetaRatio*(iTheta-1) + 1
             iTheta1 = iTheta0 + nThetaRatio - 1

             Br_II(iTheta,iPhi) = Br_II(iTheta,iPhi) &
                  + Weight * sum( Br0_II(iTheta0:iTheta1, kPhi))
          end do
       end do
    end do

    Br_II = Br_II / (nThetaRatio*nPhiRatio)

    ! remove monopole
    BrAverage = sum(Br_II)/(nThetaAll*nPhiAll)
    Br_II = Br_II - BrAverage

    deallocate(Br0_II)

  end subroutine read_magnetogram

  !============================================================================

  subroutine init_potential_field

    use ModConst, ONLY: cPi, cTwoPi

    integer :: iR, iTheta, iPhi
    real:: dR, dLogR, dTheta, dPhi, dZ, z
    !--------------------------------------------------------------------------

    ! The processor coordinate
    iProcTheta = iProc/nProcPhi
    iProcPhi   = iProc - iProcTheta*nProcPhi

    ! Calculate the nTheta, nPhi. To distribute as even as possible,
    ! there will be two different nTheta and nPhi. 
    nThetaLgr  = ceiling(real(nThetaAll)/nProcTheta)
    nThetaSml  = floor(  real(nThetaAll)/nProcTheta)
    nPhiLgr    = ceiling(real(nPhiAll)/nProcPhi)
    nPhiSml    = floor(  real(nPhiAll)/nProcPhi)

    ! Calculate the number of processors which has large/small 
    ! local number of Theta/Phi.
    nProcThetaLgr = mod(nThetaAll, nProcTheta)
    nProcThetaSml = nProcTheta - nProcThetaLgr
    nProcPhiLgr = mod(nPhiAll, nProcPhi)
    nProcPhiSml = nProcPhi - nProcPhiLgr

    ! Test if the partitioning works
    if (iProc == 0) then
       write(*,*) 'nThetaLgr = ',nThetaLgr, 'nThetaSml = ', nThetaSml
       write(*,*) 'nPhiLgr   = ', nPhiLgr,  'nPhiSml   = ', nPhiSml
       write(*,*) 'Partitioning in nThetaAll gives: ', nThetaLgr*nProcThetaLgr + &
            nThetaSml*nProcThetaSml, &
            'Actual nThetaAll is: ', nThetaAll
       write(*,*) 'Partitioning in nPhiAll gives:   ', nPhiLgr*nProcPhiLgr +  &
            nPhiSml*nProcPhiSml, &
            'Actual nPhiAll is:   ', nPhiAll
    end if

    !Both iProcTheta and iProcPhi in the large region
    if (iProcTheta < nProcThetaLgr .and. iProcPhi < nProcPhiLgr) then
       nTheta = nThetaLgr
       nPhi   = nPhiLgr
       iTheta0 = iProcTheta* nThetaLgr
       iPhi0   = iProcPhi  * nPhiLgr
    end if

    !Only iProcTheta in the large region
    if (iProcTheta < nProcThetaLgr .and. iProcPhi >= nProcPhiLgr) then
       nTheta = nThetaLgr
       nPhi   = nPhiSml
       iTheta0 = iProcTheta  * nThetaLgr
       iPhi0   = nProcPhiLgr * nPhiLgr + (iProcPhi - nProcPhiLgr)*nPhiSml
    end if

    !Only iProcPhi in the large region
    if (iProcTheta >= nProcThetaLgr .and. iProcPhi < nProcPhiLgr) then
       nTheta = nThetaSml
       nPhi   = nPhiLgr
       iTheta0 = nProcThetaLgr * nThetaLgr + (iProcTheta - nProcThetaLgr)*nThetaSml
       iPhi0   = iProcPhi      * nPhiLgr
    end if

    !Both iProcTheta and iProcPhi in the small region
    if (iProcTheta >= nProcThetaLgr .and. iProcPhi >= nProcPhiLgr) then
       nTheta = nThetaSml
       nPhi   = nPhiSml
       iTheta0 = nProcThetaLgr*nThetaLgr &
            + (iProcTheta - nProcThetaLgr)*nThetaSml
       iPhi0   = nProcPhiLgr  *nPhiLgr   &
            + (iProcPhi   - nProcPhiLgr)  *nPhiSml
    end if

    allocate( BrLocal_II(nTheta,nPhi), &
         Radius_I(0:nR+1), Theta_I(0:nTheta+1), Phi_I(0:nPhi+1), &
         dRadius_I(nR), dPhi_I(nPhi), &
         SinTheta_I(0:nTheta+1), dTheta_I(nTheta), dCosTheta_I(nTheta), &
         SinThetaNode_I(nTheta+1), dCosThetaNode_I(nTheta+1), &
         RadiusNode_I(nR+1), ThetaNode_I(nTheta+1), PhiNode_I(nPhi+1), &
         dRadiusNode_I(nR+1), dThetaNode_I(nTheta+1), dPhiNode_I(nPhi+1) , &
         Potential_C(nR,nTheta,nPhi), &
         Rhs_C(nR,nTheta,nPhi), &
         B0_DF(3,nR+1,nTheta+1,nPhi+1), &
         DivB_C(nR,nTheta,nPhi))

    ! Set BrLocal_II, this is used in set_boundary when UseBr is true
    BrLocal_II(:,:) = Br_II(iTheta0 + 1: iTheta0 + nTheta, &
         iPhi0   + 1: iPhi0   + nPhi)

    ! nR is the number of mesh cells in radial direction
    ! cell centered radial coordinate

    if(UseLogRadius)then
       dLogR = log(rMax/rMin)/nR
       do iR = 0, nR+1
          Radius_I(iR) = rMin*exp( (iR - 0.5)*dLogR )
       end do
       ! node based radial coordinate                                               
       do iR = 1, nR+1
          RadiusNode_I(iR) = rMin*exp( (iR - 1)*dLogR )
       end do

    else
       dR = (rMax - rMin)/nR
       do iR = 0, nR+1
          Radius_I(iR) = rMin + (iR - 0.5)*dR
       end do
       ! node based radial coordinate
       do iR = 1, nR+1
          RadiusNode_I(iR) = rMin + (iR - 1)*dR
       end do
    end if

    dRadius_I     = RadiusNode_I(2:nR+1) - RadiusNode_I(1:nR)
    dRadiusNode_I = Radius_I(1:nR+1) - Radius_I(0:nR)

    if(UseCosTheta)then
       dZ = 2.0/nThetaAll

       !Set Theta_I
       do iTheta = 0, nTheta+1
          z = max(-1.0, min(1.0, 1 - (iTheta + iTheta0 - 0.5)*dZ))
          Theta_I(iTheta) = acos(z)
       end do

       !Set the boundary condition of Theta_I
       if (iProcTheta == 0) &
            Theta_I(0) = -Theta_I(1)
       if (iProcTheta == nProcTheta-1) &
            Theta_I(nTheta+1) = cTwoPi - Theta_I(nTheta)

       !Set ThetaNode_I
       do iTheta = 1, nTheta + 1
          z = max(-1.0, min(1.0, 1 - (iTheta + iTheta0 -1)*dZ))
          ThetaNode_I(iTheta) = acos(z)
       end do
    else

       dTheta = cPi/nThetaAll

       !Set Theta_I
       do iTheta = 0, nTheta+1
          Theta_I(iTheta) = (iTheta  + iTheta0 - 0.5)*dTheta
       end do

       !Set ThetaNode_I
       do iTheta = 1, nTheta+1
          ThetaNode_I(iTheta) = (iTheta + iTheta0 - 1)*dTheta
       end do
    end if

    dTheta_I = ThetaNode_I(2:nTheta+1) - ThetaNode_I(1:nTheta)
    SinTheta_I = sin(Theta_I)
    SinThetaNode_I = sin(ThetaNode_I)

    if(UseCosTheta)then
       dCosTheta_I     = dZ
       dCosThetaNode_I = dZ
       dThetaNode_I    = Theta_I(1:nTheta+1) - Theta_I(0:nTheta)
    else
       dCosTheta_I(1:nTheta) = SinTheta_I(1:nTheta)*dTheta
       dCosThetaNode_I       = SinThetaNode_I*dTheta
       dThetaNode_I          = dTheta
    end if

    dPhi = cTwoPi/nPhiAll
    !Set Phi_I
    do iPhi = 0, nPhi+1
       Phi_I(iPhi) = (iPhi + iPhi0 - 1)*dPhi
    end do

    PhiNode_I = Phi_I(1:nPhi+1) - 0.5*dPhi
    dPhi_I = PhiNode_I(2:nPhi+1) - PhiNode_I(1:nPhi)
    dPhiNode_I = Phi_I(1:nPhi+1) - Phi_I(0:nPhi)

    Potential_C       =   0.0
    Rhs_C             =   0.0

  end subroutine init_potential_field

  !============================================================================

  subroutine save_potential_field

    use ModIoUnit,      ONLY: UnitTmp_
    use ModNumConst,    ONLY: cHalfPi
    use ModPlotFile,    ONLY: save_plot_file

    integer :: iR, jR, iTheta, iPhi
    real    :: r, CosTheta, SinTheta, CosPhi, SinPhi
    real    :: Br, Btheta, Bphi
    real    :: rI, rJ, rInv
    real, allocatable :: B_DX(:,:,:,:), B_DII(:,:,:)
    integer :: iError
    integer :: iStatus_I(mpi_status_size)
    integer:: nPhiOut
    !-------------------------------------------------------------------------

    ! Only the last processors in the phi direction write out the ghost cell
    if(iProcPhi == nProcPhi - 1)then
       nPhiOut = nPhi + 1
    else
       nPhiOut = nPhi
    end if

    allocate(B_DX(3,nR+1,nPhiOut,nTheta), B_DII(3,nR+1,nTheta))

    ! Average the magnetic field to the R face centers

    ! For the radial component only the theta index changes
    do iPhi = 1, nPhi; do iTheta = 1, nTheta; do iR = 1, nR+1
       B_DX(1,iR,iPhi,nTheta+1-iTheta) = B0_DF(1,iR,iTheta,iPhi)
    end do; end do; end do

    ! Use radius as weights to average Bphi and Btheta 
    ! Also swap phi and theta components (2,3) -> (3,2)
    do iPhi = 1, nPhi; do iTheta = 1, nTheta; do iR = 1, nR

       ! Use first order approximation at lower boundary. Reduces noise.
       jR = max(1,iR-1)

       rI   = Radius_I(iR)
       rJ   = Radius_I(jR)
       rInv = 0.25/RadiusNode_I(iR)

       B_DX(2,iR,iPhi,nTheta+1-iTheta) = rInv* &
            ( rI*(B0_DF(3,iR,iTheta,iPhi) + B0_DF(3,iR,iTheta,iPhi+1)) &
            + rJ*(B0_DF(3,jR,iTheta,iPhi) + B0_DF(3,jR,iTheta,iPhi+1)) )

       B_DX(3,iR,iPhi,nTheta+1-iTheta) = rInv* &
            ( rI*(B0_DF(2,iR,iTheta,iPhi) + B0_DF(2,iR,iTheta+1,iPhi)) &
            + rJ*(B0_DF(2,jR,iTheta,iPhi) + B0_DF(2,jR,iTheta+1,iPhi)) )

    end do; end do; end do

    ! set tangential components to zero at the top
    B_DX(2:3,nR+1,:,:) = 0.0
    
    ! Apply periodicity in Phi to fill the nPhi+1 ghost cell
    if (nProcPhi > 1) then
       if (iProcPhi ==0 ) then 
          b_DII = B_DX(:,:,1,:)
          call mpi_send(b_DII, 3*(nR+1)*nTheta, MPI_REAL, &
               iProcTheta*nProcPhi + nProcPhi-1, 21, iComm,  iError)
       end if
       if (iProcPhi == nProcPhi -1) then
          call mpi_recv(b_DII, 3*(nR+1)*nTheta, MPI_REAL, &
               iProcTheta*nProcPhi , 21, iComm, iStatus_I, iError)
          B_DX(:,:,nPhiOut,:) = b_DII
       end if
    else
       B_DX(:,:,nPhiOut,:) = B_DX(:,:,1,:)
    end if

    if(DoSaveField)then
       ! Note the fake processor index to be used by redistribute.pl
       write(NameFile,'(a,2i2.2,a,i3.3,a)') &
            trim(NameFileField)//'_np01', nProcPhi, nProcTheta, '_', &
            iProcPhi + (nProcTheta - 1 - iProcTheta)*nProcPhi, '.out'

       call save_plot_file(NameFile, TypeFileIn=TypeFileField, &
            StringHeaderIn = &
            'Radius [Rs] Longitude [Rad] Latitude [Rad] B [G]', &
            nameVarIn = 'Radius Longitude Latitude Br Bphi Btheta' &
            //' Ro_PFSSM Rs_PFSSM', &
            ParamIn_I = (/ rMin, rMax /), &
            nDimIn=3, VarIn_VIII=B_DX, &
            Coord1In_I=RadiusNode_I, &
            Coord2In_I=Phi_I(1:nPhiOut), &
            Coord3In_I=cHalfPi-Theta_I(nTheta:1:-1))
    end if

    if(DoSaveTecplot)then
       open(unit = UnitTmp_, file=NameFileTecplot, status='replace')

       write (UnitTmp_, '(a)') 'Title = "'     // 'PFSSM' // '"'
       write (UnitTmp_, '(a)') &
         'Variables = ' // trim ('"X [Rs]", "Y [Rs]", "Z [Rs]","Bx [G]",'// &
	' "By [G]", "Bz [G]"')
       write(UnitTmp_, '(a)') 'ZONE T="Rectangular zone"'
       write(UnitTmp_, '(a,i6,a,i6,a,i6,a)') &
            ' I = ', nR+1, ', J=', nThetaAll, ', K=', nPhiAll+1, ', ZONETYPE=Ordered'
       write(UnitTmp_, '(a)')' DATAPACKING=POINT'
       write(UnitTmp_, '(a)')' DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )'

       do iPhi = 1, nPhiAll+1; do iTheta = 1, nThetaAll; do iR = 1, nR+1
          Br     = B_DX(1,iR,iPhi,nThetaAll+1-iTheta)
          Btheta = B_DX(3,iR,iPhi,nThetaAll+1-iTheta)
          Bphi   = B_DX(2,iR,iPhi,nThetaAll+1-iTheta)
          r = RadiusNode_I(iR)
          SinTheta = SinTheta_I(iTheta)
          CosTheta = cos(Theta_I(iTheta))
          SinPhi   = sin(Phi_I(iPhi))
          CosPhi   = cos(Phi_I(iPhi))

          write (UnitTmp_,fmt="(6(E14.6))") &
               r*SinTheta*CosPhi, &
               r*SinTheta*SinPhi, &
               r*CosTheta, &
               Br*SinTheta*CosPhi + Btheta*CosTheta*CosPhi - Bphi*SinPhi, &
               Br*SinTheta*SinPhi + Btheta*CosTheta*SinPhi + Bphi*CosPhi, &
               Br*CosTheta        - Btheta*SinTheta
       end do; end do; end do

       close(UnitTmp_)
    end if

    deallocate(B_DX)

  end subroutine save_potential_field

  !============================================================================

  subroutine set_boundary(x_C, x_G)

    real, intent(in):: x_C(nR,nTheta,nPhi)
    real, intent(inout):: x_G(0:nR+1,0:nTheta+1,0:nPhi+1)

    integer:: iPhi, jPhi, shift
    integer:: status(mpi_status_size)
    integer:: SendRequest010, SendRequest020, SendRequest360, SendRequest180, &
         SendRequest12, SendRequest21, SendRequest34, SendRequest43, &
         RecvRequest010, RecvRequest020, RecvRequest360, RecvRequest180, &
         RecvRequest12, RecvRequest21, RecvRequest34, RecvRequest43
    integer:: jProc
    integer:: iError

    if (.not. allocated(tmpXPhi0_II)) allocate( &
         tmpXPhi0_II(0:nR+1,nPhiAll),              &
         tmpXPhipi_II(0:nR+1,nPhiAll),             &
         sendBC010_II(0:nR+1,nPhi),        &
         sendBC180_II(0:nR+1,nPhi),        &
         sendBC12_II(0:nR+1,0:nTheta+1),   &
         sendBC21_II(0:nR+1,0:nTheta+1),   &
         sendBC34_II(0:nR+1,nPhi),         &
         sendBC43_II(0:nR+1,nPhi),         &
         sendBC020_II(0:nR+1,0:nTheta+1),  &
         sendBC360_II(0:nR+1,0:nTheta+1),  &
         recvBCLgr010_II(0:nR+1,nPhiLgr),       &
         recvBCLgr180_II(0:nR+1,nPhiLgr),       &
         recvBCSml010_II(0:nR+1,nPhiSml),       &
         recvBCSml180_II(0:nR+1,nPhiSml),       &
         recvBC12_II(0:nR+1,0:nTheta+1),   &
         recvBC21_II(0:nR+1,0:nTheta+1),   &
         recvBC34_II(0:nR+1,nPhi),         &
         recvBC43_II(0:nR+1,nPhi),         &
         recvBC020_II(0:nR+1,0:nTheta+1),  &
         recvBC360_II(0:nR+1,0:nTheta+1))

    !--------------------------------------------------------------------------
    ! Current solution inside
    x_G(1:nR,1:nTheta,1:nPhi) = x_C

    ! The slope is forced to be Br at the inner boundary
    if(UseBr)then
       x_G(0,1:nTheta,1:nPhi) = x_C(1,:,:) - dRadiusNode_I(1)*BrLocal_II
    else
       x_G(0,1:nTheta,1:nPhi) = x_C(1,:,:)
    end if

    ! The potential is zero at the outer boundary
    x_G(nR+1,:,:) = -x_G(nR,:,:)

    ! Set tmpXPhi0_II and tmpXPhipi_II which used to store the global 
    ! boundary close to the pole to the root
    if (iProcTheta ==0) then
       sendBC010_II = x_G(:,1,1:nPhi)
       call MPI_ISEND(sendBC010_II, (nR+2)*nPhi,MPI_REAL, 0, 010, &
                      iComm, SendRequest010, iError)
    end if
    if (iProcTheta == nProcTheta-1) then
       sendBC180_II = x_G(:,nTheta,1:nPhi)
       call MPI_ISEND(sendBC180_II, (nR+2)*nPhi,MPI_REAL, 0, 180, &
                      iComm, SendRequest180, iError)
    end if

    ! The root update tmpXPhi0_II/tmpXPhipi_II from all processors
    if (iProc == 0) then
       do jProc=0, nProcPhi-1
          if (jProc < nProcPhiLgr) then
             shift = jProc * nPhiLgr
             
             call MPI_IRECV(recvBCLgr010_II, (nR+2)*nPhiLgr, MPI_REAL, &
                            jProc, &
                            010, iComm, RecvRequest010, iError)
             call mpi_wait(RecvRequest010, status, iError)
             tmpXPhi0_II(:, shift + 1: shift + nPhiLgr) = recvBCLgr010_II
                
             call MPI_IRECV(recvBCLgr180_II, (nR+2)*nPhiLgr, MPI_REAL, &
                            (nProcTheta-1)*nProcPhi+jProc, &
                            180, iComm, RecvRequest180, iError)
             call mpi_wait(RecvRequest180, status, iError)
             tmpXPhipi_II(:, shift + 1: shift + nPhiLgr) = recvBCLgr180_II
          else
             shift = nProcPhiLgr*nPhiLgr + (jProc - nProcPhiLgr)*nPhiSml

             call MPI_IRECV(recvBCSml010_II, (nR+2)*nPhiSml,  MPI_REAL, &
                            jProc, &
                            010, iComm, RecvRequest010, iError)
             call mpi_wait(RecvRequest010, status, iError)
             tmpXPhi0_II(:, shift + 1: shift + nPhiSml) = recvBCSml010_II

             call MPI_IRECV(recvBCSml180_II , (nR+2)*nPhiSml,  MPI_REAL, &
                            (nProcTheta-1)*nProcPhi+jProc, &
                            180, iComm, RecvRequest180, iError)
             call mpi_wait(RecvRequest180, status, iError)
             tmpXPhipi_II(:, shift + 1: shift + nPhiSml) = recvBCSml180_II
          end if
       end do
    end if

    call  MPI_bcast(tmpXPhi0_II,  (nR+2)*nPhiAll, MPI_REAL, 0,  iComm, iError)
    call  MPI_bcast(tmpXPhipi_II, (nR+2)*nPhiAll, MPI_REAL, 0,  iComm, iError)

    ! Symmetric in Theta but shifted by nPhiAll/2, be careful about the shift!
    if (iProcTheta == 0) then
       do iPhi = 1, nPhi
          jPhi = modulo(iPhi + iPhi0 - 1 + nPhiAll/2, nPhiAll) + 1
          x_G(:,0,iPhi)        = tmpXPhi0_II(:,jPhi)
       end do
    end if
    if (iProcTheta == nProcTheta-1) then
       do iPhi = 1, nPhi
          jPhi = modulo(iPhi + iPhi0 - 1 + nPhiAll/2, nPhiAll) + 1
          x_G(:,nTheta+1,iPhi) = tmpXPhipi_II(:,jPhi)
       end do
    end if

    !Update the local theta boundary
    if (iProcTheta /= nProcTheta-1) then
       sendBC34_II = x_G(:,nTheta,1:nPhi)
       call MPI_ISEND(sendBC34_II, (nR+2)*nPhi, MPI_REAL, &
                     (iProcTheta+1)*nProcPhi+iProcPhi, &
                     34, iComm, SendRequest34,  iError)
    end if
    if (iProcTheta /= 0) then
       sendBC43_II = x_G(:,1,1:nPhi)
       call MPI_ISEND(sendBC43_II, (nR+2)*nPhi, MPI_REAL, &
                     (iProcTheta-1)*nProcPhi+iProcPhi, &
                     43, iComm,  SendRequest43,  iError)
    end if
    if (iProcTheta /= nProcTheta-1) then
       call MPI_IRECV(recvBC43_II, (nR+2)*nPhi, MPI_REAL, &
                     (iProcTheta+1)*nProcPhi+iProcPhi, &
                     43, iComm,  RecvRequest43, iError)
       call mpi_wait(RecvRequest43, status, iError)
       x_G(:,nTheta+1,1:nPhi) = recvBC43_II
    end if
    if (iProcTheta /= 0) then
       call MPI_IRECV(recvBC34_II, (nR+2)*nPhi, MPI_REAL, &
                     (iProcTheta-1)*nProcPhi+iProcPhi, &
                     34, iComm,  RecvRequest34, iError)
       call mpi_wait(RecvRequest34, status, iError)
       x_G(:,0,1:nPhi) = recvBC34_II
    end if

    ! Periodic around phi
    ! Send boundary info
    if (iProcPhi == 0) then
       sendBC020_II = x_G(:,:,1)
       call MPI_ISEND(sendBC020_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + nProcPhi-1, &
                      020, iComm, SendRequest020, iError)
    end if
    if (iProcPhi == nProcPhi-1) then
       sendBC360_II = x_G(:,:,nPhi)
       call MPI_ISEND(sendBC360_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi ,  &
                      360, iComm, SendRequest360, iError)
    end if

    ! Update boundary info
    if (iProcPhi == 0) then
       call MPI_IRECV(recvBC360_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + nProcPhi-1, &
                      360, iComm, RecvRequest360, iError)
       call mpi_wait(RecvRequest360, status, iError)
       x_G(:,:,0) = recvBC360_II
    end if
    if (iProcPhi == nProcPhi-1) then
       call MPI_IRECV(recvBC020_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi , &
                      020, iComm, RecvRequest020, iError)
       call mpi_wait(RecvRequest020, status, iError)
       x_G(:,:,nPhi+1) = recvBC020_II
    end if

    ! Start to send and update the local boundary
    if (iProcPhi /= nProcPhi-1) then
       sendBC12_II = x_G(:,:,nPhi)
       call MPI_ISEND(sendBC12_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + iProcPhi+1, &
                      12,  iComm, SendRequest12, iError)
    end if
    if (iProcPhi /= 0) then
       sendBC21_II = x_G(:,:,1)
       call MPI_ISEND(sendBC21_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + iProcPhi-1, &
                      21,  iComm, SendRequest21, iError)
    end if
    if (iProcPhi /= nProcPhi-1) then
       call MPI_IRECV(recvBC21_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + iProcPhi+1, &
                      21,  iComm, RecvRequest21, iError)
       call mpi_wait(RecvRequest21, status, iError)
       x_G(:,:,nPhi+1) = recvBC21_II
    end if
    if (iProcPhi /= 0) then
       call MPI_IRECV(recvBC12_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + iProcPhi-1, &
                      12,  iComm, RecvRequest12, iError)
       call mpi_wait(RecvRequest12, status, iError)
       x_G(:,:,0) = recvBC12_II
    end if

  end subroutine set_boundary

end module ModPotentialField
