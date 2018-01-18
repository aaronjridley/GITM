!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
program post_idl

  ! Read a header file from STDIN, then read the data from the .idl files 
  ! written by the multiple processors.
  ! Set the coordinate and plot variable arrays (Coord_DC and PlotVar_VC).
  ! For structured output the data is put into a regular grid.
  ! For unstructured output only the first index is used (linear arrays).
  ! Average out points that have the same coordinates.
  ! Save the results into an output file.

  use ModPlotFile,       ONLY: save_plot_file
  use ModNumConst,       ONLY: cHalfPi
  use ModInterpolate,    ONLY: linear
  use ModNumConst,       ONLY: cTwoPi, cRadToDeg
  use ModCoordTransform, ONLY: rot_matrix_z
  use ModUtilities,      ONLY: lower_case, split_string, join_string
  use ModSort,           ONLY: sort_quick
  use ModReadParam,      ONLY: read_file, read_init, read_echo_set, &
       read_line, read_command, read_var, lStringLine
  use ModIoUnit,         ONLY: UnitTmp_
  use ModMpi,            ONLY: MPI_COMM_SELF

  implicit none

  ! This is copied from ModKind, because PostIDL.exe may be compiled with
  ! different precision then the rest of the codes.
  integer, parameter :: Real4_=selected_real_kind(6,30)
  integer, parameter :: Real8_=selected_real_kind(12,100)
  integer, parameter :: nByteReal = 4 + (1.00000000041 - 1.0)*10000000000.0

  ! Global variables

  character(len=20) :: TypeFile='real4'

  integer:: nCell_D(3), n1, n2, n3
  integer:: iCell, nCellCheck, nCellPlot, nPlotVar, nParamPlot, nProc, nStep
  real :: t
  real, allocatable :: Coord_DC(:,:,:,:), PlotVar_VC(:,:,:,:)
  real, allocatable :: GenCoord_DI(:,:) ! gen. coords for unstructured grids
  real, allocatable :: PlotVar_V(:), PlotParam_I(:)
  real, allocatable :: Param_I(:)

  real(Real4_)              :: DxCell4, Xyz4_D(3)
  real(Real4_), allocatable :: PlotVar4_V(:)
  real(Real8_)              :: DxCell8, Xyz8_D(3)
  real(Real8_), allocatable :: PlotVar8_V(:)

  ! Coordinates, sizes, indices
  real, dimension(3) :: Xyz_D, CoordMin_D, CoordMax_D
  real, dimension(3) :: CellSizePlot_D, dCoordPlot_D, dCoordMin_D
  real, dimension(3) :: GenCoord_D
  real::    CellSize_D(3), CellSize1, dCoord1Plot
  real::    Fraction
  real(selected_real_kind(12))  :: VolumeCheck, VolumePlot
  integer:: Ijk_D(3), i, j, k, iVar
  integer:: IjkMin_D(3), IjkMax_D(3), iMin, iMax, jMin, jMax, kMin, kMax

  ! Variables related to sorting and averaging unstructured data
  logical:: DoSort3D = .false.           ! do sort unstructured 3D files?
  integer, allocatable :: iSort_I(:)     ! indirect index for sorting
  real, allocatable :: Sort_I(:)         ! sorting criterion
  integer:: nSum                         ! number of points averaged
  real, allocatable :: StateSum_V(:)     ! sum of states
  real, allocatable :: CellSizeMin_D(:)  ! min distance among different points

  integer :: iDim, iDimCut_D(3), nDim, nParamExtra, l1, l2, l3, l4
  real    :: ParamExtra_I(3)
  character(len=5):: NameCoord_D(3) = (/'x    ','y    ','z    '/)
  character(len=5):: NameCoordPlot_D(3)

  logical :: IsStructured, DoReadBinary=.false.
  character (len=100) :: NameFile, NameFileHead, NameCoord
  character (len=lStringLine) :: NameVar, NameUnit
  integer :: l, iProc

  ! Variables for the 2D lookup table
  integer :: iDim1, iDim2
  integer :: iDim0 ! the ignored dimension
  integer :: iError
  real    :: CoordShift1, CoordShift2

  ! Variables for checking binary compatibility
  integer            :: nByteRealRead

  ! Variables for generalized coordinates
  character (len=79) :: TypeGeometry='cartesian'
  logical            :: UseDoubleCut = .false.

  ! Logarithmic radial coordinate
  logical:: IsLogRadius = .false.

  ! Stretched radial coordinates
  logical:: IsGenRadius = .false.
  integer:: nRgen = -1
  real, allocatable:: LogRgen_I(:)

  ! Roundcube parameters
  real:: rRound0    ! fully Cartesian distance
  real:: rRound1    ! fully round distance
  real:: SqrtNDim   ! sqrt 2 or sqrt 3

  integer :: nDimSim ! Dimension of simulation domain.

  ! Periodicity of the simulation domain. Not needed here
  logical:: IsPeriodic_D(3)

  ! Control verbose output
  logical:: IsVerbose = .false.

  ! Tecplot related variables for reading tecplot .dat to save as idl file
  logical :: DoReadTecplot = .false.
  integer :: nHeaderTec=0, nNodeTec=0, nCellTec=0
  integer :: iHeaderTec, iNodeTec, iCellTec, iNeiTec
  character(len=500) :: StringTecHeader
  real :: xMinTec, xMaxTec
  real, allocatable :: XyzTec_DI(:,:), PlotVarTec_VI(:,:)
  integer, allocatable:: iNodeTec_II(:,:)
  
  character (len=lStringLine) :: NameCommand, StringLine
  !---------------------------------------------------------------------------

  write(*,'(a)')'PostIDL (G.Toth 2000-) starting'

  ! No file name given, so read from the header file from STDIN
  ! Use MPI_COMM_SELF to indicate serial execution
  call read_file(iCommIn = MPI_COMM_SELF)

  ! Initialize ModReadParam
  call read_init 
  
  call read_echo_set(.true.)

  ! Read the information from the header file
  READPARAM: do
     if(.not.read_line(StringLine) )then
        EXIT READPARAM
     end if

     if(.not.read_command(NameCommand)) CYCLE READPARAM

     select case(NameCommand)
     case('#VERBOSE')
        call read_var('IsVerbose', IsVerbose)

     case('#SORT3D')
        call read_var('DoSort3D', DoSort3D)

     case('#HEADFILE')
        call read_var('HeadNameFile', NameFileHead)
        call read_var('nProc',nProc)
        call read_var('SaveBinary', DoReadBinary)
        if(DoReadBinary) then
           call read_var('nByteReal', nByteRealRead)
           if(nByteRealRead==nByteReal)then
              write(*,*)'nByteReal=',nByteReal
           else if(nByteRealRead < nByteReal)then
              write(*,*)'!!! Warning: PostIDL was compiled with ',&
                   nByteReal,' byte reals but file contains nByteReal=', &
                   nByteRealRead
           end if
        endif
        
        ! Get rid of the directory part
        NameFileHead = NameFileHead( &
             index(NameFileHead,'/',BACK=.true.)+1:len(NameFileHead))

     case('#NDIM')
        call read_var('nDimSim', nDimSim)

     case('#NSTEP')
        call read_var('nStep', nStep)

     case('#TIMESIMULATION')
        call read_var('TimeSimulation', t)

     case('#PLOTRANGE')
        CoordMin_D = 0; CoordMax_D = 0        
        do i = 1, nDimSim
           call read_var('CoordMin', CoordMin_D(i))
           call read_var('CoordMax', CoordMax_D(i))
        enddo

     case('#PLOTRESOLUTION')
        CellSizePlot_D = 1
        do i = 1, nDimSim
           call read_var('CellSizePlot', CellSizePlot_D(i))
        enddo

     case('#CELLSIZE')
        dCoordMin_D = 1
        do i = 1, nDimSim
           call read_var('CellSizeMin', dCoordMin_D(i))
        enddo

     case('#NCELL')
        call read_var('nCellPlot', nCellPlot)

     case('#PLOTVARIABLE')
        call read_var('nPlotVar', nPlotVar)
        call read_var('NameVar',  NameVar)
        call read_var('NameUnit', NameUnit)

     case('#SCALARPARAM')
        call read_var('nParam', nParamPlot)
        allocate(PlotParam_I(nParamPlot))
        do i = 1, nParamPlot
           call read_var('Param',PlotParam_I(i))
        enddo

     case('#GRIDGEOMETRYLIMIT')
        call read_var('TypeGeometry', TypeGeometry)
        IsLogRadius = index(TypeGeometry,'lnr')  > 0
        IsGenRadius = index(TypeGeometry,'genr') > 0
        if(IsGenRadius)then
           call read_var('nRgen', nRGen)
           allocate(LogRgen_I(nRgen))
           do i = 1, nRgen
              call read_var('LogRgen', LogRgen_I(i))
           end do
        end if

        if(TypeGeometry == 'roundcube')then
           call read_var('rRound0', rRound0)
           call read_var('rRound1', rRound1)
           call read_var('SqrtNDim',SqrtNDim)
        endif

     case('#PERIODIC')
        do i = 1, nDimSim
           call read_var('IsPeriodic',IsPeriodic_D(i))
        enddo
        
     case('#OUTPUTFORMAT')
        call read_var('TypeOutPutFormat', TypeFile)

     case('#ROOTBLOCK', '#GRIDBLOCKSIZE')
        ! Ignore these commands
        
     case('#TECPLOTCONVERT')
        call read_var('DoReadTecplot', DoReadTecplot)
        call read_var('nHeaderTec', nHeaderTec)
        call read_var('nNodeTec',  nNodeTec)
        call read_var('nCellTec',  nCellTec)

     case default
        write(*,*) 'WARNING: unknown command ', NameCommand
     end select

  enddo READPARAM
  
  if(NameFileHead(1:3) == 'sph')then
     NameCoord_D    = (/'r    ','theta','phi  '/)
  elseif(TypeGeometry == 'rz' .or. TypeGeometry == 'xr')then
     TypeGeometry = 'cartesian'
     NameCoord_D    = (/'x    ','r    ','phi  '/)
  elseif(NameFileHead(1:3) == 'cut')then
     select case(TypeGeometry(1:5))
     case('spher')
        ! In reality this is the r-lon-lat coordinate system
        NameCoord_D = (/'r    ','lon  ','lat  '/)
     case('cylin')
        NameCoord_D = (/'r    ','phi  ','z    '/)
     case('round')
        NameCoord_D = (/'xi   ','eta  ','zeta '/)
     end select
  end if

  ! Save input CellSizePlot_D into dCoordPlot_D that may get overwritten
  dCoordPlot_D = CellSizePlot_D

  ! First component is used a lot, save into a scalar
  dCoord1Plot  = dCoordPlot_D(1)

  ! Unstructured grid is indicated with negative plot resolution
  IsStructured = dCoord1Plot >= 0.0

  ! If negative or zero plot resolution is set, 
  ! use the smallest grid cell size
  if(dCoord1Plot <= 0.0) dCoordPlot_D = dCoordMin_D

  ! Calculate structured grid size
  nCell_D = max(1, nint((CoordMax_D - CoordMin_D)/dCoordPlot_D))

  write(*,*)'plot area size=', nCell_D

  ! Calculate dimensionality of the cut and add parameters if needed
  nDim=0
  nParamExtra=0
  iDimCut_D=0
  do iDim = 1, 3
     if(nCell_D(iDim) > 1)then
        nDim = nDim + 1
        iDimCut_D(nDim) = iDim
     else
        iDimCut_D(3) = iDim
        nParamExtra = nParamExtra + 1
        ParamExtra_I(nParamExtra) = 0.5*(CoordMax_D(iDim) + CoordMin_D(iDim))
        NameVar=trim(NameVar)//' cut'//trim(NameCoord_D(iDim))
     end if
  end do

  if(IsVerbose)then
     write(*,*) 'dCoordPlot_D=', dCoordPlot_D
     write(*,*) 'IsStructured=', IsStructured
     write(*,*) 'NameCoord_D =', NameCoord_D
     write(*,*) 'nDim        =', nDim
     write(*,*) 'iDimCut_D   =', iDimCut_D
  end if

  if(nDim==2 .and. NameFileHead(1:3) /= 'cut')then
     ! For double cut 
     iDim1 = iDimCut_D(1)
     iDim2 = iDimCut_D(2)
     iDim0 = iDimCut_D(3)
     CoordShift1 = CoordMax_D(iDim1) - CoordMin_D(iDim1)
     CoordShift2 = CoordMax_D(iDim2) - CoordMin_D(iDim2)

     ! Sph/cyl. X=0 and Y=0 cuts require doubling the plot size (+/- r)
     if(iDim0 == 2)then
        if(TypeGeometry(1:3)=='sph' .and. CoordMax_D(iDim2) > cHalfPi) then
           ! Use LatMin < Lat' < 2*LatMax-LatMin as generalized coordinate
           UseDoubleCut = .true.; 
           ! nCell_D(3) = 2*nCell_D(3)
        elseif(TypeGeometry(1:3)=='cyl' .and. CoordMax_D(iDim2) > cHalfPi)then
           ! Use rMin < r' < 2*rMax - rMin as generalized coordinate
           UseDoubleCut = .true.;
           ! nCell_D(1) = 2*nCell_D(1)
        end if

        ! Do not attempt to use structured grid for double cuts
        if(UseDoubleCut)then
           IsStructured = .false.
           CellSizePlot_D = -1.0
        end if
     end if

     if(IsVerbose)then
        write(*,*) 'iDim1, iDim2, iDim0=', iDim1, iDim2, iDim0
        write(*,*) 'CoordShift1, Shift2=', CoordShift1, CoordShift2
        write(*,*) 'UseDoubleCut       =', UseDoubleCut
     end if

  endif

  ! For unstructured grid make the Coord_DC and PlotVar_VC arrays linear
  if(.not.IsStructured)then
     nCell_D(1)=nCellPlot
     nCell_D(2:3)=1
  end if

  ! Set scalars for plot size
  n1 = nCell_D(1);   n2 = nCell_D(2);   n3 = nCell_D(3)

  ! Allocate PlotVar_VC and Coord_DC, the arrays of variables and coordinates
  allocate( &
       PlotVar_V(nPlotVar), &
       PlotVar_VC(nPlotVar,n1,n2,n3), &
       Coord_DC(nDim,n1,n2,n3), &
       STAT=iError)
  if(iError /= 0) stop 'PostIDL.exe ERROR: could not allocate arrays'

  if(.not.IsStructured)then     
     allocate(GenCoord_DI(nDim,nCellPlot), STAT=iError)
     if(iError /= 0) stop 'PostIDL.exe ERROR: could not allocate enCoord_DI'
  end if

  if(DoReadBinary .and. nByteRealRead==4) allocate(PlotVar4_V(nPlotVar))
  if(DoReadBinary .and. nByteRealRead==8) allocate(PlotVar8_V(nPlotVar))
  
  ! Initialize PlotVar_VC and Coord_DC
  PlotVar_VC = 0.0
  Coord_DC = 0.0

  ! Reorder coordinate names for cuts
  NameCoordPlot_D(1:nDim) = NameCoord_D(iDimCut_D(1:nDim))

  ! Fix first coordinate name for non-cartesian cut along phi=90,270 deg
  ! that is a cut in the SECOND generalized coordinate.
  if(NameFileHead(1:3) == 'x=0') NameCoordPlot_D(1) = 'y'

  ! Create complete space separated list of "variable" names
  call join_string(NameCoordPlot_D(1:nDim), NameCoord)
  NameVar = trim(NameCoord)//' '//trim(NameVar)

  ! For Tecplot format NameUnit contains coordinate and variable names
  ! The coordinate names have to be fixed (except for 3D cartesian)

  if(NameUnit(1:11) == 'VARIABLES =')then

     l1 = index(NameUnit, '[')  ! start of dist unit for "X"
     l2 = index(NameUnit, ']')  ! end of dist unit for "X"
     l3 = index(NameUnit, '"Z') ! last coordinate is "Z"

     if(l3 > 11)then
        l4 = index(NameUnit(l3:len(NameUnit)), ',') + l3 ! end of coord names

        ! Reconstruct coordinate names in Tecplot format
        NameCoord = 'VARIABLES ="'
        do iDim = 1, nDim
           NameCoord = trim(NameCoord)//trim(NameCoordPlot_D(iDim))
           select case(NameCoordPlot_D(iDim))
           case('x', 'y', 'z', 'r')
              if(l1 > 1 .and. l2 > l1) &
                   NameCoord = trim(NameCoord)//NameUnit(l1-1:l2)
           case('theta', 'phi', 'lon', 'lat')
              NameCoord = trim(NameCoord)//' [deg]'
           end select
           NameCoord = trim(NameCoord)//'", "'
        end do

        NameUnit=NameCoord(1:len_trim(NameCoord)-2)//NameUnit(l4:len(NameUnit))
     end if
  end if

  ! Logic for Tecplot conversion
  if(DoReadTecplot)then
     l=len_trim(NameFileHead)
     write(NameFile,'(a)')NameFileHead(1:l-2)//".dat"
     write(*,*)'reading file=',trim(NameFile),' ...'

     allocate( &
          XyzTec_DI(3,nNodeTec), &
          PlotVarTec_VI(nPlotVar,nNodeTec), &
          iNodeTec_II(8,nCellTec), &
          STAT=iError)
     if(iError /= 0) stop 'PostIDL.exe ERROR: could not allocate arrays'

     open(UnitTmp_, file=NameFile, status='old', iostat=iError)

     ! Read Tecplot file header
     do iHeaderTec = 1, nHeaderTec
        read(UnitTmp_,*) StringTecHeader
     end do
     ! Read the coordinates and variables at the nodes
     do iNodeTec = 1, nNodeTec
        read(UnitTmp_,*) XyzTec_DI(:,iNodeTec), PlotVarTec_VI(:,iNodeTec)
     end do
     ! Read the node indexes of the nodes surrounding each cell center
     do iCellTec = 1, nCellTec
        read(UnitTmp_,*) iNodeTec_II(:,iCellTec)
     end do

     close(UnitTmp_)
  end if

  ! Collect info from all files and put it into PlotVar_VC and Coord_DC
  VolumeCheck = 0.0
  iCell = 0
  nCellCheck = 0
  l=len_trim(NameFileHead)
  do iProc = 0, nProc-1
     if(.not.DoReadTecplot)then
        if(    nProc > 100000)then
           write(NameFile,'(a,i6.6,a)')NameFileHead(1:l-2)//"_pe",iProc,'.idl'
        elseif(nProc > 10000)then
           write(NameFile,'(a,i5.5,a)')NameFileHead(1:l-2)//"_pe",iProc,'.idl'
        else
           write(NameFile,'(a,i4.4,a)')NameFileHead(1:l-2)//"_pe",iProc,'.idl'
        end if

        if(iProc==0)write(*,*)'reading files=',trim(NameFile),&
             '...',nProc-1,'.idl'

        if(DoReadBinary)then
           open(UnitTmp_, file=NameFile, status='old', form='unformatted', &
                iostat=iError)
        else
           open(UnitTmp_, file=NameFile, status='old', iostat=iError)
        end if

        ! Assume that missing files were empty. 
        if(iError /=0) CYCLE
     end if

     ! Read file
     do
        if(DoReadTecplot)then
           ! Collect information from the Tecplot arrays already read above
           iCellTec = nCellCheck + 1
           if(iCellTec > nCellTec) EXIT
           xMinTec   =  1e30
           xMaxTec   = -1e30
           Xyz_D     = 0.
           PlotVar_V = 0.
           ! Loop over nodes surrounding the cell center,
           ! average out the coordinates and the plot variables
           ! and calculate the cell size (only works for Cartesian)
           do iNeiTec = 1, 8
              iNodeTec  = iNodeTec_II(iNeiTec,iCellTec)
              xMinTec   = min(xMinTec, XyzTec_DI(1,iNodeTec))
              xMaxTec   = max(xMaxTec, XyzTec_DI(1,iNodeTec))
              Xyz_D     = Xyz_D + XyzTec_DI(:,iNodeTec)
              PlotVar_V = PlotVar_V + PlotVarTec_VI(:,iNodeTec)
           end do
           CellSize1 = xMaxTec - xMinTec
           Xyz_D     = Xyz_D/8.
           PlotVar_V = PlotVar_V/8.
        else if(DoReadBinary)then
           if(nByteRealRead == 4)then
              read(UnitTmp_,ERR=999,END=999) DxCell4, Xyz4_D, PlotVar4_V
              CellSize1 = DxCell4; Xyz_D = Xyz4_D; PlotVar_V = PlotVar4_V
           else
              read(UnitTmp_,ERR=999,END=999) DxCell8, Xyz8_D, PlotVar8_V
              CellSize1 = DxCell8; Xyz_D = Xyz8_D; PlotVar_V = PlotVar8_V
           end if
        else
           read(UnitTmp_,*,ERR=999,END=999) CellSize1, Xyz_D, PlotVar_V
        end if

        nCellCheck = nCellCheck + 1

        ! Set cell size based on 1st component and similarity with 
        ! the shape of the smallest cell
        CellSize_D = CellSize1*dCoordMin_D/dCoordMin_D(1)

        ! Set GenCoord_D, possibly modify Xyz_D too
        call set_gen_coord

        if(.not.IsStructured)then

           ! Simply put data into array. 
           ! Sorting and averaging will be done at the end.

           iCell = iCell + 1
           PlotVar_VC(:,iCell,1,1) = PlotVar_V
           do iDim = 1, nDim
              Coord_DC(iDim,iCell,1,1) = Xyz_D(iDimCut_D(iDim))
              GenCoord_DI(iDim,iCell)  = GenCoord_D(iDimCut_D(iDim))
           end do

           ! We are finished with unstructured
           CYCLE

        endif

        if(CellSize1 < dCoord1Plot + 1.e-6)then
           ! Cell has the correct size or finer
           Ijk_D = max(1, nint(( GenCoord_D - CoordMin_D)/dCoordPlot_D + 0.5))
           i = Ijk_D(1); j = Ijk_D(2); k = Ijk_D(3)

           if(CellSize1 < dCoord1Plot - 1.e-6)then
              ! Cell is finer, calculate VolumePlot fraction
              Fraction = (CellSize1/dCoord1Plot)**nDim
           else
              Fraction = 1.0
           end if

           PlotVar_VC(:,i,j,k) = PlotVar_VC(:,i,j,k) + Fraction*PlotVar_V
           do iDim = 1, nDim
              Coord_DC(iDim,i,j,k) = Coord_DC(iDim,i,j,k) &
                   + Fraction*Xyz_D(iDimCut_D(iDim))
           end do
           VolumeCheck = VolumeCheck + Fraction
        else
           ! Cell is coarser than required resolution
           IjkMin_D = min(nCell_D, max(1, &
                nint((GenCoord_D - 0.5*CellSize_D - CoordMin_D)/dCoordPlot_D &
                + 1)))
           IjkMax_D = min(nCell_D, max(1, &
                nint((GenCoord_D + 0.5*CellSize_D - CoordMin_D)/dCoordPlot_D)))
           iMin = IjkMin_D(1); jMin  = IjkMin_D(2); kMin = IjkMin_D(3)
           iMax = IjkMax_D(1); jMax  = IjkMax_D(2); kMax = IjkMax_D(3)

           ! First order prolongation
           do iVar = 1, nPlotVar
              PlotVar_VC(iVar,iMin:iMax,jMin:jMax,kMin:kMax)= &
                   PlotVar_VC(iVar,iMin:iMax,jMin:jMax,kMin:kMax) &
                   + PlotVar_V(iVar)
           end do
           do iDim = 1, nDim
              Coord_DC(iDim,iMin:iMax,jMin:jMax,kMin:kMax) = &
                   Coord_DC(iDim,iMin:iMax,jMin:jMax,kMin:kMax) &
                   + Xyz_D(iDimCut_D(iDim))
           end do

           if(iMax < iMin .or. jMax < jMin .or. kMax < kMin)&
                write(*,*)'!!! Empty box for cell CellSize_D(1), GenCoord_D=',&
                CellSize_D(1), GenCoord_D

           VolumeCheck = VolumeCheck &
                + (iMax - iMin + 1)*(jMax - jMin + 1)*(kMax - kMin + 1)
        end if
     end do ! read file

999  continue

     if(.not.DoReadTecplot) close(UnitTmp_)
  end do ! iProc

  if(IsVerbose)write(*,*)'nCellCheck=', nCellCheck, &                  
       ' nCellPlot=', nCellPlot

  if(nCellCheck /= nCellPlot)&
       write(*,*)'!!! Discrepancy: nCellCheck=', nCellCheck, &
       ' nCellPlot=', nCellPlot,' !!!'

  if(IsStructured)then
     VolumePlot = product(real(nCell_D))
     if(ndim==1 .and. abs(VolumeCheck/VolumePlot - 4.0) < 0.0001)then
        PlotVar_VC = 0.25*PlotVar_VC
        Coord_DC = 0.25*Coord_DC
        write(*,*)'Averaged 1D structured file everywhere'
     elseif(abs(VolumeCheck/VolumePlot - 2.0) < 0.0001)then
        PlotVar_VC = 0.5*PlotVar_VC
        Coord_DC = 0.5*Coord_DC
        write(*,*)'Averaged structured file everywhere'
     elseif(abs(VolumeCheck/VolumePlot - 1.0) > 0.0001)then
        write(*,*)'!!! Discrepancy in structured file:',&
             'filled VolumeCheck=',VolumeCheck,' VolumePlot=',VolumePlot,' !!!'
     end if
  else
     if(iCell /= nCellPlot) &
          write(*,*)'!!! Error: nCellPlot=',nCellPlot,' /= iCell=',iCell,' !!!'

     n1 = iCell
     nCell_D(1) = iCell
  end if

  if(NameFileHead(1:3) == 'cut')then

     ! Convert radians to degrees
     do iDim = 1, nDim
        select case(NameCoordPlot_D(iDim))
        case('r')
           ! Convert generalized coordinates to radius and degrees
           if(IsLogRadius)then
              Coord_DC(iDim,1:n1,:,:) = exp(Coord_DC(iDim,1:n1,:,:))
           elseif(IsGenRadius)then
              do k = 1, n3; do j = 1, n2; do i = 1, n1
                 Coord_DC(iDim,i,j,k) = exp(linear(LogRgen_I, 0, nRgen-1, &
                      Coord_DC(iDim,i,j,k)*(nRgen-1), DoExtrapolate=.true.) )
              end do; end do; end do
           end if

        case('phi', 'theta', 'lat')
           Coord_DC(iDim,1:n1,:,:) = Coord_DC(iDim,1:n1,:,:)*cRadToDeg
        end select
     end do

  endif

  if(TypeFile == 'tec')then
     NameFile = NameFileHead(1:l-2)//'.dat'
  else
     NameFile = NameFileHead(1:l-2)//'.out'
  end if
  write(*,*)'writing file =',trim(NameFile)

  ! Param_I is the combination of eqpar and ParamExtra_I
  allocate(Param_I(nParamPlot+nParamExtra))
  do i = 1, nParamPlot
     Param_I(i) = PlotParam_I(i)
  end do
  do i = 1, nParamExtra
     Param_I(i + nParamPlot) = ParamExtra_I(i)
  end do

  if(IsVerbose)then
     write(*,*)'nParamPlot, nParamExtra=', nParamPlot, nParamExtra
     write(*,*)' Param_I=',  Param_I
  end if

  if(.not.IsStructured .and. (nDim < 3 .or. DoSort3D))then
     if(IsVerbose)write(*,*)'Sorting unstructured grid points'

     ! Sort points based on (generalized) coordinates
     allocate(Sort_I(n1), iSort_I(n1), STAT=iError)
     if(iError /= 0) stop 'PostIDL.exe ERROR: could not allocate sort arrays'

     ! Form sorting function from the generalized coordinates
     Sort_I = GenCoord_DI(1,:)
     if(nDim > 1) Sort_I = Sort_I + exp(1.0)*GenCoord_DI(2,:)
     if(nDim > 2) Sort_I = Sort_I + exp(2.0)*GenCoord_DI(3,:)

     ! Sort points according to the sorting function
     call sort_quick(n1, Sort_I, iSort_I)

     Coord_DC(:,:,1,1) = Coord_DC(:,iSort_I,1,1)
     PlotVar_VC(:,:,1,1) = PlotVar_VC(:,iSort_I,1,1)

     if(IsVerbose)write(*,*)'Sorting is done'

     ! Average out coinciding points
     if(nDim < 3) then
        if(IsVerbose)write(*,*)'Averaging coinciding points'
        GenCoord_DI = GenCoord_DI(:,iSort_I)

        allocate(StateSum_V(nPlotVar), CellSizeMin_D(nDim))
        CellSizeMin_D = dCoordMin_D(iDimCut_D(1:nDim))
        i = 1
        k = 1
        do while(i < n1)
           StateSum_V = PlotVar_VC(:,i,1,1)
           nSum       = 1
           j = i + 1
           do while( sum(abs(GenCoord_DI(:,j) - GenCoord_DI(:,i)) &
                /CellSizeMin_D) < 0.01)
              StateSum_V = StateSum_V + PlotVar_VC(:,j,1,1)
              nSum = nSum + 1
              j = j + 1
              if(j > n1) EXIT
           end do
           if(j > i+1) then
              ! Put average value into i-th element
              PlotVar_VC(:,i,1,1) = StateSum_V/nSum              
           end if
           ! Save the index for the unique coordinates 
           iSort_I(k) = i
           k = k + 1
           i = j 
        end do
        deallocate(StateSum_V, CellSizeMin_D, GenCoord_DI)

        ! Special judgement for the last point
        if(j == n1) then
           iSort_I(k) = j
           n1 = k
        else
           n1 = k - 1
        end if
        
        ! move the elements after finding out all the coinciding ones 
        Coord_DC(:,1:n1,1,1) = Coord_DC(:,iSort_I(1:n1),1,1)
        PlotVar_VC(:,1:n1,1,1) = PlotVar_VC(:,iSort_I(1:n1),1,1)
     
        if(IsVerbose)write(*,*)'Averaging done'
     end if
  
     deallocate(Sort_I, iSort_I)

     if(IsVerbose)then
        write(*,*)'After sorting and averaging n1=', n1
     end if
  end if

  if(IsVerbose)then
     write(*,*)'shape(Coord_DC)  =', shape(Coord_DC)
     write(*,*)'shape(PlotVar_VC)=', shape(PlotVar_VC)
  end if

  ! the sizes of Coord_DC and PlotVar_VC may be modified by cell averaging 
  ! in unstructured grids. Only the first dimension (1:n1) needs to be set
  call save_plot_file(NameFile,&
       TypeFileIn = TypeFile, &
       StringHeaderIn = NameUnit, &
       nStepIn = nStep, TimeIn = t, &
       ParamIn_I = Param_I, &
       NameVarIn = NameVar, &
       IsCartesianIn = TypeGeometry=='cartesian' .and. IsStructured,&
       nDimIn = nDim,&
       CoordIn_DIII = Coord_DC(:,1:n1,:,:), & 
       VarIn_VIII = PlotVar_VC(:,1:n1,:,:))

  deallocate(Coord_DC, PlotVar_VC, Param_I, PlotVar_V)
  if(allocated(PlotParam_I))  deallocate(PlotParam_I)
  if(allocated(PlotVar4_V))   deallocate(PlotVar4_V)
  if(allocated(PlotVar8_V))   deallocate(PlotVar8_V)

  write(*,'(a)')'PostIDL finished'

contains
  !===========================================================================
  subroutine set_gen_coord

    ! Calculate the generalized coordinates GenCoord_D from Xyz_D
    ! This routine is similar but not exactly the same as 
    ! BATL_geometry:xyz_to_coord

    real:: rCyl ! distance from axis 

    ! Toroidal variables

    ! Rotation matrix for rotated Cartesian grid
    real, allocatable, save:: GridRot_DD(:,:)

    ! List of vector variables for rotated Cartesian coordinates
    integer:: nVector = 0
    integer, allocatable, save:: iVarVector_I(:)

    ! Temporary variables to process the variable name string
    integer:: iVector, iVar, l
    character(len=10), allocatable:: NameVar_V(:)
    character(len=10)             :: NameVector1, NameVector

    ! Variables for roundcube coordinate transformation
    real :: r2, Dist1, Dist2, Coef1, Coef2
    !---------------------------------------------------------------------
    if(TypeGeometry == 'cartesian' .or. NameFileHead(1:3) == 'cut')then
       GenCoord_D = Xyz_D
       
       RETURN
    end if

    if(TypeGeometry == 'roundcube')then

       r2 = sum(Xyz_D**2)
       if (r2 > 0.0) then
          ! L1 and L2 distance
          Dist1 = maxval(abs(Xyz_D))
          Dist2 = sqrt(r2)
          if (rRound1 > rRound0 ) then
             ! The rounded (distorted) grid is outside the non-distorted part
             if (Dist1 > rRound0) then
                ! Outside the undistorted region
                ! We have to solve a quadratic equation of w. 
                ! w^2 - Coef1*w- Coef2 = 0
                Coef1 = -1 + &
                     rRound0/(rRound1-rRound0)*(dist1*SqrtNDim/Dist2 - 1)
                Coef2 = Dist1/(rRound1-rRound0)*(SqrtNDim*Dist1/Dist2 - 1)
                GenCoord_D = Xyz_D/(-Coef1 + sqrt(Coef1**2 + 4*Coef2))*2
             else
                ! No distortion
                GenCoord_D = Xyz_D
             end if
         
          else
             ! The rounded (distorted) grid is inside of the non-distorted part
             if (Dist2 < rRound1) then
                ! Solving w^2 - w + Coef1 = 0
                Coef1 = Dist1/rRound1*(1 - Dist1/Dist2)
                GenCoord_D = Xyz_D / (1 + sqrt(1-4*Coef1))*2
             else
                ! Solving w^2 + Coef1*w + Coef2 = 0
                Coef1 = -1 + (1 - Dist1/Dist2)/(rRound0-rRound1)*rRound0
                Coef2 = -(1 - Dist1/Dist2)/(rRound0 - rRound1)*Dist1
                Coef2 = (-Coef1 + sqrt(Coef1**2 - 4*Coef2))*0.5
                GenCoord_D = Xyz_D / Coef2
             end if
          end if
           
       else
          GenCoord_D = 0.0
       end if
       RETURN
    end if

    if(TypeGeometry(1:7)=='rotated')then

       if(.not.allocated(GridRot_DD))then
          ! Setup rotation matrix
          allocate(GridRot_DD(3,3))

          ! The rotation matrix should be the same as in BATL_geometry
          GridRot_DD = rot_matrix_z(0.6,0.8)

          ! Find vectors
          allocate(NameVar_V(nDim + nPlotVar + nParamPlot))
          call split_string(NameVar, NameVar_V)

          ! Make this array large enough
          allocate(iVarVector_I(nPlotVar/3))

          NameVector = ' '
          do iVar = 1, nPlotVar - 2
             ! Add nDim to skip the coordinate names
             NameVector1 = NameVar_V(nDim+iVar)
             call lower_case(NameVector1)

             l = len_trim(NameVector1)

             ! Identify vectors as 3 strings ending with x, y, z
             if(NameVector1(l:l) /= 'x') CYCLE

             ! Prospective vector component
             NameVector = NameVector1(1:l-1)

             ! Check the next two names
             NameVector1 = NameVar_V(nDim+iVar+1)
             call lower_case(NameVector1)
             if(NameVector1 /= trim(NameVector)//'y') CYCLE

             NameVector1 = NameVar_V(nDim+iVar+2)
             call lower_case(NameVector1)
             if(NameVector1 /= trim(NameVector)//'z') CYCLE

             nVector = nVector + 1
             iVarVector_I(nVector) = iVar
          end do
          deallocate(NameVar_V)

          !write(*,*)'nVector, iVarVector_I=', nVector, iVarVector_I(1:nVector)

       end if

       ! Unrotate the coordinates for comparison with Cartesian runs
       Xyz_D = matmul(Xyz_D, GridRot_DD)
       GenCoord_D = Xyz_D

       ! Unrotate vectors
       do iVector = 1, nVector
          iVar = iVarVector_I(iVector)
          PlotVar_V(iVar:iVar+2) = matmul(PlotVar_V(iVar:iVar+2), GridRot_DD)
       end do

       RETURN
    end if

    rCyl = sqrt(Xyz_D(1)**2 + Xyz_D(2)**2)

    ! Calculate phi
    if (rCyl == 0.0) then
       GenCoord_D(2) = 0.0
    else
       GenCoord_D(2) = &
            modulo(atan2(Xyz_D(2), Xyz_D(1)) - CoordMin_D(2), cTwoPi) &
            + CoordMin_D(2)
    end if

    select case(TypeGeometry)
    case('cylindrical', 'cylindrical_lnr', 'cylindrical_genr')
       GenCoord_D(1) = rCyl
       GenCoord_D(3) = Xyz_D(3)

       if(nDim==2)then
          ! Set the 'X-Y' coordinates for plotting a 2D cut
          select case(iDim0)
          case(1)
             ! This is R=const slice, use longitude [deg] vs height
             Xyz_D(2)=GenCoord_D(2)*cRadToDeg
          case(2)
             ! This is x=0 or y=0 plane, use signed radius vs Z
             Xyz_D(1) = sign(1.0, Xyz_D(1)+Xyz_D(2))*rCyl
             ! Radial distance
             GenCoord_D(1) = rCyl
             ! The generalized coordinate runs from rMin to 2*rMax-rMin
             if(UseDoubleCut .and. Xyz_D(1) < 0.0) &
                  GenCoord_D(1) = GenCoord_D(1) + CoordShift1
          end select
       end if
    case('spherical', 'spherical_lnr', 'spherical_genr')
       GenCoord_D(1) = sqrt(rCyl**2 + Xyz_D(3)**2)

       if(nDim==2)then
          ! Set the 'X-Y' coordinates for plotting a 2D cut
          select case(iDim0)
          case(1)
             ! This is R=const slice, use longitude vs latitude in degs.
             Xyz_D(2:3)=GenCoord_D(2:3)*cRadToDeg
          case(2)
             ! This is x=0 or y=0 plane, use axial radius vs Z
             Xyz_D(1) = sign(1.0, Xyz_D(1) + Xyz_D(2))*rCyl
          case(3)
             ! This is the z=0 plane
             ! Stretch X and Y with rSph/rCyl instead of simply
             ! projecting the points down to the X-Y plane
             Xyz_D(1:2) = Xyz_D(1:2)*GenCoord_D(1)/rCyl
          end select
       end if
       ! Latitude
       GenCoord_D(3) = asin(Xyz_D(3)/GenCoord_D(1))
       ! Shift by width of latitude range for the left half 
       if(UseDoubleCut .and. Xyz_D(1) < 0.0) & 
            GenCoord_D(3) = GenCoord_D(3) + CoordShift2

    case default
       write(*,*)'Unknown TypeGeometry='//TypeGeometry
       stop
    end select

    if(IsLogRadius .or. IsGenRadius) GenCoord_D(1) = log(GenCoord_D(1))

    if(IsGenRadius)then
       i = min(nRgen-1, count(LogRgen_I < GenCoord_D(1)))
       GenCoord_D(1) = ( i -1  &
            + (GenCoord_D(1) - LogRgen_I(i)) &
            / (LogRgen_I(i+1) - LogRgen_I(i)) )&
            / (nRgen - 1)
    end if

  end subroutine set_gen_coord

end program post_idl

!=============================================================================

subroutine CON_stop(String)

  ! This routine is needed for ModPlotFile

  implicit none

  character(len=*), intent(in):: String
  write(*,*) 'ERROR in PostIDL: '//String
  stop

end subroutine CON_stop

!=============================================================================
