!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPlotFile

  ! Save or read VAC/IDL type plotfiles from 1 up to 3 dimensions.
  ! ASCII, single (real4), or double (real8) precision binary file formats 
  ! can be used.
  ! The plot file contains 5 header lines: 
  !
  !    Header
  !    nStep Time nDim nParam nVar
  !    n1 .. nNDim
  !    Param1 .. ParamNParam
  !    NameVar
  !
  ! Header   (string) describes the plot. It is up to 500 characters.
  ! nStep    (integer) is the number of time steps/iterations etc.
  ! Time     (real) is the simulation time.
  ! nDim     (integer) number of dimensions. Negative for non-Cartesian grids.
  ! nParam   (integer) number of parameters (for example adiabatic index)
  ! n1 ..    (nDim integers) grid sizes in the nDim dimensions
  ! Param1.. (nParam reals) parameter values
  ! NameVar  (string) space separated list of names for the 
  !                   nDim coordinates, nVar variables and nParam parameters
  !
  ! The header is followed by the coordinate/variable data. 
  ! In ASCII files each line contains the coordinates+variables for one cell.
  ! The cells are ordered such that the first coordinate index changes fastest.
  ! In binary files the coordinates are saved in a single array, followed by
  ! the variables, saved as one record per variable.
  !
  ! For save_plot_file the coordinates can be either given as full arrays,
  ! or for Cartesian grids they can be given as 1D arrays for each dimension,
  ! or for uniform Cartesian grids they can be given with min and max values.
  ! The number of dimensions, the size of the grid and the number of the 
  ! variables is determined from the size of the variable array.
  !
  ! For read_plot_file the number of dimensions, variables, parameters, grid
  ! size are optionaal parameters.

  use ModIoUnit,    ONLY: UnitTmp_
  use ModKind,      ONLY: Real4_
  use ModHdf5Utils, ONLY: save_hdf5_file
  implicit none

  private ! except

  public:: read_plot_file
  public:: save_plot_file
  public:: test_plot_file

  integer, parameter :: MaxDim = 3

contains

  !=========================================================================

  subroutine save_plot_file(NameFile, TypePositionIn, &
       TypeFileIn, StringHeaderIn, nStepIn, TimeIn, &
       ParamIn_I, NameVarIn, NameVarIn_I, NameUnitsIn,&
       IsCartesianIn, &
       nDimIn,&
       CoordMinIn_D, CoordMaxIn_D, &
       Coord1In_I, Coord2In_I, Coord3In_I, &
       CoordIn_I, CoordIn_DII, CoordIn_DIII,&
       VarIn_I,  VarIn_II,  VarIn_III,  &
       VarIn_VI, VarIn_VII, VarIn_VIII, &
       VarIn_IV, VarIn_IIV, VarIn_IIIV, iCommIn)

    use ModUtilities, ONLY: split_string, join_string

    character(len=*),           intent(in):: NameFile       ! Name of plot file
    character(len=*), optional, intent(in):: TypePositionIn !asis/rewind/append
    character(len=*), optional, intent(in):: TypeFileIn     ! ascii/real8/real4
    character(len=*), optional, intent(in):: StringHeaderIn ! header line
    integer,          optional, intent(in):: nStepIn        ! number of steps
    real,             optional, intent(in):: TimeIn         ! simulation time  
    real,             optional, intent(in):: ParamIn_I(:)   ! parameters
    character(len=*), optional, intent(in):: NameVarIn      ! list of names
    character(len=*), optional, intent(in):: NameVarIn_I(:) ! list of names 
    character(len=*), optional, intent(in):: NameUnitsIn    ! list of units
    logical,          optional, intent(in):: IsCartesianIn  ! Cartesian grid?
    integer,          optional, intent(in):: nDimIn         ! grid dimensions
    real,             optional, intent(in):: CoordIn_I(:)   ! coords in 1D
    real,             optional, intent(in):: CoordIn_DII(:,:,:)       ! 2D
    real,             optional, intent(in):: CoordIn_DIII(:,:,:,:)    ! 3D
    real,             optional, intent(in):: Coord1In_I(:)  ! coords for axis 1
    real,             optional, intent(in):: Coord2In_I(:)  ! coords for axis 2
    real,             optional, intent(in):: Coord3In_I(:)  ! coords for axis 3
    real,             optional, intent(in):: CoordMinIn_D(:)! min coordinates
    real,             optional, intent(in):: CoordMaxIn_D(:)! max coordinates
    real,             optional, intent(in):: VarIn_I(:)     ! variable  in 1D
    real,             optional, intent(in):: VarIn_II(:,:)               ! 2D
    real,             optional, intent(in):: VarIn_III(:,:,:)            ! 3D
    real,             optional, intent(in):: VarIn_VI(:,:)  ! variables in 1D
    real,             optional, intent(in):: VarIn_VII(:,:,:)            ! 2D
    real,             optional, intent(in):: VarIn_VIII(:,:,:,:)         ! 3D
    real,             optional, intent(in):: VarIn_IV(:,:)  ! variables in 1D
    real,             optional, intent(in):: VarIn_IIV(:,:,:)            ! 2D
    real,             optional, intent(in):: VarIn_IIIV(:,:,:,:)         ! 3D
    integer,          optional, intent(in):: iCommIn ! MPI communicator for HDF

    character(len=10)  :: TypePosition
    character(len=10)  :: TypeStatus
    character(len=20), allocatable  :: NameVar_I(:)
    character(len=20)  :: TypeFile
    character(len=500) :: StringHeader
    character(len=500) :: NameVar,NameUnits
    integer :: nStep, nDim, nParam, nVar, n1, n2, n3
    integer :: nCellsPerBlock(3), iBlk, nBlocks
    integer :: nBlocksXYZ(3), iG,jG,kG
    integer, allocatable:: MinimumBlockIjk(:,:)
    real             :: Time, Coord, ResultMod, iBlockDxDyDz(3)
    logical          :: IsCartesian
    real, allocatable:: Param_I(:), Coord_ID(:,:), Var_IV(:,:)
    real, allocatable:: VarHdf5Output(:,:,:,:,:), XYZMinMax(:,:,:)
    real(Real4_), allocatable:: Param4_I(:), Coord4_ID(:,:), Var4_I(:)

    integer :: n_D(0:MaxDim),ii,jj,kk
    integer :: i, j, k, i_D(3), iDim, iVar, n, nDimOut, iError
    logical  :: IsSplitSuccessfull

    character(len=*), parameter:: NameSub = 'save_plot_file'
    !---------------------------------------------------------------------
    ! either write a new file (remove old one if any)
    ! or append to an existing file
    TypePosition = 'rewind'
    if(present(TypePositionIn))TypePosition = TypePositionIn
    TypeStatus = 'replace'
    if(TypePosition == 'append')TypeStatus = 'unknown'

    TypeFile = 'ascii'
    if(present(TypeFileIn)) TypeFile = TypeFileIn
    StringHeader = 'No header info'
    if(present(StringHeaderIn)) StringHeader = StringHeaderIn
    nStep = 0
    if(present(nStepIn)) nStep = nStepIn
    Time = 0.0
    if(present(TimeIn)) Time = TimeIn

    if(present(ParamIn_I))then
       nParam = size(ParamIn_I)
       allocate(Param_I(nParam))
       Param_I = ParamIn_I
    else
       nParam = 0
    end if
    ! Figure out grid dimensions and number of variables. Default is 1.
    n_D = 1
    if(present(VarIn_I))then
       nDim = 1
       n_D(1:1) = shape(VarIn_I)
    elseif(present(VarIn_II)) then
       nDim = 2
       n_D(1:2) = shape(VarIn_II)
    elseif(present(VarIn_III)) then
       nDim = 3
       n_D(1:3) = shape(VarIn_III)
    elseif(present(VarIn_VI))then
       nDim = 1
       n_D(0:1) = shape(VarIn_VI)
    elseif(present(VarIn_VII))then
       nDim = 2
       n_D(0:2) = shape(VarIn_VII)
    elseif(present(VarIn_VIII))then
       nDim = 3
       n_D(0:3) = shape(VarIn_VIII) 
       ! For IV, IIV, IIIV types 
    elseif(present(VarIn_IV))then
       nDim = 1
       n_D(0:1) = shape(VarIn_IV)
       n_D(0:1) = cshift(n_D(0:1), -1)   ! shift nVar/n_D(1) to n_D(0)  
    elseif(present(VarIn_IIV))then
       nDim = 2
       n_D(0:2) = shape(VarIn_IIV)
       n_D(0:2) = cshift(n_D(0:2), -1)   ! shift nVar/n_D(2) to n_D(0)
    elseif(present(VarIn_IIIV))then
       nDim = 3
       n_D(0:3) = shape(VarIn_IIIV)
       n_D = cshift(n_D, -1)        ! shift nVar/n_D(3) to n_D(0)
    else
       call CON_stop(NameSub // &
            ': none of the VarIn_* variables are present')
    endif
    ! Extract information
    nVar = n_D(0)
    n1   = n_D(1)
    n2   = n_D(2)
    n3   = n_D(3)

    ! The plot dimension may be different from the dimensionality of VarIn
    if(present(nDimIn))then
       nDim = nDimIn
       if(n1 == 1 .and. n2 == 1)then
          n_D(1:3) = (/ n3, 1, 1/)
       elseif(n1 == 1)then
          n_D(1:3) = (/ n2, n3, 1/)
       elseif(n2 == 1)then
          n_D(1:3) = (/ n1, n3, 1/)
       end if
    end if
    IsCartesian = .true.
    if(present(IsCartesianIn)) IsCartesian = IsCartesianIn

    ! nDim is saved with a negative sign for non-Cartesian grid
    nDimOut = nDim
    if(.not. IsCartesian) nDimOut = -nDim

    ! Set variable names
    if(present(NameVarIn))then
       NameVar = NameVarIn
    else if(present(NameVarIn_I)) then
       call join_string(NameVarIn_I, NameVar)
    else
       ! Create some arbitrary variable names
       NameVar = 'x1'
       do i = 2, nDim
          write(NameVar, "(a, i1)") trim(NameVar) // ' x', i
       end do
       do i = 1, nVar
          write(NameVar, "(a, i2.2)") trim(NameVar) // ' v', i
       end do
       do i = 1, nParam
          write(NameVar, "(a, i2.2)") trim(NameVar) // ' p', i
       end do
    end if

    ! Create a variable name array
    allocate(NameVar_I(nDim + nVar + nParam))
    call split_string(NameVar, NameVar_I, i, UseArraySyntaxIn=.true.)
    if(i /= nDim + nVar + nParam)then
       write(*,*) NameSub,': NameFile=', trim(NameFile)
       write(*,*) NameSub,': NameVar=', trim(NameVar)
       write(*,*) NameSub,': number of substrings=', i
       write(*,*) NameSub,': nDim, nVar, nParam, sum=',&
            nDim, nVar, nParam, nDim + nVar + nParam
       call CON_stop(NameSub// &
            ': number of names in NameVar does not match nDim+nVar+nParam !')
    end if

    ! Allocate arrays with a shape that is convenient for saving data
    if(TypeFile == 'hdf5') then
       ! VisIt is much much faster if you give it blocks so it can parallelize
       ! and on some machines hdf5 is faster in parallel.
       ! this routine can be called in serial or parallel for hdf5. 
       ! Just calling it with iproc = 0 seems to be the best thing for now.

       nCellsPerBlock = 1
       nBlocksXYZ = 1
       do i = 1, nDim
          if (n_D(i) > 1) then
             do j= 4, 25
                ResultMod = mod(n_D(i), j)
                if (ResultMod == 0) then 
                   IsSplitSuccessfull = .true.
                else
                   IsSplitSuccessfull = .false.
                end if
                if(IsSplitSuccessfull) exit
             end do
             if(.not. IsSplitSuccessfull ) then
                do j=3,2,-1
                   ResultMod = mod(n_D(i), j)
                   if (ResultMod == 0) then 
                      IsSplitSuccessfull = .true.
                   else
                      IsSplitSuccessfull = .false.
                   end if
                   if(IsSplitSuccessfull) exit
                end do
             end if

             if(IsSplitSuccessfull ) then
                nCellsPerBlock(i) = j
             else 
                nCellsPerBlock(i) = n_D(i)
             end if
          else
             nCellsPerBlock(i) = 1
          end if
       end do

       do i = 1, nDim
          nBlocksXYZ(i)=n_D(i)/nCellsPerBlock(i)
          iBlockDxDyDz(i) = (CoordMaxIn_D(i) - CoordMinIn_D(i))/nBlocksXYZ(i)
       end do
       nBlocks = product(nBlocksXYZ(1:nDim))

       allocate(VarHdf5Output(nCellsPerBlock(1),nCellsPerBlock(2),&
            nCellsPerBlock(3),nBlocks, nVar))
       allocate(XYZMinMax(2,nDim,nBlocks))
       allocate(MinimumBlockIjk(nDim, nBlocks))
       iBlk = 0
       ! do kk=1, nBlocksXYZ(3); do jj=1, nBlocksXYZ(2); do ii=1,nBlocksXYZ(1);
       do kk=1,1; do jj=1, nBlocksXYZ(2); do ii=1,nBlocksXYZ(1);
          iBlk = iBlk +1
          do k = 1, 1; do j = 1, nCellsPerBlock(2); do i = 1,nCellsPerBlock(1)
             ! do k = 1, nCellsPerBlock(3); do j = 1, nCellsPerBlock(2); do i = 1,nCellsPerBlock(1)
             iG = (ii - 1)*nCellsPerBlock(1) + i
             jG = (jj - 1)*nCellsPerBlock(2) + j
             kG = (kk - 1)*nCellsPerBlock(3) + k
             if(present(VarIn_I)) then
                VarHdf5Output(i,1,1,iBlk,1)      = VarIn_I(iG)
             elseif(present(VarIn_II)) then
                VarHdf5Output(i,j,1,iBlk,1)      = VarIn_II(iG,jG)
             elseif(present(VarIn_III)) then
                VarHdf5Output(i,j,k,iBlk,1)      = VarIn_III(iG,jG,kG)
             elseif(present(VarIn_VI)) then
                VarHdf5Output(i,1,1,iBlk,1:nVar) = VarIn_VI(1:nVar,iG)
             elseif(present(VarIn_VII)) then
                VarHdf5Output(i,j,1,iBlk,1:nVar) = VarIn_VII(1:nVar,iG,jG)
             elseif(present(VarIn_VIII)) then
                VarHdf5Output(i,j,k,iBlk,1:nVar)= VarIn_VIII(1:nVar,iG,jG,kG)
             elseif(present(VarIn_IV)) then
                VarHdf5Output(i,1,1,iBlk,1:nVar) = VarIn_IV(iG,1:nVar)
             elseif(present(VarIn_IIV)) then
                VarHdf5Output(i,j,1,iBlk,1:nVar) = VarIn_IIV(iG,jG,1:nVar)
             elseif(present(VarIn_IIIV)) then
                VarHdf5Output(i,j,k,iBlk,1:nVar)= VarIn_IIIV(iG,jG,kG,1:nVar)
             endif
          end do; end do; end do;
          do n = 1, nDim
             if(n==1) then
                MinimumBlockIjk(n, iBlk) = (ii-1)*nCellsPerBlock(n) 
                XYZMinMax(1,n,iBlk) = iBlockDxDyDz(n)*(ii-1) + CoordMinIn_D(n)
                XYZMinMax(2,n,iBlk) = iBlockDxDyDz(n)*ii+ CoordMinIn_D(n)
             else if(n==2) then
                MinimumBlockIjk(2, iBlk) = (jj-1)*nCellsPerBlock(n) 
                XYZMinMax(1,n,iBlk) = iBlockDxDyDz(n)*(jj-1)+ CoordMinIn_D(n)
                XYZMinMax(2,n,iBlk) = iBlockDxDyDz(n)*jj+ CoordMinIn_D(n)
             else if(n==3) then 
                MinimumBlockIjk(n, iBlk) = (kk-1)*nCellsPerBlock(n)
                XYZMinMax(1,n,iBlk) = iBlockDxDyDz(n)*(kk-1)+ CoordMinIn_D(n)
                XYZMinMax(2,n,iBlk) = iBlockDxDyDz(n)*kk+ CoordMinIn_D(n)
             end if
          end do

       end do; end do; end do;
    else
       allocate(Coord_ID(n1*n2*n3,nDim), Var_IV(n1*n2*n3,nVar))
       ! Fill in the Coord_ID coordinate array using the available information
       do iDim = 1, nDim
          n = 0
          do k = 1, n3; do j = 1, n2; do i = 1, n1
             n = n + 1
             Coord = huge(1.0)
             if(present(CoordIn_I))    Coord = CoordIn_I(i)
             if(present(CoordIn_DII))  Coord = CoordIn_DII(iDim,i,j)
             if(present(CoordIn_DIII)) Coord = CoordIn_DIII(iDim,i,j,k)
             if(present(Coord1In_I) .and. iDim==1) Coord = Coord1In_I(i)
             if(present(Coord2In_I) .and. iDim==2) Coord = Coord2In_I(j)
             if(present(Coord3In_I) .and. iDim==3) Coord = Coord3In_I(k)
             if(present(CoordMinIn_D)) then
                i_D = (/i, j, k/)
                Coord = CoordMinIn_D(iDim) + (i_D(iDim)-1)* &
                     ((CoordMaxIn_D(iDim) - CoordMinIn_D(iDim))/(n_D(iDim)-1))
             end if
             Coord_ID(n, iDim) = Coord
          end do; end do; end do; 
       end do

       ! Check if all coordinates were set
       if(any(Coord_ID == huge(1.0))) call CON_stop(NameSub // & 
            ' coordinates were not defined')

       ! Fill in the Var_IV variable array using the available information
       Var_IV = huge(1.0)
       do iVar = 1, nVar
          n = 0
          do k = 1, n3; do j = 1, n2; do i = 1,n1
             n = n + 1
             if(present(VarIn_I))    Var_IV(n,iVar) = VarIn_I(i)
             if(present(VarIn_II))   Var_IV(n,iVar) = VarIn_II(i,j)
             if(present(VarIn_III))  Var_IV(n,iVar) = VarIn_III(i,j,k)
             if(present(VarIn_VI))   Var_IV(n,iVar) = VarIn_VI(iVar,i)
             if(present(VarIn_VII))  Var_IV(n,iVar) = VarIn_VII(iVar,i,j)
             if(present(VarIn_VIII)) Var_IV(n,iVar) = VarIn_VIII(iVar,i,j,k)
             if(present(VarIn_IV))   Var_IV(n,iVar) = VarIn_IV(i,iVar)
             if(present(VarIn_IIV))  Var_IV(n,iVar) = VarIn_IIV(i,j,iVar)
             if(present(VarIn_IIIV)) Var_IV(n,iVar) = VarIn_IIIV(i,j,k,iVar)
          end do; end do; end do; 
       end do

       ! Check if all variables were set
       if(any(Var_IV == huge(1.0))) call CON_stop(NameSub // & 
            ' variables were not defined')
    end if

    select case(TypeFile)
    case('hdf5')
       if (present(NameUnitsIn)) then
          NameUnits = NameUnitsIn
       else
          NameUnits = ''
          do iVar=1, nVar
             NameUnits = NameUnits//'normalized '
          end do
       end if
       call save_hdf5_file(NameFile,TypePosition, TypeStatus, StringHeader,&
            nStep, nBlocks, Time, nDim, nParam, nVar,&
            nCellsPerBlock(1:nDim), NameVar_I(nDim+1:nDim+nVar), NameUnits, &
            MinimumBlockIjk, XYZMinMax, VarHdf5Output, iComm=iCommIn, &
            CoordMin=CoordMinIn_D, CoordMax=CoordMaxIn_D)
       deallocate(VarHdf5Output)
       deallocate(XYZMinMax)
       deallocate(MinimumBlockIjk)
    case('tec')
       open(UnitTmp_, file=NameFile, &
            position=TypePosition, status=TypeStatus, iostat=iError)
       if(iError /= 0)call CON_stop(NameSub // &
            ' could not open tecplot file=' // trim(NameFile))
       write(UnitTmp_, "(a)", ADVANCE="NO") 'VARIABLES='
       if(nDim == 3) write(UnitTmp_, "(a)", ADVANCE="NO") '"K", '
       if(nDim >= 2) write(UnitTmp_, "(a)", ADVANCE="NO") '"J", '
       write(UnitTmp_, "(a)", ADVANCE="NO") '"I", '
       call join_string(NameVar_I(1:nDim+nVar), NameVar, '", "')
       write(UnitTmp_, "(a)") '"'//trim(NameVar)//'"'
       write(UnitTmp_,'(a,i6,a,i6,a,i6,a)') &
            'ZONE T="'//trim(StringHeader)// &
            '", I=',n1,', J=',n2,', K=',n3,', F=POINT'
       write(UnitTmp_,'(a,i8,a)') 'AUXDATA ITER="', nStep, '"'
       write(UnitTmp_,'(a,es18.10,a)') 'AUXDATA TIMESIM="', Time, '"'
       do i = 1, nParam
          write(UnitTmp_,'(a,100es18.10,a)') &
               'AUXDATA '//trim(NameVar_I(nDim+nVar+i))//'="', Param_I(i)
       end do
       select case(nDim)
       case(1)
          do i = 1, n1
             write(UNITTMP_,'(i8,100es18.10)') &
                  i, Coord_ID(n,:), Var_IV(n,:)
          end do
       case(2)
          do j = 1, n2; do i = 1, n1
             n = i + n1*(j-1)
             write(UNITTMP_,'(2i6,100es18.10)') &
                  j, i, Coord_ID(n,:), Var_IV(n, :)
          end do; end do
       case(3)
          do k = 1, n3; do j = 1, n2; do i = 1, n1
             n = i + n1*(j-1) + n1*n2*(k-1)
             write(UNITTMP_,'(3i6,100es18.10)') &
                  k, j, i, Coord_ID(n,:), Var_IV(n, :)
          end do; end do; end do
       end select
       close(UnitTmp_)
    case('formatted', 'ascii')
       open(UnitTmp_, file=NameFile, &
            position = TypePosition, status=TypeStatus, iostat=iError)
       if(iError /= 0)call CON_stop(NameSub // &
            ' could not open ascii file=' // trim(NameFile))

       write(UnitTmp_, "(a)")             trim(StringHeader)
       write(UnitTmp_, "(i7,es18.10,3i3)") nStep, Time, nDimOut, nParam, nVar
       write(UnitTmp_, "(3i8)")           n_D(1:nDim)
       if(nParam > 0) write(UnitTmp_, "(100es18.10)")     Param_I
       write(UnitTmp_, "(a)")             trim(NameVar)

       where(abs(Var_IV) < 1d-99) Var_IV = 0.0

       ! write out coordinates and variables line by line
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          write(UnitTmp_, "(100es18.10)") Coord_ID(n,:), Var_IV(n, :) 
       end do; end do; end do
       close(UnitTmp_)
    case('real8')
       open(UnitTmp_, file=NameFile, form='unformatted', &
            position=TypePosition, status=TypeStatus, iostat=iError)
       if(iError /= 0)call CON_stop(NameSub // &
            ' could not open real8 file=' // trim(NameFile))
       write(UnitTmp_) StringHeader
       write(UnitTmp_) nStep, Time, nDimOut, nParam, nVar
       write(UnitTmp_) n_D(1:nDim)
       if(nParam > 0)write(UnitTmp_) Param_I
       write(UnitTmp_) NameVar
       write(UnitTmp_) Coord_ID
       ! write out variables 1 by 1 to avoid segmentation fault 
       ! for very large Var_IV array
       do iVar = 1, nVar
          write(UnitTmp_) Var_IV(:,iVar)
       end do
       close(UnitTmp_)
    case('real4')
       open(UnitTmp_, file=NameFile, form='unformatted', &
            position = TypePosition, status=TypeStatus, iostat=iError)
       if(iError /= 0)call CON_stop(NameSub // &
            ' could not open real4 file=' // trim(NameFile))

       write(UnitTmp_) StringHeader
       write(UnitTmp_) nStep, real(Time, Real4_), nDimOut, nParam, nVar
       write(UnitTmp_) n_D(1:nDim)
       if(nParam > 0)then
          allocate(Param4_I(nParam))       
          Param4_I = Param_I
          write(UnitTmp_) Param4_I
          deallocate(Param4_I)
       end if
       write(UnitTmp_) NameVar
       ! Copy into single precision arrays to avoid compiler issues.
       allocate(Coord4_ID(n1*n2*n3, nDim))
       Coord4_ID = Coord_ID
       write(UnitTmp_) Coord4_ID
       deallocate(Coord4_ID)
       allocate(Var4_I(n1*n2*n3))
       do iVar = 1, nVar
          Var4_I = Var_IV(:,iVar)
          write(UnitTmp_) Var4_I
       end do
       deallocate(Var4_I)
       close(UnitTmp_)
    case default
       call CON_stop(NameSub // ' unknown TypeFile =' // trim(TypeFile))
    end select

    if(allocated(Param_I))   deallocate(Param_I)
    if(allocated(NameVar_I)) deallocate(NameVar_I)
    if(allocated(Coord_ID))  deallocate(Coord_ID)
    if(allocated( Var_IV))   deallocate(Var_IV)

  end subroutine save_plot_file

  !=========================================================================

  subroutine read_plot_file(NameFile, iUnitIn,         &
       TypeFileIn, StringHeaderOut,                    &
       nStepOut, TimeOut, nDimOut, nParamOut, nVarOut, &
       IsCartesianOut,                                 &
       n1Out, n2Out, n3Out, nOut_D,                    &
       ParamOut_I, NameVarOut,                         &
       CoordMinOut_D, CoordMaxOut_D,                   &
       CoordOut_DI,                                    &
       Coord1Out_I, Coord2Out_I, Coord3Out_I,          &
       CoordOut_I, CoordOut_DII, CoordOut_DIII,        &
       VarOut_I,  VarOut_II,  VarOut_III,              &
       VarOut_VI, VarOut_VII, VarOut_VIII,             &
       VarOut_IV, VarOut_IIV, VarOut_IIIV,             &
       iErrorOut)

    ! Both VarOut_VI and CoordOut_DI can be used in 1D, 2D, and 3D

    character(len=*),           intent(in) :: NameFile
    integer,          optional, intent(in) :: iUnitIn
    character(len=*), optional, intent(in) :: TypeFileIn
    character(len=*), optional, intent(out):: StringHeaderOut
    character(len=*), optional, intent(out):: NameVarOut
    real,             optional, intent(out):: TimeOut  
    integer,          optional, intent(out):: nStepOut
    integer,          optional, intent(out):: nDimOut   ! number of dimensions
    integer,          optional, intent(out):: nParamOut ! number of parameters
    integer,          optional, intent(out):: nVarOut   ! number of variables
    integer,          optional, intent(out):: n1Out, n2Out, n3Out ! grid size
    integer,          optional, intent(out):: nOut_D(:) ! grid size array
    logical,          optional, intent(out):: IsCartesianOut ! Cartesian grid?
    real,             optional, intent(out):: ParamOut_I(:)  ! parameters
    real,             optional, intent(out):: CoordMinOut_D(:)
    real,             optional, intent(out):: CoordMaxOut_D(:)
    real,             optional, intent(out):: CoordOut_DI(:,:) ! for 1D,2D,3D
    real,             optional, intent(out):: Coord1Out_I(:)
    real,             optional, intent(out):: Coord2Out_I(:)
    real,             optional, intent(out):: Coord3Out_I(:)
    real,             optional, intent(out):: CoordOut_I(:)          ! 1D
    real,             optional, intent(out):: CoordOut_DII(:,:,:)    ! 2D
    real,             optional, intent(out):: CoordOut_DIII(:,:,:,:) ! 3D
    real,             optional, intent(out):: VarOut_I(:)    ! variable  in 1D
    real,             optional, intent(out):: VarOut_II(:,:)       !        2D
    real,             optional, intent(out):: VarOut_III(:,:,:)    !        3D
    real,             optional, intent(out):: VarOut_VI(:,:) ! variables in 1D
    real,             optional, intent(out):: VarOut_VII(:,:,:)    !        2D
    real,             optional, intent(out):: VarOut_VIII(:,:,:,:) !        3D
    real,             optional, intent(out):: VarOut_IV(:,:)       !        1D
    real,             optional, intent(out):: VarOut_IIV(:,:,:)    !        2D
    real,             optional, intent(out):: VarOut_IIIV(:,:,:,:) !        3D
    integer,          optional, intent(out):: iErrorOut            ! I/O error

    integer            :: iUnit
    character(len=20)  :: TypeFile
    logical            :: DoReadHeader = .true.
    character(len=500) :: StringHeader
    character(len=500) :: NameVar
    integer            :: nStep, nDim, nParam, nVar, n1, n2, n3, n_D(MaxDim)
    real               :: Time, Coord
    real(Real4_)       :: Time4
    logical            :: IsCartesian
    real(Real4_), allocatable:: Param4_I(:), Coord4_ID(:,:), Var4_IV(:,:)
    real,         allocatable:: Param_I(:),  Coord_ID(:,:),  Var_IV(:,:)

    integer :: i, j, k, iDim, iVar, n

    ! Remember these values after reading header
    save :: nDim, nVar, n1, n2, n3, TypeFile, iUnit

    character(len=*), parameter:: NameSub = 'read_plot_file'
    !---------------------------------------------------------------------
    iUnit = UnitTmp_
    if(present(iUnitIn)) iUnit = iUnitIn

    TypeFile = 'ascii'
    if(present(TypeFileIn)) TypeFile = TypeFileIn

    if(present(iErrorOut)) iErrorOut = 0
    
    if(DoReadHeader) call read_header
    DoReadHeader = .false.

    ! No data is read. Leave file open !
    if(.not. (present(VarOut_I) .or. present(VarOut_II) &
         .or. present(VarOut_III) &
         .or. present(VarOut_VI) .or. present(VarOut_VII) &
         .or. present(VarOut_VIII).or. present(VarOut_IV)&
         .or. present(VarOut_IIV).or. present(VarOut_IIIV))) RETURN

    if((present(VarOut_I) .or. present(VarOut_II) .or. present(VarOut_III)) &
         .and. nVar /= 1)then
       write(*,*) NameSub,': the number of variables is ', nVar, &
            ' (larger than 1) in file ', NameFile
       call CON_stop(NameSub//' called with scalar variable argument')
    end if

    ! If data is read, next header needs to be read
    DoReadHeader = .true.

    ! Read coordinates and variables into suitable 2D arrays
    allocate(Coord_ID(n1*n2*n3, nDim), Var_IV(n1*n2*n3, nVar))
    select case(TypeFile)
    case('ascii', 'formatted')
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          read(iUnit, *, ERR=77, END=77) Coord_ID(n, :), Var_IV(n, :)
       end do; end do; end do

    case('real8')
       read(iUnit, ERR=77, END=77) Coord_ID
       do iVar = 1, nVar
          read(iUnit, ERR=77, END=77) Var_IV(:, iVar)
       end do

    case('real4')
       allocate(Coord4_ID(n1*n2*n3, nDim), Var4_IV(n1*n2*n3, nVar))
       read(iUnit, ERR=77, END=77) Coord4_ID
       Coord_ID = Coord4_ID
       do iVar = 1, nVar
          read(iUnit, ERR=77, END=77) Var4_IV(:, iVar)
       end do
       Var_IV = Var4_IV
       deallocate(Coord4_ID, Var4_IV)
    end select

    ! if iUnitIn is passed, keep file connected
    if(.not.present(iUnitIn)) close(iUnit) 

    if(present(CoordMinOut_D)) CoordMinOut_D(1:nDim) = minval(Coord_ID, DIM=1)
    if(present(CoordMaxOut_D)) CoordMaxOut_D(1:nDim) = maxval(Coord_ID, DIM=1)

    ! Fill in output coordinate arrays
    do iDim = 1, nDim
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          Coord = Coord_ID(n, iDim)
          if(present(CoordOut_DI))   CoordOut_DI(iDim,n)       = Coord
          if(present(CoordOut_I))    CoordOut_I(i)             = Coord
          if(present(CoordOut_DII))  CoordOut_DII(iDim,i,j)    = Coord
          if(present(CoordOut_DIII)) CoordOut_DIII(iDim,i,j,k) = Coord
          if(present(Coord1Out_I) .and. iDim==1 .and. j==1 .and. k==1) &
               Coord1Out_I(i) = Coord
          if(present(Coord2Out_I) .and. iDim==2 .and. i==1 .and. k==1) &
               Coord2Out_I(j) = Coord
          if(present(Coord3Out_I) .and. iDim==3 .and. i==1 .and. j==1) &
               Coord3Out_I(k) = Coord
       end do; end do; end do
    end do

    ! Fill in output variable arrays
    do iVar = 1, nVar
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          if(present(VarOut_I))    VarOut_I(n)             = Var_IV(n,iVar)
          if(present(VarOut_II))   VarOut_II(i,j)          = Var_IV(n,iVar)
          if(present(VarOut_III))  VarOut_III(i,j,k)       = Var_IV(n,iVar)
          if(present(VarOut_VI))   VarOut_VI(iVar,n)       = Var_IV(n,iVar)
          if(present(VarOut_VII))  VarOut_VII(iVar,i,j)    = Var_IV(n,iVar)
          if(present(VarOut_VIII)) VarOut_VIII(iVar,i,j,k) = Var_IV(n,iVar)
          if(present(VarOut_IV))   VarOut_IV(i,iVar)       = Var_IV(n,iVar)
          if(present(VarOut_IIV))  VarOut_IIV(i,j,iVar)    = Var_IV(n,iVar)
          if(present(VarOut_IIIV)) VarOut_IIIV(i,j,k,iVar) = Var_IV(n,iVar)

       end do; end do; end do
    end do

    deallocate(Coord_ID, Var_IV)

    RETURN

77  if(.not.present(iErrorOut)) call CON_stop(NameSub // &
         ' could not read data from file=' // trim(NameFile))

    iErrorOut = 3
    close(iUnit)
    if(allocated(Coord_ID))  deallocate(Coord_ID, Var_IV)
    if(allocated(Coord4_ID)) deallocate(Coord4_ID, Var4_IV)

  contains
    !==========================================================================
    subroutine read_header

      n_D = 1
      select case(TypeFile)
      case('ascii', 'formatted')
         open(iUnit, file=NameFile, status='old', ERR=66)

         read(iUnit, '(a)', ERR=77, END=77) StringHeader
         read(iUnit, *    , ERR=77, END=77) nStep, Time, nDim, nParam, nVar
         read(iUnit, *    , ERR=77, END=77) n_D(1:abs(nDim))
         if(nParam > 0)then
            allocate(Param_I(nParam))
            read(iUnit, *    , ERR=77, END=77) Param_I
         end if
         read(iUnit, '(a)', ERR=77, END=77) NameVar
      case('real8')
         open(iUnit, file=NameFile, status='old', form='unformatted', ERR=66)

         read(iUnit, ERR=77, END=77) StringHeader       
         read(iUnit, ERR=77, END=77) nStep, Time, nDim, nParam, nVar
         read(iUnit, ERR=77, END=77) n_D(1:abs(nDim))
         if(nParam > 0)then
            allocate(Param_I(nParam))
            read(iUnit, ERR=77, END=77) Param_I
         end if
         read(iUnit, ERR=77, END=77) NameVar
      case('real4')
         open(iUnit, file=NameFile, status='old', form='unformatted', ERR=66)

         read(iUnit, ERR=77, END=77) StringHeader
         read(iUnit, ERR=77, END=77) nStep, Time4, nDim, nParam, nVar
         Time = Time4
         read(iUnit, ERR=77, END=77) n_D(1:abs(nDim))
         if(nParam > 0)then
            allocate(Param_I(nParam), Param4_I(nParam))
            read(iUnit, ERR=77, END=77) Param4_I
            Param_I = Param4_I
         end if
         deallocate(Param4_I)
         read(iUnit, ERR=77, END=77) NameVar
      case default
         call CON_stop(NameSub // ' unknown TypeFile =' // trim(TypeFile))
      end select

      IsCartesian = nDim > 0
      nDim = abs(nDim)
      n1 = n_D(1); n2 = n_D(2); n3 = n_D(3) 

      if(present(StringHeaderOut)) StringHeaderOut = trim(StringHeader)
      if(present(NameVarOut))      NameVarOut      = trim(NameVar)
      if(present(TimeOut))         TimeOut         = Time
      if(present(nStepOut))        nStepOut        = nStep
      if(present(nDimOut))         nDimOut         = nDim
      if(present(nParamOut))       nParamOut       = nParam
      if(present(nVarOut))         nVarOut         = nVar
      if(present(n1Out))           n1Out           = n1
      if(present(n2Out))           n2Out           = n2
      if(present(n3Out))           n3Out           = n3
      if(present(nOut_D))          nOut_D(1:nDim)  = n_D(1:nDim)
      if(present(IsCartesianOut))  IsCartesianOut  = IsCartesian
      if(present(ParamOut_I) .and. nParam > 0) &
           ParamOut_I(1:nParam) = Param_I

      if(allocated(Param_I)) deallocate(Param_I)

      RETURN

66    if(.not.present(iErrorOut)) call CON_stop(NameSub // &
           ' could not open '//trim(TypeFile)//' file=' // trim(NameFile))
      
      iErrorOut = 1
      RETURN

77    if(.not.present(iErrorOut)) call CON_stop(NameSub // &
           ' could not read header from file=' // trim(NameFile))

      iErrorOut = 2
      close(iUnit)
      if(allocated(Param_I)) deallocate(Param_I)
      if(allocated(Param4_I)) deallocate(Param4_I)
      RETURN

    end subroutine read_header

  end subroutine read_plot_file

  !=========================================================================

  subroutine test_plot_file

    ! Set up a hydro shock tube initial condition on a 2D Cartesian grid
    ! Save plot file then read it and check consistency
    ! Do this multiple times with various settings

    character(len=*), parameter:: StringHeaderIn = "test_hd22"
    real,    parameter :: TimeIn = 25.0
    integer, parameter :: nStepIn = 10, nDimIn = 2, nParamIn = 2, nVarIn = 4
    integer, parameter :: n1In= 10, n2In = 2
    real,    parameter :: CoordMinIn_D(nDimIn) = (/ 0.5, -0.5 /)
    real,    parameter :: CoordMaxIn_D(nDimIn) = (/ 9.5,  0.5 /)
    real,    parameter :: ParamIn_I(nParamIn) = (/ 1.667, 2.5 /)
    character(len=*), parameter:: NameVarIn = "x y rho ux uy p gamma rbody"
    real    :: CoordIn_DII(nDimIn, n1In, n2In)
    real    :: VarIn_VII(nVarIn, n1In, n2In)
    real    :: CoordIn_DIII(nDimIn, n1In, 1, n2In)
    real    :: VarIn_VIII(nVarIn, n1In, 1, n2In)
    real    :: VarIn_IIV(n1In, n2In, nVarIn)

    ! Do tests with ascii/real8/real4 files, 
    ! Cartesian/non-Cartesian coordinates
    ! 2D/3D input arrays
    integer, parameter:: nTest = 15
    character(len=5)  :: TypeFileIn_I(nTest) = &
         (/ 'ascii', 'real8', 'real4', 'ascii', 'real8', 'real4', &
         'ascii', 'real8', 'real4', 'ascii', 'real8', 'real4', &
         'ascii', 'real8', 'real4' /)
    logical           :: IsCartesianIn_I(nTest) = &
         (/ .true.,   .true., .true.,  .false.,   .false., .false.,&
         .true.,   .true., .true.,  .false.,   .false., .false.,&
         .true., .false., .false. /)

    ! Input and output of tests
    character(len=80)    :: NameFile
    character(len=20)    :: TypeFileIn
    character(len=100)   :: StringHeaderOut
    real                 :: TimeOut
    integer              :: nStepOut, nDimOut, nParamOut, nVarOut
    integer              :: nOut_D(3)
    real                 :: ParamOut_I(100)
    character(len=100)   :: NameVarOut
    logical              :: IsCartesianIn, IsCartesianOut
    real                 :: CoordMinOut_D(nDimIn), CoordMaxOut_D(nDimIn)
    real                 :: Coord1Out_I(n1In), Coord2Out_I(n2In)
    real                 :: CoordOut_DII(nDimIn, n1In, n2In)
    real                 :: VarOut_VII(nVarIn, n1In, n2In)

    ! Tolerance for errors
    real :: Eps

    ! Indexes
    integer :: i, j, iTest

    character(len=*), parameter:: NameSub = 'test_plot_file'
    !----------------------------------------------------------------------

    ! Initialize coordinates and variables: shock tube on a 2D uniform grid
    do j = 1, n2In; do i = 1, n1In
       CoordIn_DII(1, i, j) = CoordMinIn_D(1) &
            + (i-1)*((CoordMaxIn_D(1)-CoordMinIn_D(1))/(n1In - 1))
       CoordIn_DII(2, i, j) = CoordMinIn_D(2) &
            + (j-1)*((CoordMaxIn_D(2)-CoordMinIn_D(2))/(n2In - 1))
       CoordIn_DIII(1, i, 1, j) = CoordIn_DII(1,i,j)
       CoordIn_DIII(2, i, 1, j) = CoordIn_DII(2,i,j)

       if(i <= n1In/2)then
          VarIn_VII(:, i, j) = (/ 1.0, 0.0, 0.0, 1.0 /)
          VarIn_IIV(i, j, :) = (/ 1.0, 0.0, 0.0, 1.0 /)
          VarIn_VIII(:, i, 1, j) = (/ 1.0, 0.0, 0.0, 1.0 /)
       else
          VarIn_VII(:, i, j) = (/ 0.1, 0.0, 0.0, 0.125 /)
          VarIn_IIV(i, j, :) = (/ 0.1, 0.0, 0.0, 0.125 /)
          VarIn_VIII(:, i, 1, j) = (/ 0.1, 0.0, 0.0, 0.125 /)  
       end if
    end do; end do

    ! Test ascii, real8 and real4 files
    do iTest = 1, nTest 
       write(NameFile, '(a,i2.2,a)') 'test_plot_file',iTest,'.out'
       write(*,*) NameSub, ' writing file=', trim(NameFile)

       TypeFileIn    = TypeFileIn_I(iTest)
       IsCartesianIn = IsCartesianIn_I(iTest)

       if(TypeFileIn == 'real4')then
          Eps = 1e-5
       else
          Eps = 1e-12
       end if

       ! Test saving it
       select case(iTest)
       case(1)
          ! Use coordinate ranges
          call save_plot_file(NameFile,      &
               TypeFileIn     = TypeFileIn,     &
               StringHeaderIn = StringHeaderIn, &
               nStepIn        = nStepIn,        &
               TimeIn         = TimeIn,         &
               ParamIn_I      = ParamIn_I,      &
               NameVarIn      = NameVarIn,      &
               IsCartesianIn  = IsCartesianIn,  &
               CoordMinIn_D   = CoordMinIn_D,   &
               CoordMaxIn_D   = CoordMaxIn_D,   &
               VarIn_IIV      = VarIn_IIV)
       case(2)
          ! Use 1D coordinate arrays
          call save_plot_file(NameFile,      &
               TypeFileIn     = TypeFileIn,     &
               StringHeaderIn = StringHeaderIn, &
               nStepIn        = nStepIn,        &
               TimeIn         = TimeIn,         &
               ParamIn_I      = ParamIn_I,      &
               NameVarIn      = NameVarIn,      &
               IsCartesianIn  = IsCartesianIn,  &
               Coord1In_I     = CoordIn_DII(1,:,1), &
               Coord2In_I     = CoordIn_DII(2,1,:), &
               VarIn_VII      = VarIn_VII)
       case(3:6)
          ! Use full coordinate array
          call save_plot_file(NameFile,      &
               TypeFileIn     = TypeFileIn,     &
               StringHeaderIn = StringHeaderIn, &
               nStepIn        = nStepIn,        &
               TimeIn         = TimeIn,         &
               ParamIn_I      = ParamIn_I,      &
               NameVarIn      = NameVarIn,      &
               IsCartesianIn  = IsCartesianIn,  &
               CoordIn_DII    = CoordIn_DII,    &
               VarIn_VII      = VarIn_VII)
       case default
          ! Test 3D input array
          ! Use full coordinate array
          call save_plot_file(NameFile,      &
               TypeFileIn     = TypeFileIn,     &
               StringHeaderIn = StringHeaderIn, &
               nStepIn        = nStepIn,        &
               TimeIn         = TimeIn,         &
               ParamIn_I      = ParamIn_I,      &
               NameVarIn      = NameVarIn,      &
               IsCartesianIn  = IsCartesianIn,  &
               nDimIn         = nDimIn,         &
               CoordIn_DIII   = CoordIn_DIII,   &
               VarIn_VIII     = VarIn_VIII)

       end select

       if(iTest == 13 .or. iTest == 14 .or. iTest == 15)then
          ! Read header and data separately
          call read_plot_file(NameFile,        &
               TypeFileIn      = TypeFileIn,      &
               StringHeaderOut = StringHeaderOut, &
               nStepOut        = nStepOut,        &
               TimeOut         = TimeOut,         &
               nDimOut         = nDimOut,         &
               nParamOut       = nParamOut,       &
               nVarOut         = nVarOut,         &
               ParamOut_I      = ParamOut_I,      &
               NameVarOut      = NameVarOut,      &
               IsCartesianOut  = IsCartesianOut)
          call read_plot_file(NameFile,        &
               TypeFileIn      = TypeFileIn,      &
               CoordOut_DII    = CoordOut_DII,    &
               Coord1Out_I     = Coord1Out_I,     &
               Coord2Out_I     = Coord2Out_I,     &
               CoordMinOut_D   = CoordMinOut_D,   &
               CoordMaxOut_D   = CoordMaxOut_D,   &
               VarOut_VII      = VarOut_VII)
       else
          call read_plot_file(NameFile,        &
               TypeFileIn      = TypeFileIn,      &
               StringHeaderOut = StringHeaderOut, &
               nStepOut        = nStepOut,        &
               TimeOut         = TimeOut,         &
               nDimOut         = nDimOut,         &
               nParamOut       = nParamOut,       &
               nVarOut         = nVarOut,         &
               ParamOut_I      = ParamOut_I,      &
               NameVarOut      = NameVarOut,      &
               IsCartesianOut  = IsCartesianOut,  &
               CoordOut_DII    = CoordOut_DII,    &
               Coord1Out_I     = Coord1Out_I,     &
               Coord2Out_I     = Coord2Out_I,     &
               CoordMinOut_D   = CoordMinOut_D,   &
               CoordMaxOut_D   = CoordMaxOut_D,   &
               VarOut_VII      = VarOut_VII)

       end if

       if(nStepOut /= nStepIn)then
          write(*,*)'nStepIn=', nStepIn,' nStepOut=', nStepOut
          call CON_stop(NameSub)
       end if

       if(abs(TimeOut - TimeIn) > Eps)then
          write(*,*)'TimeIn=', TimeIn,' TimeOut=', TimeOut
          call CON_stop(NameSub)
       end if

       if(nDimOut /= nDimIn)then
          write(*,*)'nDimIn=', nDimIn,' nDimOut=', nDimOut
          call CON_stop(NameSub)
       end if

       if(nParamOut /= nParamIn)then
          write(*,*)'nParamIn=', nParamIn,' nParamOut=', nParamOut
          call CON_stop(NameSub)
       end if

       if(nVarOut /= nVarIn)then
          write(*,*)'nVarIn=', nVarIn,' nVarOut=', nVarOut
          call CON_stop(NameSub)
       end if

       if(any(abs(ParamOut_I(1:nParamIn) - ParamIn_I) > Eps))then
          write(*,*)'ParamIn=', ParamIn_I,' ParamOut=', ParamOut_I(1:nParamIn)
          call CON_stop(NameSub)
       end if

       if(IsCartesianOut .neqv. IsCartesianIn)then
          write(*,*)'IsCartesianIn, Out=', IsCartesianIn, IsCartesianOut
          call CON_stop(NameSub)
       end if

       if(NameVarOut /= NameVarIn)then
          write(*,*)'NameVarIn=', NameVarIn,' NameVarOut=', NameVarOut
          call CON_stop(NameSub)
       end if

       !To simplify, replace the 3D input array with 2D 
       if(iTest > 6)then
          CoordIn_DII = CoordIn_DIII(:,:,1,:)
          VarIn_VII = VarIn_VIII(:,:,1,:)
       end if
       do j = 1, n2In; do i = 1, n1In
          if(any(abs(CoordIn_DII(:,i,j) - CoordOut_DII(:,i,j)) > Eps))then
             write(*,*)'i,j=', i, j
             write(*,*)'CoordIn =', CoordIn_DII(:,i,j)
             write(*,*)'CoordOut=', CoordOut_DII(:,i,j)
             call CON_stop(NameSub)
          end if
          if(abs(CoordIn_DII(1,i,j) - Coord1Out_I(i)) > Eps )then
             write(*,*)'i,j=', i, j
             write(*,*)'CoordIn(1)=', CoordIn_DII(1,i,j)
             write(*,*)'Coord1Out =', Coord1Out_I(i)
             call CON_stop(NameSub)
          end if
          if(abs(CoordIn_DII(2,i,j) - Coord2Out_I(j)) > Eps )then
             write(*,*)'i,j=', i, j
             write(*,*)'CoordIn(2)=', CoordIn_DII(2,i,j)
             write(*,*)'Coord2Out =', Coord2Out_I(j)
             call CON_stop(NameSub)
          end if
          if(any(abs(VarIn_VII(:,i,j) - VarOut_VII(:,i,j)) > Eps))then
             write(*,*)'i,j=', i, j
             write(*,*)'VarIn =', VarIn_VII(:,i,j)
             write(*,*)'VarOut=', VarOut_VII(:,i,j)
             call CON_stop(NameSub)
          end if
          if(any(abs(VarIn_IIV(i,j,:) - VarOut_VII(:,i,j)) > Eps))then
             write(*,*)'i,j=', i, j
             write(*,*)'VarIn =', VarIn_IIV(i,j,:)
             write(*,*)'VarOut=', VarOut_VII(:,i,j)
             call CON_stop(NameSub)
          end if
       end do; end do

       if(abs(CoordMinOut_D(1) -  minval(CoordIn_DII(1,:,:))) >Eps)then
          write(*,*)'CoordMinOut_D(1)     =',CoordMinOut_D(1)
          write(*,*)'minval(CoordIn_DII(1)=',minval(CoordIn_DII(1,:,:))
          call CON_stop(NameSub)
       end if
       if(abs(CoordMinOut_D(2) -  minval(CoordIn_DII(2,:,:))) > Eps)then
          write(*,*)'CoordMinOut_D(2)     =',CoordMinOut_D(2)
          write(*,*)'minval(CoordIn_DII(2)=',minval(CoordIn_DII(2,:,:))
          call CON_stop(NameSub)
       end if
       if(abs(CoordMaxOut_D(1) -  maxval(CoordIn_DII(1,:,:))) > Eps )then
          write(*,*)'CoordMaxOut_D(1)     =',CoordMaxOut_D(1)
          write(*,*)'maxval(CoordIn_DII(1)=',maxval(CoordIn_DII(1,:,:))
          call CON_stop(NameSub)
       end if
       if(abs(CoordMaxOut_D(2) - maxval(CoordIn_DII(2,:,:))) >Eps)then
          write(*,*)'CoordMaxOut_D(2)     =',CoordMaxOut_D(2)
          write(*,*)'maxval(CoordIn_DII(2)=',maxval(CoordIn_DII(2,:,:))
          call CON_stop(NameSub)
       end if

    end do

    ! Test using defaults for 2D input array
    NameFile = 'test_plot_file16.out'       
    call save_plot_file(NameFile, VarIn_VII=VarIn_VII, CoordIn_DII=CoordIn_DII)

    call read_plot_file(NameFile, &
         StringHeaderOut=StringHeaderOut, NameVarOut=NameVarOut, &         
         nDimOut=nDimOut, nVarOut=nVarOut, nParamOut=nParamOut, &
         IsCartesianOut=IsCartesianOut, nOut_D=nOut_D)

    if( StringHeaderOut /= 'No header info')call CON_stop(NameSub // &
         ' incorrect value for default StringHeaderOut='//StringHeaderOut)
    if( NameVarOut /= 'x1 x2 v01 v02 v03 v04') call CON_stop(NameSub // &
         ' incorrect value for default NameVarOut='//NameVarOut)
    if(.not.IsCartesianOut) call CON_stop(NameSub // &
         ' incorrect value for IsCartesianOut: should be true')
    if( any(nOut_D(1:nDimOut) /= (/ n1In, n2In /)) ) then
       write(*,*) 'n1In, n2In=', n1In, n2In
       write(*,*) 'nOut_D    =',nOut_D
       call CON_stop(NameSub)
    end if

    ! Test using defaults for 3D input array
    NameFile = 'test_plot_file17.out'       
    call save_plot_file(NameFile,nDimIn = nDimIn, VarIn_VIII = VarIn_VIII,&
         CoordIn_DIII = CoordIn_DIII)

    ! Read header info
    call read_plot_file(NameFile, &
         StringHeaderOut = StringHeaderOut, NameVarOut = NameVarOut, & 
         nDimOut = nDimOut, nVarOut = nVarOut, nParamOut = nParamOut, &
         IsCartesianOut = IsCartesianOut, nOut_D = nOut_D)

    if( StringHeaderOut /= 'No header info')call CON_stop(NameSub // &
         ' incorrect value for default StringHeaderOut='//StringHeaderOut)
    if( NameVarOut /= 'x1 x2 v01 v02 v03 v04') call CON_stop(NameSub // &
         ' incorrect value for default NameVarOut='//NameVarOut)
    if(.not.IsCartesianOut) call CON_stop(NameSub // &
         ' incorrect value for IsCartesianOut: should be true')
    if( any(nOut_D(1:nDimOut) /= (/ n1In, n2In /)) ) then
       write(*,*) 'n1In, n2In=', n1In, n2In
       write(*,*) 'nOut_D    =',nOut_D
       call CON_stop(NameSub)
    end if


    ! Now that we have the dimensions, we could allocate coordinate and
    ! variable arrays and read them

  end subroutine test_plot_file

end module ModPlotFile
