!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHdf5Utils

  use hdf5
  use ModMpiOrig

  implicit none
  private!except
  save

! Explanation of BATL File Format:
!     
! Arrays:
!     
! Axis Labels
!   Axis label names.  VisIt/HDF5 doesn’t like character strings that
!   aren’t null terminated so make sure to fill any trailing spaces with
!   null characters.  This can be done with the pad_string_with_null
!   subroutine in ModHdf5Utils.
! 
! plotVarUnits
!   The unit names for all the plot variables.  .  VisIt/HDF5 doesn’t like
!   character strings that aren’t null terminated so make sure to fill any
!   trailing spaces with null characters.  This can be done with the
!   pad_string_with_null subroutine in ModHdf5Utils.
! 
! plotVarNames
!   The names of all the plot variables.  .  VisIt/HDF5 doesn’t like
!   character strings that aren’t null terminated so make sure to fill any
!   trailing spaces with null characters.  This can be done with the
!   pad_string_with_null subroutine in ModHdf5Utils
! 
! 
! Integer Plot Metadata
!   Passes necessary/potentially useful integer metadata to VisIt in the
!   following order:
!   1.      File Format Version – Currently 1
!   2.      Time-step number
!   3.      Number of plot dimensions
!   4.      Number of AMR dimensions
!   5.      Total number of blocks used
!   6.      Number of processors
!   7.      Number of refinement levels
!   8.      Cells per block i
!   9.      Cells per block j
!   10.     Cells per block k
!   11.     Geometry type: See plot geometry handles in ModHdf5Utils.
!     Basically anything less than 2 is Cartesian, anything .ge. 2 is
!     non-Cartesian.  12 is a topologically 2D surface in 3D space.
!   12.     1 if periodic in i dimension else 0
!   13.     1 if periodic in j dimension else 0
!   14.     1 if periodic in k dimension else 0
!   15.     0 if you have Morton tree data else 1
!   16.     Number of plot variables
! 
! Integer Sim Metadata (Optional)
!   You could put anything you want here.  It won’t affect VisIt at all.
!   Currently, for plots containing this array the data is as follows:
!   1.      Cells per block i
!   2.      Cells per block j
!   3.      Cells per block k
!   4.      Number of dimensions
!   5.      Number of AMR dimensions
!   6.      Number of processors
!   7.      Number of refinement levels
!   8.      Time step number
! 
! Real Plot Metadata:
!   1.      Simulation time
!   2.      Xmin
!   3.      Xmax
!   4.      Ymin
!   5.      Ymax
!   6.      Zmin
!   7.      Zmax
! 
! Real Simulation Metadata (optional)
!   Currently the same as real plot metadata except that it shows the
!   size of the computational domain rather than the plot boundaries.
! 
! MinLogicalExtents
!   Minimum block index.  See p.210 of “GettingDataIntoVisit2.0.0” on
!   VisIt’s website for an explanation of VisIt zone indexing.  This is
!   used with the info given in the Integer Plot Metadata array to give
!   visit the indices for each block.
! 
! bounding box
!   Minimum and Maximum coordinates in Cartesian space
! 
! Nodes(X,Y,Z)
!   (X,Y,Z) Cartesian grid node coordinates for non-Cartesian plot files.
!   If the plot is 2d you only need to use NodesX and NodesY.
! 
! Processor Number
!   Processor number for each block.  This is required for VisIt 2.5.1 but
!   will be optional in VisIt 2.5.2.
! 
! refine level
!   The refinement level for each block.  Required for all files in VisIt
!   2.5.1.  In VisIt 2.5.2 it may be left out if all blocks are at the
!   same refinement level.
! 
! iMortonNode_A
!   Required if the presence of a Morton curve is indicated in the integer
!   plot metadata, otherwise not necessary.  Gives the Morton index of the
!   block.
! 
! coordinates
!   Required in VisIt 2.5.1, only required when iMortonNode_A is present
!   in VisIt 2.5.2.  Gives the block center location for drawing the
!   Morton curve.
! 
! (Variable Name)
!   Cell centered variable data.  Has minimum and maximum values for the
!   dataset stored in attributes named minimum and maximum.
! 
! (Variable Name_Ext) (optional but sometimes useful for VisIt)
!   Minimum and maximum values for each block of a plot variable.


  public:: save_hdf5_file
  public:: open_hdf5_file
  public:: write_hdf5_data
  public:: write_hdf5_attribute
  public:: pad_string_with_null
  public:: close_hdf5_file

  logical :: DoCollectiveWrite
    
  !plot geometry type handles 
  integer, public, parameter :: CartesianPlot_=0
  integer, public, parameter :: RzPlot_=1
  integer, public, parameter :: RoundCubePlot_=2
  integer, public, parameter :: CylindricalPlot_=3
  integer, public, parameter :: CylindricalLnrPlot_=4
  integer, public, parameter :: CylindricalGenrPlot_=5
  integer, public, parameter :: SphericalPlot_=6
  integer, public, parameter :: SphericalLnrPlot_=7
  integer, public, parameter :: SphericalGenrPlot_=8
  integer, public, parameter :: rLonLatPlot_=9
  integer, public, parameter :: rLonLatLnrPlot_=10
  integer, public, parameter :: rLonLatGenrPlot_=11
  integer, public, parameter :: SphShellPlot_=12

  integer, public, parameter:: lNameVar = 10
  integer, public, parameter:: lNameH5 = lNameVar + 1

contains 
  !===========================================================================
  subroutine save_hdf5_file(FileName,TypePosition, TypeStatus, StringHeader,&
            nStep, NumberOfBlocksUsed, Time, nDimOut, nParam, nVar,&
            n_D, NameVar, NameUnits, MinimumBlockIjk, XYZMinMax, PlotVarBlk,&
            iComm, CoordMin, CoordMax)
    use ModUtilities, only: split_string

    integer, intent(in) :: nDimOut, nParam, nVar, n_D(3)
    integer, intent(in) :: nStep, NumberOfBlocksUsed
    character (len=*), intent(in) :: FileName
    character (len=*), intent(in) :: TypePosition, NameVar(nVar)
    character (len=*), intent(in) ::  NameUnits
    character (len=10), intent(in) :: TypeStatus
    character (len=500), intent(in) :: StringHeader
    integer, intent(in) :: MinimumBlockIjk(:,:)
    real, intent(in) :: Time, PlotVarBlk(:,:,:,:,:), XYZMinMax(:,:,:)

    integer, optional, intent(in):: iComm
    real, optional, intent(in) :: CoordMin(:), CoordMax(:)

    character (len=501) :: HeaderString
    integer (HID_T) :: FileID, iInteger4
    integer(HSIZE_T) :: iDimension1D(1), iInteger8
    integer :: IntegerMetaData(16), iLen, iVar, iProc, nProc, iCommOpen, iError
    integer :: LengthOfString, numCells,i,j,k,n, iBlk
    integer, parameter :: FFV = 1

    real :: RealMetaData(7),varMin, varMax,xMin,xMax,yMin,yMax,zMin,zMax
    real, allocatable :: coordinates(:,:)
    integer, allocatable :: procNumAndLevel(:)
    integer(HSIZE_T), parameter :: lnameh5 = 11
    integer(HSIZE_T) :: CellsPerBlock(3), nBlkUsed, nPlotDim, iData,AnotheriInteger8
    character (len=lnameh5), allocatable  :: UnknownNameArray(:), UnknownNameArray2(:), UnknownNameArray3(:)

    !-------------------------------------------------------------------------
    if(present(iComm))then
       iCommOpen = iComm
       call MPI_comm_size(iComm, nProc, iError)
       call MPI_comm_rank(iComm, iProc, iError)
    else
       ! Serial write is the default
       iCommOpen = MPI_COMM_SELF
       iProc = 0
       nProc = 1
    end if


    nBlkUsed = NumberOfBlocksUsed
    CellsPerBlock(1:3) = n_D(1:3)
    nPlotDim = nDimOut
    
    FileID = -1
    call open_hdf5_file(FileID, fileName, iCommOpen)
    if(FileID == -1) then
       write (*,*)  "Error: unable to initialize file"
       call CON_stop("unable to initialize hdf5 file")
    end if

    xMin = 0; yMin = 0; zMin = 0
    if (present(CoordMin)) then
        xMin = CoordMin(1)
        if (nDimOut .ge. 2) then
            yMin = CoordMin(2)
            if (nDimOut == 3) &
                zMin = CoordMin(3)
        end if
    end if  ! else make a grid?
    
    xMax = 0; yMax = 0; zMax = 0
    if (present(CoordMax)) then
        xMax = CoordMax(1)
        if (nDimOut .ge. 2) then
            yMax = CoordMax(2)
            if (nDimOut == 3) &
                zMax = CoordMax(3)
        end if
    end if  ! else make a grid

 
    
    allocate(UnknownNameArray(nVar))
    call split_string(NameUnits,nVar,UnknownNameArray(1:nVar), i)
    iInteger4 = lNameh5
    call pad_string_with_null(nVar, iInteger4, UnknownNameArray, UnknownNameArray)
    iInteger8=nVar
    call write_hdf5_data(FileID, "plotVarUnits", 1,  (/iInteger8/),&
    CharacterData=UnknownNameArray, nStringChars=lNameH5)

    call pad_string_with_null(nVar, iInteger4, NameVar, UnknownNameArray)
    iInteger8=nVar
    call write_hdf5_data(FileID, "plotVarNames", 1,  (/iInteger8/),&
    CharacterData=UnknownNameArray, nStringChars=lNameH5)


   
         
    do iVar = 1, nVar
       varMin = minVal(PlotVarBlk(:,:,:,:,iVar))
       varMax = maxVal(PlotVarBlk(:,:,:,:,iVar))

        call write_hdf5_data(FileID, UnknownNameArray(iVar), 4, (/CellsPerBlock(1), CellsPerBlock(2),&
            CellsPerBlock(3),nBlkUsed/), Rank4RealData=PlotVarBlk(:,:,:,:,iVar),&
            RealAttribute1=VarMin, RealAttribute2=VarMax, &
            NameRealAttribute1="minimum", NameRealAttribute2="maximum")
    end do
  deallocate(UnknownNameArray)

        iInteger8 = 2
       call  write_hdf5_data(FileID, "bounding box", 3,&
      (/iInteger8, nPlotDim,nBlkUsed/),Rank3RealData=XYZMinMax(:,1:nPlotDim,:))

    allocate(coordinates(nPlotDim, nBlkUsed))
    do iBlk=1,nBlkUsed
       coordinates(1:nPlotDim, iBlk) = .5*(XYZMinMax(1,1:nPlotDim, iBlk) + XYZMinMax(2,1:nPlotDim, iBlk))
    end do
       call  write_hdf5_data(FileID, "coordinates", 2,&
      (/nPlotDim,nBlkUsed/),Rank2RealData=coordinates)

    deallocate(coordinates)
    
    allocate(UnknownNameArray(nPlotDim)) 
    UnknownNameArray(1) = "Y-Axis"
    if (nPlotDim .GE. 2) then
    UnknownNameArray(2) = "Z-Axis"   
    if (nPlotDim == 3) &
        UnknownNameArray(3) = "Z-Axis"
    end if
    iInteger4=lNameh5
    call pad_string_with_null(nVar,iInteger4 , UnknownNameArray, UnknownNameArray)
    iInteger8=nPlotDim
    call write_hdf5_data(FileID, "Axis Labels", 1,  (/iInteger8/),&
    CharacterData=UnknownNameArray, nStringChars=lNameH5)
    deallocate(UnknownNameArray)

    iLen = len(trim(StringHeader)) + 1
    call pad_string_with_null(1, iLen, StringHeader, HeaderString)
    AnotheriInteger8 = 1
    call write_hdf5_data(FileID, "Header", 1,  (/AnotheriInteger8/),&
    CharacterData=(/HeaderString/), nStringChars=iInteger8)

    call  write_hdf5_data(FileID, "MinLogicalExtents", 2,&
      (/nPlotDim,nBlkUsed/),Rank2IntegerData=MinimumBlockIjk)

    !As of VisIt 2.5.2 Processor Number and refine level are not required
    !by the plugin.  They are here so older versions of VisIt work
    allocate(ProcnumAndLevel(nBlkUsed))
    ProcnumAndLevel = 1
     call  write_hdf5_data(FileID, "refine level", 1, (/nBlkUsed/),&
        Rank1IntegerData=ProcNumAndLevel)
    procNumAndLevel = iProc
     call  write_hdf5_data(FileID, "Processor Number", 1, (/nBlkUsed/),&
        Rank1IntegerData=ProcNumAndLevel)
    deallocate(ProcnumAndLevel)
    iData = 1
    !    AttributeName(1) = 'Simulation Time'
    RealMetaData(iData) = Time
    iData = iData + 1
    RealMetaData(iData) = xMin
    iData = iData + 1
    RealMetaData(iData) = xMax
    iData = iData + 1
    RealMetaData(iData) = yMin
    iData = iData + 1
    RealMetaData(iData) = yMax
    iData = iData + 1
    if(nDimOut < 3) then
       RealMetaData(iData) = 0!coordMin(3)
       iData = iData + 1
       RealMetaData(iData) = 0!coordMax(3)   
    else
       RealMetaData(iData) = zMin
       iData = iData + 1
       RealMetaData(iData) = zMax
    end if

    !-------------------------------------------------------------------
    !write the real Metadata
    call  write_hdf5_data(FileID, "Real Plot Metadata", 1, (/iData/),&
        Rank1RealData=RealMetaData)

    iData = 1

    IntegerMetaData(iData) = FFV
    !    AttributeName(1) = 'File Format Version'
    iData = iData + 1
    IntegerMetaData(iData) = nStep 
    !    AttributeName(2) = "Time Step"
    iData = iData + 1
    IntegerMetaData(iData) = nDimOut
    !    AttributeName(3) = 'nDim'
    iData = iData + 1
    IntegerMetaData(iData) = 0
    !    AttributeName(3) = 'nDimAMR'
    iData = iData + 1
    IntegerMetaData(iData) = nBlkUsed
    !    AttributeName(4) = 'globalNumBlocks'
    iData = iData + 1    
    IntegerMetaData(iData) = nProc
    !    AttributeName(5) = 'numProcessors'
    iData = iData + 1

    IntegerMetaData(iData) = 1 ! no AMR levels
    !    AttributeName(7) = 'nLevel'
    iData = iData + 1
    IntegerMetaData(iData) = n_D(1)
    iData = iData + 1
    IntegerMetaData(iData) = n_D(2)
    iData = iData + 1
    IntegerMetaData(iData) = n_D(3)
    iData = iData + 1

    !This only does cartesian for now, all that is needid
    !for non-cartesian plots is node locations and a change
    !in this metadata item
    IntegerMetaData(iData) = CartesianPlot_
    iData = iData + 1
    !as of 2/3/2012 this is not implimented in the plugin but it probably
    !should be in the future
    do i = 1, 3
       IntegerMetaData(iData) = 0 ! none of the axis are periodic
       iData = iData + 1
    end do

    integerMetaData(iData) = 1 !Tells VisIt that this is a cut file, which to VisIt means that
    ! there is no morton curve index to be read.
    iData = iData + 1
    integerMetaData(iData) = nVar

    !-------------------------------------------------------------------
    !write the integer Metadata
    call  write_hdf5_data(FileID, "Integer Plot Metadata", 1, (/iData/),&
        Rank1IntegerData=IntegerMetaData)
    call close_hdf5_file(FileID)
  end subroutine save_hdf5_file
  
  !======================================================================
  !=====================================================================
  subroutine open_hdf5_file(FileID, Filename, iComm)

    character (len=80), intent(in) :: Filename
    integer :: iErrorHdf, AccessTemplate
    integer, optional, intent(in) :: iComm
    integer, intent(inout) :: FileID
    integer :: iHdfMajor, iHdfMinor, iHdfRelease

    call h5open_f(iErrorHdf)                    

    !create MPI info Object

    !Create file access propertty list
    call h5pcreate_f(H5P_FILE_ACCESS_F, AccessTemplate, iErrorHdf)
    if(present(iComm)) then
         CALL h5pset_fapl_mpio_f(AccessTemplate, iComm, MPI_INFO_NULL, iErrorHdf)

        !determine if we want to do collective writes. Collective write
        !is faster but all processors must call h5dwrite, even in cuts where 
        !some processors write no data.  In newer versions of hdf5 a null
        !write can be called but attempting to do so on older versions will
        !cause the code to crash. We should be able to do collective for any
        !hdf5 1.8.x where all procs write data but that didn't work for me 
        !for some reason so for now collective write mode is restricted to
        !hdf5 version 1.8.8+
        !     write (*,*) "write plot iErrorHdf 3"v
        call h5get_libversion_f(iHdfMajor, iHdfMinor, iHdfRelease, iErrorHdf)   
        if (iHdfMajor > 1) then
           DoCollectiveWrite = .true.
        else if (iHdfMajor == 1 .and. iHdfMinor > 8) then
           DoCollectiveWrite = .true.
        else if (iHdfMajor == 1 .and. iHdfMinor == 8) then
           if (iHdfRelease .ge. 8 ) then!.or. allProcsWrite) then
              DoCollectiveWrite = .true.
           else
              DoCollectiveWrite = .false.
           end if
        else 
           DoCollectiveWrite = .false.
        end if
    end if

    ! Create the file Collectively.

    CALL h5fcreate_f(Filename, H5F_ACC_TRUNC_F, FileID, iErrorHdf,&
         access_prp = AccessTemplate)
    CALL h5pclose_f(AccessTemplate, iErrorHdf)
    if (iErrorHdf == -1) &
         call CON_stop(&
         "iErrorHdf in subroutine hdf5_init_file. Error marker 1")

    if (iErrorHdf == -1) FileID = -1

  end subroutine open_hdf5_file

  !=====================================================================
  !=====================================================================
  subroutine write_hdf5_data(FileID, DatasetName, DatasetRank,  nDatasetDimension,&
    nOffsetLocal, nBlocksLocalIn, CoordArray, nCellsLocalIn, Rank4RealData, Rank3RealData, Rank2RealData,&
    Rank1RealData, Rank4IntegerData, Rank3IntegerData, Rank2IntegerData, Rank1IntegerData,&
    CharacterData, nStringChars, RealAttribute1,RealAttribute2,NameRealAttribute1,NameRealAttribute2,&
    IntegerAttribute1,IntegerAttribute2,NameIntegerAttribute1,NameIntegerAttribute2)

    character (len=*), intent(in) :: DatasetName

    integer(HID_T), intent(in) :: FileID
    integer, intent(in) :: DatasetRank

    integer(HSIZE_T), intent(in)  :: nDatasetDimension(DatasetRank)

    integer, optional, intent(in) :: nOffsetLocal
    integer, optional, intent(in) :: nBlocksLocalIn
    integer(HSIZE_T), optional, intent(in) :: CoordArray(:,:)
    integer(HSIZE_T), optional, intent(in) :: nCellsLocalIn
    integer(HSIZE_T), optional, intent(in) :: nStringChars

    ! output data
    real, optional, intent(in) :: Rank4RealData(:,:,:,:)
    real, optional, intent(in) :: Rank3RealData(:,:,:)
    real, optional, intent(in) :: Rank2RealData(:,:)
    real, optional, intent(in) :: Rank1RealData(:)
    integer, optional, intent(in) :: Rank4IntegerData(:,:,:,:)
    integer, optional, intent(in) :: Rank3IntegerData(:,:,:)
    integer, optional, intent(in) :: Rank2IntegerData(:,:)
    integer, optional, intent(in) :: Rank1IntegerData(:)
    character (len=*), optional, intent(in) :: CharacterData(:)

    real, optional, intent(in) :: RealAttribute1
    real, optional, intent(in) :: RealAttribute2
    character (len=*), optional, intent(in) :: NameRealAttribute1
    character (len=*), optional, intent(in) :: NameRealAttribute2
    
    integer, optional, intent(in) :: IntegerAttribute1
    integer, optional, intent(in) :: IntegerAttribute2
    character (len=*), optional, intent(in) :: NameIntegerAttribute1
    character (len=*), optional, intent(in) :: NameIntegerAttribute2


    integer :: iErrorHdf, nBlocksLocal
    integer(HID_T) :: DatasetID, PropertyListID
    integer(HID_T) :: DataSpaceId
    integer(HID_T) :: MemorySpaceId
    integer(HSIZE_T) :: nStart(DatasetRank)
    integer(HSIZE_T), allocatable :: nCount(:)
    integer(HID_T) :: DataType
    integer(HSIZE_T) :: iOneOrZero, nCellsLocal

    !--------------------------------------------------------------------------

    if(present(nBlocksLocalIn)) then
        nBlocksLocal = nBlocksLocalIn
        allocate(nCount(DatasetRank))
    else
        nBlocksLocal = nDatasetDimension(DatasetRank)
    end if
    if (present(nCellsLocalIn)) then
        nCellsLocal = nCellsLocalIn
        allocate(nCount(1))
    else
        nCellsLocal = nBlocksLocal
    end if

    if((.not. present(nBlocksLocalIn)) .and. (.not. present(nCellsLocalIn)))&
        allocate(nCount(DatasetRank))        

    !Set the nDatasetDimensionions of the DatasetID
    iOneOrZero = 0
    
    !Determine the hdf5 datatype
    if (present(Rank1RealData) .or. present(Rank2RealData) .or.&
        present(Rank3RealData) .or. present(Rank4RealData)) then
        DataType = H5T_NATIVE_DOUBLE
    else if (present(Rank1IntegerData) .or. present(Rank2IntegerData) .or.&
        present(Rank3IntegerData) .or. present(Rank4IntegerData)) then
        DataType = H5T_NATIVE_INTEGER
    else if (present(CharacterData)) then
        call h5tcopy_f(H5T_NATIVE_CHARACTER, DataType, iErrorHdf)
        call h5tset_size_f(DataType, nStringChars, iErrorHdf)
    end if

    call h5screate_simple_f(DatasetRank, nDatasetDimension, DataSpaceId, iErrorHdf) 
    !create the DatasetID
    call h5dcreate_f(FileID, DatasetName, DataType, DataSpaceId, DatasetID, iErrorHdf)

    call h5pcreate_f(H5P_DATASET_XFER_F, PropertyListID, iErrorHdf)

    if (iErrorHdf == -1) &
         call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 1")

    if (DoCollectiveWrite) &
       call h5pset_dxpl_mpio_f(PropertyListID, H5FD_MPIO_COLLECTIVE_F, iErrorHdf)
    
    if (Present(RealAttribute1))&
        call write_hdf5_attribute(NameRealAttribute1, DatasetID, &
            RealAtt = RealAttribute1)

    if (Present(RealAttribute2))&
        call write_hdf5_attribute(NameRealAttribute2, DatasetID, &
            RealAtt = RealAttribute2)

    if (Present(IntegerAttribute1)) &
        call write_hdf5_attribute(NameIntegerAttribute1, DatasetID, &
            IntAtt = IntegerAttribute1)

    if (Present(IntegerAttribute2)) &
        call write_hdf5_attribute(NameIntegerAttribute2, DatasetID, &
            IntAtt = IntegerAttribute2)

    if (nCellsLocal == 0) then
       if (DoCollectiveWrite) then
          call h5screate_simple_f(1, (/iOneOrZero/), MemorySpaceId, iErrorHdf)
          call h5sselect_none_f(MemorySpaceId,iErrorHdf)
          call h5sselect_none_f(DataSpaceId,iErrorHdf)
       else
        iOneOrZero = 1
        call h5screate_simple_f(1, (/iOneOrZero/), MemorySpaceId, iErrorHdf)
        call h5sclose_f(MemorySpaceId,iErrorHdf)
        call h5pclose_f(PropertyListID, iErrorHdf)
        call h5sclose_f(DataSpaceId, iErrorHdf)
        call h5dclose_f(DatasetID, iErrorHdf)
        deallocate(nCount)
        return
       end if
    else
       if(present(nOffsetLocal)) then 
           if(DatasetRank == 1) then 
                nStart(1) = nOffsetLocal
                nCount(1) = nBlocksLocal
            else
               nStart(1:DatasetRank-1) = 0
               nStart(DatasetRank) = nOffsetLocal 
               nCount(1:DatasetRank-1) = nDatasetDimension(1:DatasetRank-1)
               nCount(DatasetRank) = nBlocksLocal
            end if
           !create the hyperslab.  This will differ on the different processors
           call h5sselect_hyperslab_f(DataSpaceId, H5S_SELECT_SET_F, nStart, nCount, iErrorHdf)
           !create the memory space
           call h5screate_simple_f(DatasetRank, nCount, MemorySpaceId, iErrorHdf)
        else if(present(CoordArray)) then
               nCount(1) = nCellsLocal
               !select the cells in the coordinate list
               call h5sselect_elements_f(DataspaceID,H5S_SELECT_SET_F, DatasetRank, nCellsLocal,&
                    CoordArray, iErrorHdf)
               !         !create the memory space
               call h5screate_simple_f(1, (/nCellsLocal/), MemorySpaceID, iErrorHdf)
        else
           nCount(1:DatasetRank) = nDatasetDimension(1:DatasetRank)
           call h5sselect_all_f(DataSpaceId,iErrorHdf)
           !         !create the memory space
           call h5screate_simple_f(DatasetRank, nDatasetDimension, MemorySpaceId, iErrorHdf)
        end if
    end if
   if(DataType == H5T_NATIVE_DOUBLE) then
       if (present(Rank1RealData)) then
           call h5dwrite_f(DatasetID, DataType, Rank1RealData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 6")
       else if (present(Rank2RealData)) then
           call h5dwrite_f(DatasetID, DataType, Rank2RealData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 7")
       else if (present(Rank3RealData)) then
           call h5dwrite_f(DatasetID, DataType, Rank3RealData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 8")
       else if (present(Rank4RealData)) then
           call h5dwrite_f(DatasetID, DataType, Rank4RealData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 9")
        end if
    else if(DataType == H5T_NATIVE_INTEGER) then
       if (present(Rank1IntegerData)) then
           call h5dwrite_f(DatasetID, DataType, Rank1IntegerData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 10")
       else if (present(Rank2IntegerData)) then
           call h5dwrite_f(DatasetID, DataType, Rank2IntegerData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 11")
       else if (present(Rank3IntegerData)) then
           call h5dwrite_f(DatasetID, DataType, Rank3IntegerData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 12")
       else if (present(Rank4IntegerData)) then
           call h5dwrite_f(DatasetID, DataType, Rank4IntegerData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 13")
        end if
    else 
           call h5dwrite_f(DatasetID, DataType, CharacterData, nCount, iErrorHdf, &
                mem_space_id = MemorySpaceId, file_space_id = DataSpaceId,&
                xfer_prp = PropertyListID)
                if (iErrorHdf == -1) &
                    call CON_stop("iErrorHdf in subroutine write_hdf5_data. Error marker 14")
            call h5tclose_f(DataType,iErrorHdf)
    end if
    call h5sclose_f(MemorySpaceId,iErrorHdf)
    call h5pclose_f(PropertyListID, iErrorHdf)
    call h5sclose_f(DataSpaceId, iErrorHdf)
    call h5dclose_f(DatasetID, iErrorHdf)
    deallocate(nCount)
  end subroutine write_hdf5_data

  !=====================================================================
  !=====================================================================

  subroutine write_hdf5_attribute(AttributeName, iDatasetID, RealAtt,&
       IntAtt)
    integer(HID_T), intent(in) :: iDatasetID
    character (len=*), intent(in) :: AttributeName
    real, optional, intent(in) :: RealAtt
    integer, optional, intent(in) :: IntAtt
    integer(HID_T) :: iAttributeSpaceID, attribute, datatype
    integer(HSIZE_T) :: iDimension1D(1)
    integer :: iErrorHdf
    
    iDimension1D=(1)
    call h5screate_simple_f(1, iDimension1D, iAttributeSpaceID, iErrorHdf)

    if (Present(RealAtt)) then
        call h5acreate_f(iDatasetID, AttributeName, H5T_NATIVE_DOUBLE,&
             iAttributeSpaceID, attribute,&
             iErrorHdf, H5P_DEFAULT_F)
       call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, RealAtt, iDimension1D, iErrorHdf)
    elseif (present(IntAtt)) then
        call h5acreate_f(iDatasetID, AttributeName, H5T_NATIVE_INTEGER,&
             iAttributeSpaceID, attribute,&
             iErrorHdf, H5P_DEFAULT_F)
       call h5awrite_f(attribute, H5T_NATIVE_INTEGER, RealAtt, iDimension1D, iErrorHdf)
    endif
    call h5sclose_f(iAttributeSpaceID, iErrorHdf)
    call h5aclose_f(attribute, iErrorHdf)

  end subroutine write_hdf5_attribute
  !=====================================================================
  !=====================================================================
  subroutine close_hdf5_file(FileID) 
    integer(HID_T), intent(in) :: FileID
    integer :: error
    call h5garbage_collect_f(error)
    call h5fclose_f(FileID,error)
    !closing the hdf5 interface
    call h5close_f(error)
    if (error == -1) &
        write (*,*) "close_hdf5_file failed!"
  end subroutine
  !=====================================================================
  !=====================================================================
  
  subroutine pad_string_with_null(nStrings, nCharNeeded, StringsIn, StringsOut)
  integer,intent(in) :: nStrings
  integer, intent(in) :: nCharNeeded
  character (len=*), intent(in) :: StringsIn(nStrings)
  character (len=nCharNeeded), intent(out) :: StringsOut(nStrings)
  integer :: iStr, iLen, nCharsIn
  do iStr = 1, nStrings
       StringsOut(iStr) = StringsIn(iStr)
       nCharsIn = len_trim(StringsOut(iStr))
       do iLen = nCharsIn + 1,nCharNeeded
          StringsOut(iStr)(iLen:iLen) = CHAR(0)
       end do
    end do
  end subroutine pad_string_with_null
end module ModHdf5Utils
