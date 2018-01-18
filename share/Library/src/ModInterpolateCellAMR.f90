!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModInterpolateCellAMR
  use ModInterpolateAMR, ONLY: interpolate_extended_stencil,&
                               iPowerOf2_D, iShift_DI
  implicit none
  integer,parameter:: Coarse_=0, Fine_=1
  !\
  !Interface for tools and methods in ModInterpolateAMR to be used with
  !CELL - adaptive grids, rather than for block-adaptive ones as used in
  !BATSRUS and assumed in ModInterpolateAMR. Example of such grid
  !as used in test_interpolate_cell_amr:
  !/
contains
  !\
  ! Interpolation on cell adaptive Cartesian grid. Similar to block
  ! adaptive grid, but all blocks are 1*1*1 cells
  !/
  subroutine interpolate_amr_cell(nDim, Xyz_D, XyzCell_D, DXyzInv_D, &
       nId, iCellId_IIII, XyzGrid_DIII, iLevel_II, IsOut_II,         &
       nGridOut, Weight_I, iIndexes_II,                              &
       IsSecondOrder)
    !\
    !INPUTS:
    !/
    integer, intent(in)  :: nDim        !Dimensionality
    !\
    !Coordinates of the point at which to interpolate
    !/
    real,    intent(in)  :: Xyz_D(nDim) 
    !\
    ! Coordinates of a center of the cell to which the point Xyz_D belongs
    !/
    real,    intent(in)  :: XyzCell_D(nDim)
    !\
    ! Inverse of the cell size, for each direction
    !/
    real,    intent(in)  :: DXyzInv_D(nDim)
    !\
    ! Length of the cell ID: nID = 1 means that ID is a single integer, 
    ! nID = 2 means that the cell ID is a pair of integers etc
    ! We assume that nID<=nDim
    !/
    integer, intent(in)  :: nID
    !\
    ! The extended stencil in a structured form:
    ! A cubic 2*2*2 grid with 2*2*2 subgrids covering fine vertexes
    ! Here, we have assume the cell to be decomposed to 8(4) octants 
    ! (quadrants) and have 8(4) extended stencils stored for each octants
    ! (quadrants) into which the Xyz point may fall
    !/
    integer, intent(in):: iCellId_IIII(nId,2**nDim,2**nDim,2**nDim)
    real,    intent(in):: XyzGrid_DIII(nDim,0:2**nDim,2**nDim,2**nDim)
    integer, intent(in):: iLevel_II(2**nDim,2**nDim)
    logical, intent(in):: IsOut_II( 2**nDim,2**nDim) 
    !\
    !OUTPUTS:
    !/
    !\
    !Number of grid points involved into interpolation
    !/
    integer, intent(out):: nGridOut      
    !\
    ! Interpolation weights (only the first nGridOut values are meaningful)
    !/
    real,    intent(out):: Weight_I(2**nDim) 
    !\
    ! Cell IDs and processor number for grid points to be
    ! invilved into interpolation. iProc numbers are stored in 
    ! iIndexes_II(0,:) 
    !/
    integer, intent(out):: iIndexes_II(nID,2**nDim)
    !\
    !The following is true if stencil does not employ 
    !the out-of-grid points
    !/
    logical, intent(out), optional:: IsSecondOrder  
    !\
    ! To call the interpolation routine we need only one extended stencil
    !/
    integer            :: iCellIndexes_DII(nId,2**nDim,2**nDim)
    real               :: XyzGrid_DII(nDim,0:2**nDim,2**nDim)
    integer            :: iLevel_I(2**nDim)
    logical            :: IsOut_I(2**nDim)    
    real               :: DXyzInvInput_D(nDim)
    integer            :: iIndexesExt_II(0:nId,2**nDim)
    integer, parameter :: iBlock_I(8) = -1, iProc_I(8) = 0!Unused
    !\
    ! To find octant the point Xyz_D belongs to.
    !/
    integer :: iOctant, iDiscr_D(nDim)
    integer :: nGrid
    !-----------------
    nGrid=2**nDim
    !\
    ! Find octant to which the point Xyz_D belongs
    !/  
    iDiscr_D = nint(0.50 + SIGN(0.50,Xyz_D - XyzCell_D))
    iOctant  = sum(iDiscr_D*iPowerOf2_D(1:nDim)) + 1

    XyzGrid_DII = XyzGrid_DIII(:,:,:,iOctant)
    iCellIndexes_DII = -1 
    iCellIndexes_DII(1:nID,:,:) = iCellId_IIII(:,:,:,iOctant)
    iLevel_I = iLevel_II(:,iOctant)
    IsOut_I  =  IsOut_II(:,iOctant)
    if(any(iLevel_I==-1))then
       iLevel_I = iLevel_I +1
       DXyzInvInput_D = 0.50*DXyzInv_D
    else
       DXyzInvInput_D = DXyzInv_D
    end if
    !\
    ! call interpolation routine
    !/ 
    call interpolate_extended_stencil(&
         nDim            = nDim, &
         Xyz_D           = Xyz_D, &
         nIndexes        = nId, &
         XyzGrid_DII     = XyzGrid_DII, &
         iCellIndexes_DII= iCellIndexes_DII, & 
         iBlock_I        = iBlock_I(1:nGrid),   & 
         iProc_I         = iProc_I( 1:nGrid),   &
         iLevelSubgrid_I = iLevel_I,   & 
         IsOut_I         = IsOut_I,    & 
         DxyzInv_D       = DXyzInvInput_D,&
         nGridOut        = nGridOut,   & 
         Weight_I        = Weight_I,   & 
         iIndexes_II     = iIndexesExt_II, & 
         IsSecondOrder   = IsSecondOrder)
    iIndexes_II(:,1:nGridOut) = iIndexesExt_II(1:nId,1:nGridOut)
  end subroutine interpolate_amr_cell
  !=============================
  !\
  ! Convert non-ordered connectivity list for cell-adaptive grid 
  ! to the arrays used by interpolate_amr_cell
  !/
  subroutine order_connectivity_list(nDim, XyzCell_D, DXyz_D, nId, iCellId_I,&
       nNeighbor, XyzNeighbor_DI, DiLevelNei_I, iCellId_II, IsOut_I,        &
       DXyzInv_D, iCellId_IIII, XyzGrid_DIII, iLevel_II, IsOut_II)
    !\
    ! INPUTS:
    !/
    integer,       intent(in):: nDim 
    !\
    !Coordinates of the cell center
    !/   
    real,         intent(in) ::  XyzCell_D(nDim)  
    !\
    !Cell size along each direction
    !/
    real,         intent(in) ::  DXyz_D(nDim) 
    !\
    !Number of neighboring cells in the connectivity list
    !/      
    !\
    ! Length of the cell ID: nID = 1 means that ID is a single integer, 
    ! nID = 2 means that the cell ID is a pair of integers etc
    ! We assume that nID<=nDim
    !/
    integer,      intent(in) :: nID
    integer,      intent(in) :: iCellId_I(nId)
    integer,      intent(in) :: nNeighbor      
    !\
    ! Cell center coordinates for the neighboring cells
    !/ 
    real,         intent(in) :: XyzNeighbor_DI(nDim,nNeighbor)  
    !\
    !Refinement level of neighboring cells with respect to the cell, 
    !which includes Xyz_D point. May be equal to 
    ! -1, if the neighboring cell is coarser                \ than the cell,
    !  0, if the neighboring cell is at the same resolution |-which includes
    !  1, if the neighboring cell is finer                  / Xyz_D point
    !/
    integer,      intent(in) :: DiLevelNei_I(  nNeighbor)
    integer,      intent(in) :: iCellId_II(nId,nNeighbor)
    logical,      intent(in) :: IsOut_I(nNeighbor)  !.true. for ghost cells 
    !\
    !OUPUTS:
    !/     
    !\
    ! Inverse of the cell size, for each direction
    !/
    real,    intent(out)      :: DXyzInv_D(nDim)
    !\
    ! The extended stencil in a structured form:
    ! A cubic 2*2*2 grid with 2*2*2 subgrids covering fine vertexes
    ! Here, we have assume the cell to be decomposed to 8(4) octants 
    ! (quadrants) and have 8(4) extended stencils stored for each octants
    ! (quadrants) into which the Xyz point may fall
    !/
    integer, intent(out):: iCellId_IIII(nId,2**nDim,2**nDim,2**nDim)
    real,    intent(out):: XyzGrid_DIII(nDim,0:2**nDim,2**nDim,2**nDim)
    integer, intent(out):: iLevel_II(2**nDim,2**nDim)
    logical, intent(out):: IsOut_II( 2**nDim,2**nDim) 
    !\
    !Local variables
    !/
    !\
    !Left corner of the central cell
    !/
    real    :: XyzCorner_D(1:nDim)
    !\
    ! Loop variables
    !/
    integer :: iNeighbor, iGrid, iSubGrid
    !\
    ! Number of grid points in the stencil
    !/ 
    integer :: nGrid
    !\
    ! Ordered list of neighbors
    !/
    integer:: iLevel_III(-1:1,-1:1,2-nDim:nDim-2)
    real   :: XyzGrid_DIIII(1:nDim,0:2**nDim,-1:1,-1:1,2-nDim:nDim-2)
    integer:: iCellId_IIIII(nId,2**nDim,-1:1,-1:1,2-nDim:nDim-2)
    logical:: IsOut_III(-1:1,-1:1,2-nDim:nDim-2)
    !\
    ! Coarse extended stencil
    !/
    logical:: UseCoarseExtStencil !.true. if any neighbor is coarse
    integer:: iCoarseNeighbor     !the neighboring cell which is coarser
    integer:: iGridCoarse         !The grid number  in the coarse stencil 
    real   :: XyzCoarseCenter_D(nDim)
    integer:: iLevelCoarse_I(2**nDim)
    real   :: XyzGridCoarse_DII(1:nDim,0:2**nDim,2**nDim)
    integer:: iCellIdCoarse_III(nId,2**nDim,2**nDim)
    logical:: IsOutCoarse_I(2**nDim)


    logical:: DoRefineGrid_III(-1:1,-1:1,2-nDim:nDim-2) 
    logical:: DoRefineCoarseGrid_I(1:2**nDim)
    !\
    ! Discriminators
    !/
    integer  ::  iDiscr_D(3), iDiscrSubGrid_D(nDim)
    character(LEN=*),parameter::NameSub = 'interpolate_amr_cell_2d'
    !\
    ! Test mode
    !/
    logical, parameter:: DoTest = .true.
    !-----------------------
    !\
    !Initialize variables
    !/
    nGrid = 2**nDim
    !\
    !Invert DXyz_D and find XyzCorner_D
    !/
    DXyzInv_D   = 1/DXyz_D 
    XyzCorner_D = XyzCell_D - 0.50*DXyz_D
    DoRefineGrid_III = .true.
    iDiscr_D = 0
    !\
    ! Check if UseCoarseExtStencil is needed, initalize it, if needed
    UseCoarseExtStencil = any(DiLevelNei_I==-1.and..not.IsOut_I)
    iLevelCoarse_I    = 0
    XyzGridCoarse_DII = 0
    iCellIdCoarse_III = -1
    IsOutCoarse_I     = .false.
    if(UseCoarseExtStencil)then
       DoRefineCoarseGrid_I = .true.
       !\
       ! Find a coarse neighbor
       !/
       iCoarseNeighbor = MINLOC(DiLevelNei_I,MASK=.not.IsOut_I,DIM=1)
       iDiscr_D(1:nDim)= nint(0.50 + SIGN(0.50,&
            XyzNeighbor_DI(:,iCoarseNeighbor) - XyzCell_D))
       iGridCoarse = sum(iDiscr_D(1:nDim)*iPowerOf2_D(1:nDim)) +1
       XyzCoarseCenter_D = XyzNeighbor_DI(:,iCoarseNeighbor) +&
            DXyz_D*(1 - 2*iShift_DI(1:nDim,iGridCoarse))
       do iGrid = 1, nGrid
          XyzGridCoarse_DII(:,0,iGrid) =  XyzCoarseCenter_D + &
               DXyz_D*( 2*iShift_DI(1:nDim,iGrid) - 1)
       end do
       !\
       ! Store an info about the coarse cell
       !/
       XyzGridCoarse_DII(:,1,iGridCoarse) = XyzNeighbor_DI(:,iCoarseNeighbor)
       iCellIdCoarse_III(:,1,iGridCoarse) = iCellId_II(:,iCoarseNeighbor)
       iLevelCoarse_I(iGridCoarse) = -1
       !\
       ! Store the information about central cell into the coarse stencil
       !/
       iDiscr_D(1:nDim) = nint(0.50 + &
            SIGN(0.50, XyzCell_D - XyzCoarseCenter_D))
       iGridCoarse = sum(iDiscr_D(1:nDim)*iPowerOf2_D(1:nDim)) +1

       DoRefineCoarseGrid_I(iGridCoarse) = .false.
       !\
       ! Store coordinate information about refined grid
       !/
       iLevelCoarse_I(iGridCoarse) = 0
       do iGrid = 1, nGrid
          XyzGridCoarse_DII(:,iGrid,iGridCoarse) = &
               XyzGridCoarse_DII(:,0,iGridCoarse) + &
               DXyz_D*(-0.50 + iShift_DI(1:nDim,iGrid))
       end do
       
       iDiscrSubGrid_D = nint(0.50 + SIGN(0.50,XyzCell_D &
                  - XyzGridCoarse_DII(:,0,iGridCoarse)))
       iSubGrid = sum(iDiscrSubGrid_D*iPowerOf2_D(1:nDim)) + 1
       iCellIdCoarse_III(:,iSubGrid,iGridCoarse) = &
            iCellId_I
             IsOutCoarse_I(iGridCoarse) = .false.
    end if
    !\
    !Construct an ordered list of neighbors
    !/
    iLevel_III = -5
    XyzGrid_DIIII = 0.0
    iCellId_IIIII = -1
    IsOut_III = .false.
    !\
    ! Store the central cell info
    !/
    iLevel_III(0,0,0) = 0
    XyzGrid_DIIII(:,0,0,0,0) = XyzCell_D
    XyzGrid_DIIII(:,1,0,0,0) = XyzCell_D
    iCellId_IIIII(:,1,0,0,0) = iCellId_I

    do iNeighbor = 1, nNeighbor
       select case(DiLevelNei_I(iNeighbor))
       case(0)
          !\
          ! find the cell location relatively to the central cell
          iDiscr_D(1:nDim) = &
               floor((XyzNeighbor_DI(:,iNeighbor) - XyzCorner_D)* DXyzInv_D)
          !\
          ! Test list of neighbors
          !/
          if(DoTest)then
             if(all(iDiscr_D == 0))call CON_stop(&
                  NameSub//': neighbor is inside the central cell')
             if(any(iDiscr_D > 1.or.iDiscr_D < -1))call CON_stop(&
                  NameSub//': neighbor does not contact the central cell')
          end if
          !\
          !Initialize the grid coordinates, if needed
          !/
          !\
          ! Store the new cell info into the ordered list
          !/
          XyzGrid_DIIII(:,0,iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = &
               XyzNeighbor_DI(:,iNeighbor)
          XyzGrid_DIIII(:,1,iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = &
               XyzNeighbor_DI(:,iNeighbor)  
          iLevel_III(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = 0
          iCellId_IIIII(:,1,iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = &
               iCellId_II(:,iNeighbor)
          IsOut_III(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = IsOut_I(iNeighbor)
          if(UseCoarseExtStencil)then
             !Save the info to coarse extended stencil too
             iDiscr_D(1:nDim) = nint(0.50 + &
                  SIGN(0.50, XyzNeighbor_DI(:,iNeighbor) - XyzCoarseCenter_D))
             iGridCoarse = sum(iDiscr_D(1:nDim)*iPowerOf2_D(1:nDim)) +1

             if(DoRefineCoarseGrid_I(iGridCoarse))then
                DoRefineCoarseGrid_I(iGridCoarse) = .false.
                !\
                ! Store coordinate information about refined grid
                !/
                iLevelCoarse_I(iGridCoarse) = 0
                do iGrid = 1, nGrid
                   XyzGridCoarse_DII(:,iGrid,iGridCoarse) = &
                        XyzGridCoarse_DII(:,0,iGridCoarse) + &
                        DXyz_D*(-0.50 + iShift_DI(1:nDim,iGrid))
                end do
             end if
             iDiscrSubGrid_D = nint(0.50 + SIGN(0.50,XyzNeighbor_DI(:,iNeighbor)&
                  - XyzGridCoarse_DII(:,0,iGridCoarse)))
             iSubGrid = sum(iDiscrSubGrid_D*iPowerOf2_D(1:nDim)) + 1
             iCellIdCoarse_III(:,iSubGrid,iGridCoarse) = &
                  iCellId_II(:,iNeighbor)
             IsOutCoarse_I(iGridCoarse) = &
                  IsOutCoarse_I(iGridCoarse).or.IsOut_I(iNeighbor)
          end if
       case(1)
          iDiscr_D(1:nDim) = &
               floor((XyzNeighbor_DI(:,iNeighbor) - XyzCorner_D)* DXyzInv_D)
          !\
          ! Test list of neighbors
          !/
          if(DoTest)then
             if(all(iDiscr_D == 0))call CON_stop(&
                  NameSub//': neighbor is inside the central cell')
             if(any(iDiscr_D > 1.or.iDiscr_D < -1))call CON_stop(&
                  NameSub//': neighbor does not contact the central cell')
          end if
          if(DoRefineGrid_III(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)))then
             DoRefineGrid_III(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = .false.
             !\
             ! Store coordinate information about refined grid
             !/
             XyzGrid_DIIII(:,0,iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = &
                  XyzGrid_DIIII(:,0,0,0,0) + DXyz_D*iDiscr_D(1:nDim) 
             iLevel_III(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = 1
             do iGrid = 1, nGrid
                XyzGrid_DIIII(:,iGrid,iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = &
                     XyzGrid_DIIII(:,0,iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) + &
                     DXyz_D*0.50*(-0.50 + iShift_DI(1:nDim,iGrid))
             end do
          end if
          !\
          ! Put the neighbor info into the ordered list
          !/
          !\
          ! If any of finer neighbors is behind the boundary
          ! set IsOut to .true.
          !/
          IsOut_III(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = &
               IsOut_III(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)).or.&
               IsOut_I(iNeighbor)
          iDiscrSubGrid_D = nint(0.50 + SIGN(0.50,XyzNeighbor_DI(:,iNeighbor) &
               - XyzGrid_DIIII(:,0,iDiscr_D(1),iDiscr_D(2),iDiscr_D(3))))
          iSubGrid = sum(iDiscrSubGrid_D*iPowerOf2_D(1:nDim)) + 1
          iCellId_IIIII(:,iSubGrid,iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = &
               iCellId_II(:,iNeighbor)
       case(-1)
          !\
          ! Coarse grid point may need to be placed to several
          ! positions in the ordered array
          !/ 
          do iSubGrid = 1, nGrid
             !\
             ! Decompose the coarser cell into finer cells and
             ! calculate discriminators for each subcell center
             !/
             iDiscr_D(1:nDim) = floor( -0.5 + iShift_DI(1:nDim,iSubGrid) + &
                  (XyzNeighbor_DI(:,iNeighbor) - XyzCorner_D)* DXyzInv_D)
             !\
             ! Ignore subcells, which do not contact the central cell
             !/
             if(any(iDiscr_D > 1.or.iDiscr_D < -1))CYCLE
             !\
             ! Test list of neighbors
             !/
             if(DoTest)then
                if(all(iDiscr_D == 0))call CON_stop(&
                     NameSub//': neighbor is inside the central cell')
             end if
             iLevel_III(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3)) = -1
          end do
          if(iNeighbor == iCoarseNeighbor)CYCLE
          !\
          !Store the coarse cell information in coarse extended stencil
          !/
          iDiscr_D(1:nDim) = nint(0.50 + &
               SIGN(0.50, XyzNeighbor_DI(:,iNeighbor) - XyzCoarseCenter_D))
          iGridCoarse = sum(iDiscr_D(1:nDim)*iPowerOf2_D(1:nDim)) +1
          XyzGridCoarse_DII(:,1,iGridCoarse) = XyzNeighbor_DI(:,iNeighbor)
          iCellIdCoarse_III(:,1,iGridCoarse) = iCellId_II(:,iNeighbor)
          iLevelCoarse_I(iGridCoarse) = -1
          IsOutCoarse_I(iGridCoarse) = IsOut_I(iNeighbor)
       end select
    end do
    !\
    !Form output arrays
    !/
    do iGrid = 1, nGrid
       iDiscr_D(1:nDim) = 2*iShift_DI(1:nDIm,iGrid) - 1
       iLevel_II(:,iGrid) = reshape(iLevel_III(&
            (/MIN(0, iDiscr_D(1)), MAX(0, iDiscr_D(1))/), &
            (/MIN(0, iDiscr_D(2)), MAX(0, iDiscr_D(2))/), &
            (/MIN(0, iDiscr_D(3)), MAX(0, iDiscr_D(3))/)  ), (/nGrid/))
       IsOut_II(:,iGrid) =  reshape(IsOut_III(&
            (/MIN(0, iDiscr_D(1)), MAX(0, iDiscr_D(1))/), &
            (/MIN(0, iDiscr_D(2)), MAX(0, iDiscr_D(2))/), &
            (/MIN(0, iDiscr_D(3)), MAX(0, iDiscr_D(3))/)  ), (/nGrid/))
       if(all(iLevel_II(:,iGrid)>=0.or.IsOut_II(:,iGrid)))then

          iCellId_IIII(:,:,:,iGrid) = reshape(iCellId_IIIII(:,:,&
               (/MIN(0, iDiscr_D(1)), MAX(0, iDiscr_D(1))/), &
               (/MIN(0, iDiscr_D(2)), MAX(0, iDiscr_D(2))/), &
               (/MIN(0, iDiscr_D(3)), MAX(0, iDiscr_D(3))/)  ),&
               (/nId,nGrid,nGrid/))
          XyzGrid_DIII(:,:,:,iGrid) = reshape(XyzGrid_DIIII(:,:,&
               (/MIN(0, iDiscr_D(1)), MAX(0, iDiscr_D(1))/), &
               (/MIN(0, iDiscr_D(2)), MAX(0, iDiscr_D(2))/), &
               (/MIN(0, iDiscr_D(3)), MAX(0, iDiscr_D(3))/)  ),&
               (/nDim,nGrid+1,nGrid/))
       else
          iLevel_II(:,iGrid) = iLevelCoarse_I
          IsOut_II(:,iGrid) =  IsOutCoarse_I

          iCellId_IIII(:,:,:,iGrid) = iCellIdCoarse_III
          XyzGrid_DIII(:,:,:,iGrid) = XyzGridCoarse_DII
       end if
    end do
  end subroutine order_connectivity_list
  !
end module ModInterpolateCellAMR
