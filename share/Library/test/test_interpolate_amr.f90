!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS============================================
module ModTestInterpolateAMR
  use ModInterpolateAMR, ONLY: interpolate_amr
  implicit none
  !\
  ! Shift of the iGrid point in the stencil with respect to the
  ! first one
  !/
  integer, dimension(3,8),parameter :: iShift_DI = reshape((/&
       0, 0, 0, &
       1, 0, 0, &
       0, 1, 0, &
       1, 1, 0, &
       0, 0, 1, &
       1, 0, 1, &
       0, 1, 1, &
       1, 1, 1/),(/3,8/))

  !\
  ! For test: arrays of the refinement levels to be passed to
  ! find_test routine
  !/
  integer:: iLevelTest_I(8)
  integer,parameter:: nCell = 2
contains
  !==================================================================
  subroutine test_interpolate_amr(nDim,nSample)
    use ModRandomNumber, ONLY: random_real

    integer, intent(in)::nDim, nSample

    integer :: iIndexes_II(0:nDim+1,2**nDim)
    logical :: IsSecondOrder
    real, dimension(nDim):: DxyzDomain_D, DxyzCoarseBlock_D, &
         DxyzFineBlock_D, DxyzCoarse_D, &
         DxyzFine_D, Xyz_D,             &
         XyzCont_D,                     &
         XyzInterpolated_D, XyzCorner_D
    real    ::VarInterpolated, VarContInterpolated
    real, allocatable::Xyz_DCB(:,:,:,:,:)    
    real, allocatable, dimension(:,:,:,:) :: Var_CB
    real    :: Weight_I(2**nDim)
    !Loop variables
    integer :: iCase, iSample, iGrid, iSubGrid, i, j, k, iBlock, iDir

    integer :: nCell_D(3)  ! Cells per block
    integer :: iCellIndex_D(3)
    integer :: nIndexes
    integer:: iMisc , nGridOut

    integer:: iSeed = 1
    !-------------------------------------------------------------------
    nCell_D = 1; nCell_D(1:nDim) = nCell
    DxyzDomain_D      = 2*nCell
    DxyzCoarseBlock_D = nCell
    DxyzFineBlock_D   = 0.5*nCell
    DxyzCoarse_D      = 1
    DxyzFine_D        = 0.5
    allocate(Xyz_DCB(nDim,nCell_D(1), nCell_D(2), nCell_D(3),&
         (2**nDim)*(2**nDim+1)))
    Xyz_DCB = 0
    allocate(Var_CB(nCell_D(1), nCell_D(2), nCell_D(3),&
         (2**nDim)*(2**nDim+1)))
    Var_CB = 0
    do iGrid = 1, 2**nDim
       iBlock = iGrid
       XyzCorner_D = DxyzCoarseBlock_D*iShift_DI(1:nDim,iGrid)
       do k = 1, nCell_D(3)
          do j = 1, nCell_D(2)
             do i = 1, nCell_D(1)
                iCellIndex_D = (/i,j,k/)
                Xyz_DCB(:,i,j,k,iBlock) = XyzCorner_D +&
                     DxyzCoarse_D*(iCellIndex_D(1:nDim) - 0.50)
                Var_CB(i,j,k,iBlock) = random_real(iSeed)
             end do
          end do
       end do
       do iSubGrid = 1, 2**nDim
          iBlock = iGrid*(2**nDim)+iSubGrid
          XyzCorner_D = DxyzCoarseBlock_D*iShift_DI(1:nDim,iGrid) + &
               DxyzFineBlock_D*iShift_DI(1:nDim,iSubGrid)
          do k = 1, nCell_D(3)
             do j = 1, nCell_D(2)
                do i = 1, nCell_D(1)
                   iCellIndex_D = (/i,j,k/)
                   Xyz_DCB(:,i,j,k,iBlock) = XyzCorner_D +&
                        DxyzFine_D*(iCellIndex_D(1:nDim) - 0.50)
                   Var_CB(i,j,k,iBlock) = 0.25 + 0.50 * random_real(iSeed)
                end do
             end do
          end do
       end do
    end do

    nIndexes = nDim +1
    CASE:do iCase = 0, 2**(2**nDim) - 2
       iLevelTest_I = 0; iGrid = 0
       iMisc = iCase
       do while(iMisc > 0)
          iGrid = iGrid + 1
          iLevelTest_I(iGrid) = mod(iMisc, 2)
          iMisc = (iMisc - iLevelTest_I(iGrid))/2
       end do
       !write(*,*)'Case=',iLevelTest_I(1:2**nDim)
       !\
       ! We generated refinement, now sample points
       !/
       SAMPLE:do iSample = 1, nSample
          do iDir = 1, nDim
             Xyz_D(iDir) = (0.01 +0.98*random_real(iSeed))*DxyzDomain_D(iDir)
          end do
          !\
          ! call interpolate_amr
          !/
          call interpolate_amr(&
               nDim=nDim, &
               XyzIn_D=Xyz_D, &
               nIndexes=nDim+1,&
               find=find_test, &
               nCell_D=nCell_D(1:nDim),&
               nGridOut=nGridOut,&
               Weight_I=Weight_I,&
               iIndexes_II=iIndexes_II,&
               IsSecondOrder=IsSecondOrder)
          !\          
          !Compare with interpolated:
          !/
          XyzInterpolated_D = 0
          VarInterpolated   = 0
          do iGrid = 1, nGridOut
             iCellIndex_D = 1
             iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
             iBlock = iIndexes_II(nIndexes,iGrid)
             XyzInterpolated_D = XyzInterpolated_D + Weight_I(iGrid)*&
                  Xyz_DCB(:,iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
             VarInterpolated = VarInterpolated + &
                  Weight_I(iGrid)*&
                  Var_CB(iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
          end do
          if(any(abs(Xyz_D - XyzInterpolated_D) > 1.0e-6).and.&
               IsSecondOrder)then
             write(*,*)'Approximation test failed'
             write(*,*)'Grid:', iLevelTest_I(1:2**nDim)
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'Point=', Xyz_D
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DCB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'Xyz_D=',Xyz_D
             write(*,*)'XyzInterpolated_D=',XyzInterpolated_D
             call CON_stop('Correct code and redo test')
          end if
          !\
          ! Test continuity
          !/
          do iDir =1, nDim
             XyzCont_D(iDir) = Xyz_D(iDir) + (0.02*random_real(iSeed) - 0.01)
          end do
          !\
          ! call interpolate_amr
          !/
          call interpolate_amr(&
               nDim=nDim, &
               XyzIn_D=XyzCont_D, &
               nIndexes=nDim+1,&
               find=find_test, &
               nCell_D=nCell_D(1:nDim),&
               nGridOut=nGridOut,&
               Weight_I=Weight_I,&
               iIndexes_II=iIndexes_II,&
               IsSecondOrder=IsSecondOrder)
          !\          
          !Compare interpolated values of Var:
          !/
          VarContInterpolated = 0
          XyzInterpolated_D = 0
          do iGrid = 1, nGridOut
             iCellIndex_D = 1
             iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
             iBlock = iIndexes_II(nIndexes,iGrid)
             VarContInterpolated = VarContInterpolated + &
                  Weight_I(iGrid)*&
                  Var_CB(iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
             XyzInterpolated_D = XyzInterpolated_D + Weight_I(iGrid)*&
                  Xyz_DCB(:,iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
          end do
          if(any(abs(XyzCont_D - XyzInterpolated_D) > 1.0e-6).and.&
               IsSecondOrder)then
             write(*,*)'Approximation test failed'
             write(*,*)'Grid:', iLevelTest_I(1:2**nDim)
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'PointCont=', XyzCont_D
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DCB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'XyzCont_D=',XyzCont_D
             write(*,*)'XyzInterpolated_D=',XyzInterpolated_D
             call CON_stop('Correct code and redo test')
          end if
          if(abs(VarContInterpolated - VarInterpolated) > nDim*0.01)then
             write(*,*)'Continuity test failed'
             write(*,*)'Grid:', iLevelTest_I
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'XyzCont=', XyzCont_D
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DCB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'Xyz_D=',Xyz_D
             call interpolate_amr(&
                  nDim=nDim, &
                  XyzIn_D=Xyz_D, &
                  nIndexes=nDim+1,&
                  find=find_test, &
                  nCell_D=nCell_D(1:nDim),&
                  nGridOut=nGridOut,&
                  Weight_I=Weight_I,&
                  iIndexes_II=iIndexes_II,&
                  IsSecondOrder=IsSecondOrder)
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DCB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             call CON_stop('Correct code and redo test')
          end if
       end do SAMPLE
    end do CASE
    deallocate(Xyz_DCB, Var_CB)

  end subroutine test_interpolate_amr
  !============================
  subroutine find_test(nDim, Xyz_D, &
            iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
    integer, intent(in) :: nDim
    !\
    ! "In"- the coordinates of the point, "out" the coordinates of the
    ! point with respect to the block corner. In the most cases
    ! XyzOut_D = XyzIn_D - XyzCorner_D, the important distinction,
    ! however, is the periodic boundary, near which the jump in the
    ! stencil coordinates might occur. To handle the latter problem,
    ! we added the "out" intent. The coordinates for the stencil
    ! and input point are calculated and recalculated below with
    ! respect to the block corner. 
    !/
    real,  intent(inout):: Xyz_D(nDim)
    integer, intent(out):: iProc, iBlock !processor and block number
    !\
    ! Block left corner coordinates and the grid size:
    !/
    real,    intent(out):: XyzCorner_D(nDim), Dxyz_D(nDim)
    logical, intent(out):: IsOut !Point is out of the domain.
    real, dimension(nDim):: DxyzDomain_D, DxyzCoarseBlock_D ,&
         DxyzFineBlock_D, DxyzCoarse_D, DxyzFine_D
    integer:: iShift_D(3), iGrid, iSubGrid
    integer, dimension(0:1,0:1,0:1), parameter:: iGridFromShift_III=reshape(&
         (/1, 2, 3, 4, 5, 6, 7, 8/),(/2, 2, 2/))
    logical, dimension(nDim) :: IsAboveCenter_D
    !------------------- 
    DxyzDomain_D      = 2*nCell
    DxyzCoarseBlock_D = nCell
    DxyzFineBlock_D   = 0.5*nCell 
    DxyzCoarse_D      = 1
    DxyzFine_D        = 0.5
    iProc = 0; iBlock=0; XyzCorner_D=0.0; Dxyz_D = 0.0
    IsOut = any(Xyz_D < 0.0 .or. Xyz_D >= DxyzDomain_D)
    if(IsOut) RETURN
    !\
    ! Find into which coarse block the point fall
    !/ 
    IsAboveCenter_D = Xyz_D >= DxyzCoarseBlock_D
    iShift_D = 0
    where(IsAboveCenter_D)iShift_D(1:nDim) = 1
    XyzCorner_D = XyzCorner_D + DxyzCoarseBlock_D*iShift_D(1:nDim)
    Xyz_D       = Xyz_D       - DxyzCoarseBlock_D*iShift_D(1:nDim)
    iGrid = iGridFromShift_III(iShift_D(1),iShift_D(2),iShift_D(3))
    !\
    ! Check if the coarse block is used
    !/
    if(iLevelTest_I(iGrid)==0)then
       iBlock = iGrid
       Dxyz_D = DxyzCoarse_D
       RETURN
    end if
    !\
    ! The coarser block is refined, find into which fine block 
    ! the point falls
    !/
    IsAboveCenter_D = Xyz_D >= DxyzFineBlock_D
    iShift_D = 0
    where(IsAboveCenter_D)iShift_D(1:nDim) = 1
    XyzCorner_D = XyzCorner_D + DxyzFineBlock_D*iShift_D(1:nDim)
    Xyz_D       = Xyz_D       - DxyzFineBlock_D*iShift_D(1:nDim)
    iSubGrid = iGridFromShift_III(iShift_D(1),iShift_D(2),iShift_D(3))
    iBlock = iGrid*(2**nDim)+iSubGrid
    Dxyz_D = DxyzFine_D
  end subroutine find_test
end module ModTestInterpolateAMR

program test_interpolate_amr

  use ModTestInterpolateAMR, test => test_interpolate_amr

  implicit none

  call test(2,20000)
  call test(3,20000)

end program test_interpolate_amr

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  !----------------------------------------------------------------------------

  write(*,'(a)')StringError
  write(*,'(a)')'!!! SWMF_ABORT !!!'
  stop

end subroutine CON_stop

