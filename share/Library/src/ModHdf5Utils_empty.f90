!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHdf5Utils

  implicit none

contains
  subroutine save_hdf5_file(FileName,TypePosition, TypeStatus, StringHeader,&
       nStep, NumberOfBlocksUsed, Time, nDimOut, nParam, nVar,&
       n_D, NameVar, NameUnits, MinimumBlockIjk, XYZMinMax, PlotVarBlk,&
       iComm, CoordMin, CoordMax)

    integer, intent(in) :: nDimOut, nParam, nVar, n_D(3)
    integer, intent(in) :: nStep,NumberOfBlocksUsed
    character (len=*),   intent(in) :: FileName
    character (len=*),   intent(in) :: TypePosition, NameVar(nVar)
    character (len=10),  intent(in) :: TypeStatus
    character (len=500), intent(in) :: StringHeader
    character (len=*),   intent(in) ::  NameUnits
    integer, intent(in) :: MinimumBlockIjk(:,:)
    real, intent(in)  :: Time, PlotVarBlk(:,:,:,:,:), XYZMinMax(:,:,:)

    integer, optional, intent(in):: iComm
    real,    optional, intent(in):: CoordMin(:), CoordMax(:)

  end subroutine save_hdf5_file

  !=====================================================================
end module ModHdf5Utils

