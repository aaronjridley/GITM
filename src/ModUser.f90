!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModUserGITM

  use ModGITM, only: nBlocksmax,nLons,nLats,nAlts

  implicit none

  integer, parameter :: nUserOutputs = 100
  real :: UserData3D(-1:nLons+2,-1:nLats+2,-1:nAlts+2,nUserOutputs,nBlocksMax)=0.0
  real :: UserData2D(-1:nLons+2,-1:nLats+2,         1,nUserOutputs,nBlocksMax)=0.0
  real :: UserData1D(1,1,1:nalts,nUserOutputs)=0.0

  integer :: nVarsUser3d=0
  integer :: nVarsUser1d=0
  integer :: nVarsUser2d=0

end module ModUserGITM
