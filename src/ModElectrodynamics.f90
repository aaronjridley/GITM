!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModElectrodynamics

  use ModSizeGitm

  ! This is the divergence of the Neutral wind current
  real, dimension(-1:nLons+2,-1:nLats+2,1:nAlts) :: DivJu

  ! This is the height integral of the divergence
  real, dimension(-1:nLons+2,-1:nLats+2) :: DivJuAlt

  ! This is the field-line integral of the conductance and divergence
  real, dimension(-1:nLons+2,-1:nLats+2) :: &
       HallFieldLine, PedersenFieldLine, DivJuFieldLine, LengthFieldLine

  ! This is the field-line integral of the conductance and divergence
  real, dimension(-1:nLons+2,-1:nLats+2) :: &
       SigmaPP, SigmaLL, SigmaHH, SigmaCC, SigmaPL, SigmaLP, &
       KDpm, KDlm, Kpm, Klm
  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: ed1, ed2, je1, je2

  ! This is the field aligned integral in magnetic coordinates
  real, dimension(:,:), allocatable :: DivJuAltMC

  ! These are the conductances in magnetic coordinates
  real, dimension(:,:), allocatable :: SigmaHallMC
  real, dimension(:,:), allocatable :: SigmaPedersenMC

  real, dimension(:,:), allocatable :: SigmaPPMC
  real, dimension(:,:), allocatable :: SigmaLLMC
  real, dimension(:,:), allocatable :: SigmaHHMC
  real, dimension(:,:), allocatable :: SigmaCCMC
  real, dimension(:,:), allocatable :: SigmaPLMC
  real, dimension(:,:), allocatable :: SigmaLPMC
  real, dimension(:,:), allocatable :: AverageMC

  real, dimension(:,:), allocatable :: KDpmMC, kpmMC
  real, dimension(:,:), allocatable :: KDlmMC, klmMC
! New parameters
  real, dimension(:,:), allocatable :: SigmaCowlingMC
  real, dimension(:,:), allocatable :: dSigmaCowlingdpMC
  real, dimension(:,:), allocatable :: dSigmaLLdpMC
  real, dimension(:,:), allocatable :: dKDlmdpMC
  real, dimension(:,:), allocatable :: Ed1new
  real, dimension(:,:), allocatable :: Ed2new
! 
  ! These are the magnetic coordinates
  real, dimension(:,:), allocatable :: MagLatMC
  real, dimension(:,:), allocatable :: MagLocTimeMC
  real, dimension(:,:), allocatable :: MagLonMC
  real, dimension(:,:), allocatable :: GeoLonMC
  real, dimension(:,:), allocatable :: GeoLatMC

  real, dimension(:,:), allocatable :: MagBufferMC
  real, dimension(:,:), allocatable :: LengthMC

  real, dimension(:,:), allocatable :: &
       solver_a_mc, solver_b_mc, solver_c_mc, solver_d_mc, solver_e_mc, &
       solver_s_mc, deltalmc, deltapmc, &
       dSigmaLLdlMC, dSigmaLPdlMC, dSigmaPLdpMC, dSigmaPPdpMC, &
       dKDpmdpMC, dKDlmdlMC, DynamoPotentialMC, &
       dKpmdpMC, dKlmdlMC

  real, dimension(:,:), allocatable :: oldpotmc
  
   real, dimension(:), allocatable :: & 
        x,y,rhs,b,d_I,e_I,e1_I,f_I,f1_I

  real, dimension(:,:), allocatable :: &
       SmallMagLocTimeMC, SmallMagLatMC, SmallPotentialMC

  integer :: nMagLats = 140  ! 1 degrees
  integer :: nMagLons = 90  ! 4 degrees
  real :: MagLatRes = 0.5
  real :: MagLonRes = 4.0
  
  !----------------------------------------------------------------------
  ! These are in geographic coordinates : 

  real, dimension(-1:nLons+2,-1:nLats+2, nBlocksMax) :: &
       HallConductance, PedersenConductance

  real, dimension(-1:nLons+2,-1:nLats+2, -1: nAlts+2) :: &
       Sigma_0, Sigma_Pedersen, Sigma_Hall

  real, dimension(-1:nLons+2,-1:nLats+2, -1: nAlts+2, 3) :: &
       UxB, Ju

  real :: SigmaR(-1:nLons+2,-1:nLats+2, -1: nAlts+2, 3, 3)

end module ModElectrodynamics
