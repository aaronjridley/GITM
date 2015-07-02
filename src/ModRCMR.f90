!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!------------------------------------------------------------------------------
! $Id: ModRCMR.f90,v 1.3 2013/10/24 18:40:31 agburr Exp $
!
! Author: Asad
!
! Comments: Common variables needed across GITM to perform RCMR data
!           assimilation, as well as subroutines to initialize the variables
!
! AGB 3/31/13: Removed unnecessary variables, changed names to be more
!              descriptive, added variables to allow flagging from UAM.in,
!              added initialization routines to allow flagging from UAM.in
! AGB 10/23/13: Adapted to allow driving of photoelectron heating efficiency
!------------------------------------------------------------------------------

module ModRCMR

  use ModGITM, only:iProc, nProcs, Sat_Loc
  use ModSatellites, only: SatAltDat
  use ModInputs, only: iCharLen_
	
  implicit none

  logical :: RCMRFlag = .false.

  integer :: row, col, max_rows, max_cols, print_i, print_j, N, M, mm, dbuffer
  integer :: C_on, ustep, f_o_count, Pc, Pcc, lu2, l_dim, AllocateStatus, idty
  integer :: DeAllocateStatus, control_on, TRUTH_or_ID, lat_p, lon_p, alt_p
  integer :: Nc, lz, lu, ly, s_dhat, ii

  double precision :: eta, reg_val, Dts, Measure_Dts, scatter

  character (len=50) :: filename
  character (len=iCharLen_) :: RCMRInType, RCMROutType
	
  integer, dimension(1,1) :: dhat

  double precision, dimension(1,1) :: u_out, usum, UB, lambda, inp, y_k

  double precision, dimension(:), allocatable :: gathered, scattered
  double precision, dimension(:), allocatable :: gathered_sza, Sat_Proc
  double precision, dimension(:,:), allocatable :: P1, R2, T, theta1
  double precision, dimension(:,:), allocatable :: y_out, w, u, up, y0s, y, y0
  double precision, dimension(:,:), allocatable :: z, zp, zav
  double precision, dimension(:,:), allocatable :: diagn, y_mat, u_mat, z_mat
contains 

subroutine alloc_rcmr

  ! --------------------------------------------------------------
  ! Set the dimensional RCMR parameters
  ! --------------------------------------------------------------

  max_cols    = 500000
  max_rows    = max_cols
  Nc          = 100
  lz          = 1
  ly          = 1
  lu          = 1
  l_dim       = 600

  ! --------------------------------------------------------------
  ! Allocate space for the RCMR variables
  ! --------------------------------------------------------------
  
  allocate(gathered(nProcs), STAT = AllocateStatus)
  allocate(gathered_sza(nProcs), STAT = AllocateStatus)
  allocate(scattered(nProcs), STAT = AllocateStatus)
  allocate(Sat_Proc(nProcs), STAT = AllocateStatus)
  allocate(T(lz,lu), STAT = AllocateStatus)
  allocate(y0s(nProcs,max_cols), STAT = AllocateStatus)
  allocate(w(1,max_rows), STAT = AllocateStatus)
  allocate(u(lu,max_rows), STAT = AllocateStatus)
  allocate(y_out(max_rows,1), STAT = AllocateStatus)
  allocate(y(lu,max_rows), STAT = AllocateStatus)
  allocate(y0(lz,max_rows), STAT = AllocateStatus)
  allocate(z(lz,max_rows), STAT = AllocateStatus)
  allocate(zp(lz,max_rows), STAT = AllocateStatus)
  allocate(up(lz,max_rows), STAT = AllocateStatus)
  allocate(zav(lz,max_rows), STAT = AllocateStatus)
  allocate(diagn(1,max_rows), STAT = AllocateStatus)
  allocate(z_mat(lz,l_dim), STAT = AllocateStatus)
  allocate(u_mat(lu,l_dim), STAT = AllocateStatus)
  allocate(y_mat(ly,l_dim), STAT = AllocateStatus)
  allocate(R2(lz,lz), STAT = AllocateStatus)
  allocate(P1(lu*Nc+ly*(Nc),lu*Nc+ly*(Nc)), STAT = AllocateStatus)
  allocate(theta1(1,lu*Nc+ly*(Nc)), STAT = AllocateStatus)

end subroutine alloc_rcmr

subroutine init_markov_matrix

  ! T is a Markov parameter matrix that must be tuned for each type of
  ! assimilation

  if(RCMRInType == "RHO") then
     ! Markov matrix for terrestrial neutral mass density
     ! Assumes only one data input source

     if(RCMROutType == "F107") then
        T = reshape((/ 0.15 /), shape(T))
     else if(RCMROutType == "PHOTOELECTRON") then
        ! AGB 7/17/13: This is not settled yet
        write (*,*) "AGB RCMR WARNING: this is a test matrix"
        T = reshape((/ 0.15 /), shape(T))
     else
        write (*,*) "No Markov matrix for this output type: ", RCMROutType
        RCMRFlag = .false.
     end if
  else
     ! Markov matrix has not been established
     write (*,*) "No Markov matrix for this output type: ", RCMROutType
     RCMRFlag = .false.
  end if
end subroutine init_markov_matrix

subroutine init_rcmr

  SatAltDat = -1e32    ! AGB: changed from 1
  Sat_Loc   = 1

  ! --------------------------------------------------------------
  ! Set these RCMR parameters
  ! --------------------------------------------------------------

  TRUTH_or_ID = 1
  Dts         = 2.0
  Measure_Dts = 60.0

  col         = 1
  eta         = 0.0
  reg_val     = 100.0
  ! C_on Ensures that estimates are not made before the model has settled
  C_on        = 1460
  dbuffer     = 1440
  lambda(1,1) = 0.9999
  s_dhat      = 1
  dhat(1,1)   = 1
  ustep       = 1

  do mm = 0,nProcs	
     if (iProc == mm) then
        lat_p = 9
        lon_p = 1
        alt_p = 36
     end if
  end do

  u_out      = 0.0
  Sat_Proc   = 0.0
  y_mat      = 0.0
  z_mat      = 0.0
  u_mat      = 0.0
  theta1     = 0.0
  control_on = 0

  u(:,:)       = 0.0
  gathered(:)  = 0.0
  scattered(:) = 0.0
  usum         = 0.0
  P1           = 0.0

  do idty = 1, Nc*(ly+lu)
     P1(idty,idty) = reg_val
  end do

  R2 = 0.0
  do idty = 1, lz
     R2(idty,idty) = 1.0
  end do

end subroutine init_rcmr

end module ModRCMR
