!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!------------------------------------------------------------------------------
! $Id: ModRCMR.f90,v 1.5 2016/08/01 12:42:00 ridley Exp $
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
  use ModTime
  implicit none

  logical :: RCMRFlag = .false.

  integer :: row, col, max_rows, max_cols, print_i, print_j, N, M, mm, dbuffer
  integer :: C_on, ustep, f_o_count, Pc, Pcc, lu2, l_dim, AllocateStatus, idty
  integer :: DeAllocateStatus, control_on, TRUTH_or_ID, lat_p, lon_p, alt_p
  integer :: Nc, lz, lu, ly, s_dhat, ii, kkkk

  double precision :: eta, reg_val, Dts, Measure_Dts, scatter

  character (len=50) :: filename
  character (len=iCharLen_) :: RCMRInType, RCMROutType, RCMRrun
	
  integer, dimension(1,1) :: dhat

  double precision, dimension(1,1) :: u_out, usum, UB, lambda, inp, y_k

  double precision, dimension(:), allocatable :: gathered, scattered
  double precision, dimension(:), allocatable :: gathered_sza, Sat_Proc
  double precision, dimension(:,:), allocatable :: P1, R2, T, theta1
  double precision, dimension(:,:), allocatable :: y_out, w, u, up, y0s, y, y0
  double precision, dimension(:,:), allocatable :: z, zp, zav
  double precision, dimension(:,:), allocatable :: diagn, y_mat, u_mat, z_mat


 
  !ANKIT: Variables below this line are for EDC RCMR
  double precision, dimension(:,:), allocatable :: TEC_true, TEC_lon, TEC_lat
  double precision, dimension(:,:), allocatable :: TEC_currentTime
  integer, dimension(:,:), allocatable :: TEC_step


  integer :: TEC_read_IOStatus = 0
  integer :: kkk = 1                   !counter for reading TEC data
  integer :: TimeArrayDummy = 0        !dummy variable to store timearray entries
!  integer :: iError = 0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! RCAC One step variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  INTEGER :: lu, ly, lz, lv, RegZ
  INTEGER :: lv, RegZ
  INTEGER :: ltheta, lphi
!  INTEGER :: Nc
  integer :: dummy_int  ! to store tuning setting
  !DOUBLE precision :: dummy_real
  DOUBLE precision, dimension(20) :: dummy_real



  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Rz, Ru, Ruf, Rtheta
  DOUBLE PRECISION :: W_Rz, W_Ru, W_Ruf, W_Rtheta

  ! Filter settings
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Nu, Du, Nz, Dz, Nphi, Dphi
  INTEGER :: nf_Nu, nf_Du, nf_Nz, nf_Dz, nf_Nphi, nf_Dphi
  INTEGER :: nf_z, nf_u
  ! Xbar, Xfbar are Buffers for filters
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Ubar, Ufbar, Zbar, Zfbar
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHIbar, PHIfbar
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_filtered, z_filtered, phi_filtered
  
  ! u_h, z_h and y_h are needed to construct phi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_h, y_h, z_h
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_reg
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_reg_row
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Vphi, Uphi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: VphiVec, UphiVec
  
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: theta_h
  
  DOUBLE PRECISION :: lambda1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PP, Tau, Ginv
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: XX, Zp1, Rp, Rbarinv
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Iltheta
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PpXxpTauInv, XxTheta, Rzp, PRdelta, dTheta
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: thetaout, uout         ! Theta(k) and u(k)
  
  ! Lapack inversion
  INTEGER :: INFO = 0 
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV, WORK
  double precision :: alpha = 1.0          !ANKIT: multiplies z
  integer :: uFiltLength = 90              !Ankit:25Jan15 rcmr filter length

  ! Variables to read and write TEC data - ANKIT 1 January 2015
  integer :: points = 0
  real, allocatable , dimension(:,:):: TEC_location
  character (len=16), dimension(12) :: TEC_file
  integer :: TEC_proc = 0                  !Ankit19Feb2015: proc which outputs TEC data
  character(LEN = 15) :: filename1         !Ankit24Feb2015: name of tec output file
  character(LEN = 20) :: filenameRCMR1     !Ankit24Feb2015: name of tec output file
  !integer :: lz_tec = 1                    !ANkit25Feb2015: number of tec stations to be used in RCMR
  INTEGER, DIMENSION(:), ALLOCATABLE :: index_lz_tec !Ankit25Feb2015: Index of tec stations to be used in RCMR
  double precision :: Gf_gain = 1

contains 

subroutine alloc_rcmr

  ! --------------------------------------------------------------
  ! Set the dimensional RCMR parameters
  ! --------------------------------------------------------------
  
  max_cols    = 1000000   !20 days * 24 hours * 60 min * 30 two seconds = 864000
  max_cols    = 2000000   !20 days * 24 hours * 60 min * 30 two seconds = 864000
  max_rows    = max_cols
  l_dim       = 600
  
  lz = 1
  ly = lz  
  lu = 1
  allocate(index_lz_tec(lz))  !index_lz_tec picks the correct file for TEC assimilation
  do ii = 1, lz
     index_lz_tec(ii) = 1
  end do
  
  W_Rz  = 1.0
  W_Ru  = 0.0
  W_Ruf = 0.0
  ! W_Rtheta is read from UAM.in
  IF (W_Ru == 0.0) THEN
     lv = lz
  ELSE
     lv = lz + lu
  ENDIF

  nf_Nz   = 3
  nf_Dz   = 3
  nf_Nu   = 3 
  nf_Du   = 3
  nf_Nphi = 3
  nf_Dphi = 3

  !Ankit 14 March 2015: Expanding filter size to implement Markov Parameters
  nf_z = 20
  nf_u = 20

  open(unit = 223, file = 'markov_parameters.in', action = 'Read')
  read(223, *, IOSTAT=TEC_read_IOStatus) nf_u
  
  ALLOCATE(Nz(lz,   lz*(nf_z+1) ))
  ALLOCATE(Nu(lz,   lu*nf_u     ))
  ALLOCATE(Nphi(lz, lu*(nf_u+1) ))

  ALLOCATE(Dz(lz,   lz*nf_z))
  ALLOCATE(Du(lz,   lz*nf_u)) 
  ALLOCATE(Dphi(lz, lz*nf_u))

  Nz    = 0 ! RESHAPE( (/1, 0, 0, 1, 0,0,0,0, 0,0,0,0/), shape(Nz) )
  Nu    = 0 ! RESHAPE( (/dummy_real, dummy_real,0, 0,0,0/), shape(Nu) )
  Nphi  = 0 ! RESHAPE( (/0, 0, dummy_real, dummy_real, 0, 0, 0, 0 /), shape(Nphi) )
  Dz    = 0 !RESHAPE( (/0,0,0/), shape(Dz) )
  Du    = 0 !RESHAPE( (/0,0,0/), shape(Du) )
  Dphi  = 0 !RESHAPE( (/0,0,0/), shape(Dphi) )

  Do ii = 1,lz
     Nz(ii,ii) = 1
  enddo
  
  Do kkkk = 1, nf_u
     read(223, *, IOSTAT=TEC_read_IOStatus) Nu(1, kkkk)
     Nphi(1,kkkk+1) = Nu(1,kkkk)
  end Do
  write(*,*) 'Gf used in RCMR is', Nu
  close(223)


  !! size calculations
  IF (RegZ == 1) THEN
     ltheta  = Nc * lu * (lu + lz)
     lphi    = Nc * (lu + lz)
  ELSE
     ltheta  = Nc * lu * (lu + ly)
     lphi    = Nc * (lu + ly)
  ENDIF
  
  !! Regressor buffer initialization
  ALLOCATE(u_h        (lu, NC))
  ALLOCATE(z_h        (lz, NC+1))
  ALLOCATE(y_h        (ly, NC+1))
  ALLOCATE(theta_h    (ltheta, 2))
  ALLOCATE(phi_reg    (lphi, 1))
  ALLOCATE(Uphi       (lu, NC))
  ALLOCATE(UphiVec    (lu* NC, 1))
  ALLOCATE(PHI        (lu, lu*lphi))
  IF (RegZ == 1) THEN
     ALLOCATE(Vphi       (lz, NC))
     ALLOCATE(VphiVec    (lz* NC, 1))
  ELSE
     ALLOCATE(Vphi       (ly, NC))
     ALLOCATE(VphiVec    (ly* NC, 1))
  ENDIF
  u_h     = 0
  z_h     = 0
  y_h     = 0
  theta_h = 0
  phi_reg = 0
  Uphi    = 0
  Vphi    = 0
  UphiVec = 0
  VphiVec = 0

 !! Buffer values
  ALLOCATE(Ubar        (lu, nf_u))
  ALLOCATE(Ufbar       (lz, nf_u))
  
  ALLOCATE(Zbar        (lz, nf_z+1))
  ALLOCATE(Zfbar       (lz, nf_z))
  
  ALLOCATE(PHIbar      (lu*(nf_u+1), ltheta))
  ALLOCATE(PHIfbar     (lz*nf_u, ltheta))
  
  ALLOCATE(u_filtered  (lz, 1))
  ALLOCATE(z_filtered  (lz, 1))
  ALLOCATE(phi_filtered(lz, ltheta))
  
  Ubar    = 0
  Ufbar   = 0
  Zbar    = 0
  Zfbar   = 0
  PHIbar  = 0
  PHIfbar = 0
  u_filtered = 0
  z_filtered = 0
  phi_filtered = 0
  
  !! Initialize Cost function Weights
  ALLOCATE(Rz        (lz, lz))
  ALLOCATE(Ruf       (lz, lz))
  ALLOCATE(Ru        (lu, lu))
  ALLOCATE(Rtheta    (ltheta, ltheta))

  call identity(Rz, lz, W_Rz)
  call identity(Ru, lu, W_Ru)
  call identity(Ruf, lz, W_Ruf)
  call identity(Rtheta, ltheta, W_Rtheta)

  ALLOCATE(XX     (lv, ltheta))
  ALLOCATE(zp1    (lv, 1))
  ALLOCATE(Rp     (lv, lv))
  ALLOCATE(Rbarinv(lv, lv))
  ALLOCATE(Tau    (lv, lv))
  ALLOCATE(Ginv   (lv, lv))
  ALLOCATE(PpXxpTauInv (ltheta, lv))
  ALLOCATE(XxTheta (lv, 1))
  ALLOCATE(Rzp     (lv, 1))
  ALLOCATE(IPIV    (lv))
  ALLOCATE(WORK    (lv))
  ALLOCATE(PRdelta (ltheta, ltheta))
  ALLOCATE(dTheta  (ltheta, 1))
  ALLOCATE(Iltheta    (ltheta, ltheta))
  call identity(Iltheta, ltheta, 1.0_8)
  XX = 0     
  zp1 = 0
  Rp = 0
  Tau = 1
  Ginv = 0
  PpXxpTauInv = 0
  XxTheta = 0
  Rzp = 0
  
  !! Covariance Initialization
  ALLOCATE(PP     (ltheta, ltheta))
  call identity(PP, ltheta, 1/W_Rtheta)
  
  !! theta(k) and u(k) initialization
  ALLOCATE(thetaout  (ltheta, 1))
  ALLOCATE(uout      (lu,1))
  thetaout    = 0
  uout        = 0  

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

  !Ankit25Feb2015
  allocate(TEC_true(max_cols,points))
  allocate(TEC_step(max_cols,points))
  allocate(TEC_Lon(max_cols,points))
  allocate(TEC_Lat(max_cols,points))
  allocate(TEC_currentTime(max_cols,points))

  
end subroutine alloc_rcmr


! sub routine to generate EYE(size)
subroutine identity(R, L, weight)
  implicit none    
  integer :: ii
  integer, intent(in) :: L
  DOUBLE PRECISION , intent(in)   :: weight
  DOUBLE PRECISION, DIMENSION(L,L) :: R
  R = 0
  DO ii = 1,L
     R(ii,ii) = weight
  END DO
end subroutine identity

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
  
  elseif(RCMRInType == "TEC") then
     if(RCMROutType == "EDC") then
        !!write (*,*) "AGB RCMR WARNING: this is a test matrix"
        T = reshape((/ 0.15 /), shape(T))
     endif
  else
     ! Markov matrix has not been established
     write (*,*) "No Markov matrix for this output type: ", RCMROutType
     RCMRFlag = .false.
  end if
end subroutine init_markov_matrix

subroutine init_rcmr
  !SatAltDat = -1e32    ! AGB: changed from 1  !Ankit24May16: Commented out. Dont know where sataltdat is allocated :(
  Sat_Loc   = 1
  ! --------------------------------------------------------------
  ! Set these RCMR parameters
  ! --------------------------------------------------------------

  TRUTH_or_ID = 1
  !  Dts         = 2.0  !ANKIT: Set from RCMR_tuning.in
  !  Measure_Dts = 60.0

  col         = 1
  eta         = 0.0
  reg_val     = 100.0
  !  C_on Ensures that estimates are not made before the model has settled
  !  C_on        = 120!1460
  dbuffer     = C_on-20 !100!1440
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
  ! Everything above this line was used by Asad's 2 step code 
  
  do idty = 1, Nc*(ly+lu)
     P1(idty,idty) = reg_val
  end do

  R2 = 0.0
  do idty = 1, lz
     R2(idty,idty) = 1.0
  end do

  !ANKIT: Variables below this line are for EDC RCMR
  TEC_true = 0
  TEC_step = 0
  TEC_lon  = 0
  TEC_lat  = 0
  TEC_currentTime = 0

  !ANKIT: This piece opens the following files, and 
  !       reads the data in to the arrays
  if ((iproc == 0 ) .and. (RCMROutType == "EDC") .and. (RCMRrun == 'ID') ) then
     TEC_file(1) = "tec_data_1.dat"
     TEC_file(2) = "tec_data_2.dat"
     TEC_file(3) = "tec_data_3.dat"
     TEC_file(4) = "tec_data_4.dat"
     TEC_file(5) = "tec_data_5.dat"
     TEC_file(6) = "tec_data_6.dat"
     TEC_file(7) = "tec_data_7.dat"
     TEC_file(8) = "tec_data_8.dat"
     TEC_file(9) = "tec_data_9.dat"
     TEC_file(10) = "tec_data_10.dat"
     TEC_file(11) = "tec_data_11.dat"
     TEC_file(12) = "tec_data_12.dat"

     do ii = 1,points
        write(*,*) "Reading TEC data from file ", iproc, ii, TEC_location(ii,2)
        open(unit = 22, file = TEC_file(ii), action = 'Read')
        do
           read(22, * ,IOSTAT=TEC_read_IOStatus) TEC_Step(kkk,ii), TEC_currentTime(kkk,ii), TimeArrayDummy, &
                TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, &
                TEC_lon(kkk,ii), TEC_lat(kkk,ii), TEC_true(kkk,ii)
           if (TEC_read_IOStatus < 0) exit
           kkk = kkk + 1
        enddo
        close(22)
        write(*,*) "TEC data read from ", TEC_file(ii)
        kkk = 1
     end do
  endif
  
end subroutine init_rcmr


!Ankit:09May2015 This routine reads data from tec_data_locations.in
!      Lat and Lon are stored in degrees, this routine converts them
!      to radians while storing them in TEC_Location(:,:)
subroutine read_tec_locations
  !use ModRCMR
  implicit none
  !integer, intent(in) :: RCMRFlag1
  integer :: read_stat = 0
  integer :: iii = 1

  open(unit = 111111, file = "tec_data_locations.in", action = 'read')
  read(111111,*) TEC_proc
  read(111111,*) points
  allocate(TEC_location(points,2))
  do iii = 1,points !while (iii < points+1)
     read(111111, *, IOSTAT = read_stat) TEC_location(iii,1), TEC_location(iii,2)
     TEC_location(iii,1) = TEC_location(iii,1) * 3.141592653589793/180
     TEC_location(iii,2) = TEC_location(iii,2) * 3.141592653589793/180
  end do
  close(111111)
  write(*,*) "> Read locations for TEC calculations"


!  if (RCMRrun=='TRUTH') then
!     iii = 1
!     do while (iii < points+1)
!        if (iii < 10) then
!           write(filenameRCMR1, "(A14,I1,A4)") "tec_data_RCMR_", iii, '.dat'
!        else
!           write(filenameRCMR1, "(A14,I2,A4)") "tec_data_RCMR_", iii, '.dat'
!        end if
!        open(unit = 111111+iii, file = trim(filenameRCMR1), action = 'write')
!        iii = iii +1
!     end do
!  elseif (RCMRrun=='ID') then
!     do while (iii < points+1)
!        if (iii < 10) then
!           write(filename1, "(A9,I1,A4)") "tec_data_", iii, '.dat'
!        else
!           write(filename1, "(A9,I2,A4)") "tec_data_", iii, '.dat'
!        end if
!        open(unit = 111111+iii, file = trim(filename1), action = 'write')
!        iii = iii +1
!     end do
!  end if

end subroutine read_tec_locations

subroutine write_TEC_data
  implicit none
  logical:: exist
  real :: VTEC_interp, sza_test
  integer :: iii, ii_tec
  !Ankit:09May2015: This piece can't be put in a subroutine, as the file
  !      pointer is lost with the subroutine. It can be done better, but 
  !      I dont know how to do it better.
  !      Do not waste time thinking why I couldnt put it in a subroutine
  !Ankit23May16: Put the TEC writing in a subroutine
  if (RCMRrun == 'ID') Then
     do iii=1,points! while (iii < points+1)
        if (iii < 10) then
           write(filenameRCMR1, "(A14,I1,A4)") "tec_data_RCMR_", iii, '.dat'
        else
           write(filenameRCMR1, "(A14,I2,A4)") "tec_data_RCMR_", iii, '.dat'
        end if
 
        inquire(file=filenameRCMR1, exist=exist)
        if (exist) then
           open(unit = 111111+iii, file = trim(filenameRCMR1), status="old", position="append", action = 'write')
        else
           open(unit = 111111+iii, file = trim(filenameRCMR1), status="new", action = 'write')
        end if
        
     end do
  elseif (RCMRrun == 'TRUTH') then
     do iii=1,points! while (iii < points+1)
        if (iii < 10) then
           write(filename1, "(A9,I1,A4)") "tec_data_", iii, '.dat'
        else
           write(filename1, "(A9,I2,A4)") "tec_data_", iii, '.dat'
        end if
        inquire(file=filename1, exist=exist)
        if (exist) then
           open(unit = 111111+iii, file = trim(filename1), status="old", position="append", action = 'write')
        else
           open(unit = 111111+iii, file = trim(filename1), status="new", action = 'write')
        end if
           
     end do
  end if
  
  
   !! Ankit 24Jan2015 - Added TEC writing at preset locations
  if ((iproc == TEC_proc) .AND. .true.) then
     do ii_tec = 1,points !while (ii_tec < points+1)
        call calc_single_vtec_interp(TEC_location(ii_tec,1), TEC_location(ii_tec,2), VTEC_interp)
        call get_sza(TEC_location(ii_tec,1), TEC_location(ii_tec,2), sza_test)
        write(111111+ii_tec , &
             "(I7, 1X, F15.4, 1X, I4, 1X, I2, 1X, I2, 1X, I2, 1X, I2, 1X,I2, 1X, I3, 1X, F9.3, 1X, F9.3, 1X, F12.8, 1X, F12.8)") &
             iStep, CurrentTime, iTimeArray(1), iTimeArray(2), iTimeArray(3), iTimeArray(4), iTimeArray(5), iTimeArray(6), &
             iTimeArray(7), TEC_location(ii_tec,1)*180/3.141592653589793, &
             TEC_location(ii_tec,2)*180/3.141592653589793 , VTEC_interp, sza_test
     end do
  end if
  
  
  !  close TEC files - ANKIT
  do while (ii_tec < points+1)
     close(111111+ii_tec)
     ii_tec = ii_tec+1
  end do
  

end subroutine write_TEC_data

end module ModRCMR
