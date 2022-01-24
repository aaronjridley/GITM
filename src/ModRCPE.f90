!----------------------------------------------------------------------------
! $Id: RCAC_routines.f95,v 2.1 2018/05/14$
!                       ,v 2.2 2021/10/26$
!
! Author: ANKIT GOEL, University of Michigan, 10/2021
!
! Comments: Contains definitions of the variables used to run RCPE algorithm 
!
!-----------------------------------------------------------------------------


Module ModRCPE
  use ModInputs, only: iCharLen_
    
	IMPLICIT NONE
        INTEGER :: lu, ly, lz, lv
        INTEGER :: kk = 0
	INTEGER :: ltheta, lphi
	INTEGER :: Nc
	DOUBLE PRECISION :: W_Rtheta, lambda
	
        INTEGER :: nf 					! Order of the FIR filter to be optimized

        DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: Y, Yhat, Z
        DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: mu
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Rtheta
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: g_int 						! Integrator. Size lz by 1
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi, phi_kron 				! Regressor phi. Size Nc(lu+lz)+lz by 1. Size lu by lu (Nc(lu+lz)+lz).
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI_window					! PHI buffer. pn*lu by ltheta
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_window 					! Control buffer. lu by pn 
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: z_window 					! Performance buffer. lz by pc
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: N_filt, theta_star 		! Filter and controller gains. lz by nf*lu. ltheta by 1.
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_RLS
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: I_lu

	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Phi_b_rr
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U_b_rr
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI_filt
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  U_filt
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PtimesPHI_filt
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Gamma_P
        double precision, dimension(:), allocatable :: scattered
 

        !BP: Extra variables from old ModRCMR.f90 that are needed 
        logical :: rcmrFlag
        character (len=iCharLen_) :: RCMRInType, RCMROutType, coefficientToEstimate
        
        !ANK 12/04/17 Multiple estimates updated by RCMR
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RCMR_est, RCMR_est0
        real, DIMENSION(:), ALLOCATABLE :: Nf_Array
        DOUBLE PRECISION :: alpha
        integer :: C_on
        DOUBLE PRECISION :: lambda1
        double precision :: Dts, Measure_Dts, scatter
        integer :: rcmrStepToTurnOff = 1000000000
        integer :: regz
        integer :: uFiltLength = 90   !Ankit:25Jan15 rcmr filter length


        real, dimension(:), allocatable :: satLon, satLat, satAlt, satMeanDensity, &
                                     sat2Lon, sat2Lat, sat2Alt, sat2Density
        integer, dimension(:), allocatable :: satIndex, satYear, satMonth, satDay, &
                                        satHour, satMinute, satSecond
        character(len = 8), dimension(:), allocatable :: time
        character(len = 10), dimension(:), allocatable :: date
        
 

contains
   
!----------------------------------------------------------------------------
! $Id: RCAC_routines.f95,v 2.1 2018/05/14$
!                       ,v 2.2 2021/10/26$
!
! Author: ANKIT GOEL, University of Michigan, 10/2021
!
! Comments: Following Routines to run RCAC one step algorithm are defined
!   initialize_RCPE  ::  Initializes all the variables used by RCAC one step 
 !                    algorithm. Variables are defined in RCAC_variables.F95
 !   identity     ::  Creates the identity matrix
 !   weights      ::  Creates the weighting matrices
 !   kronecker_AB ::  Calculates C = A kron B
 !   RCPE   ::  Uses RCAC algorithm to find u(k)
!
!-----------------------------------------------------------------------------


subroutine initialize_RCPE(lu_val, ly_val, nf_val, lambda_val,W_Rtheta_val)
  implicit none

  INTEGER, intent(in) :: lu_val, ly_val, nf_val
  DOUBLE PRECISION, intent(in) :: W_Rtheta_val, lambda_val
  INTEGER :: ii = 0
  lu = lu_val
  ly = ly_val
  lz = ly
  nf = nf_val
  lambda = lambda_val
  W_Rtheta = W_Rtheta_val
  ! print *, lambda, W_Rtheta
  Nc = 0

  lphi = Nc * (lu + lz) + lz ! phi includes the integrator
  ltheta = lu * ( Nc * (lu + lz) + lz ) ! phi includes the integrator
  
  
  ! 16 Jan 2017: Following variables are allocated for MIMO Gf Optimization
  ALLOCATE(Y(ly,1))
  ALLOCATE(Yhat(ly,1))
  ALLOCATE(Z(lz,1))
  ALLOCATE(mu(lu,1))
  
  ALLOCATE(g_int (lz,1) )
  ALLOCATE(phi (Nc*(lu+lz)+lz,1) )
  ALLOCATE(phi_kron (lu, lu*Nc*(lu+lz)+lu*lz) )
  
  ALLOCATE(PHI_window (2*(Nc+1)*lu, ltheta))
  ALLOCATE(u_window (lu, 2*(Nc+1)))
  ALLOCATE(z_window (lz, 2*(Nc+1)))
  
  ALLOCATE(N_filt (lz, nf * lu))
  ALLOCATE(theta_star (ltheta, 1))
  
  ALLOCATE(P_RLS (ltheta, ltheta))
  ALLOCATE(I_lu (lu, lu))
  ALLOCATE(Phi_b_rr (nf*lu, ltheta))
  ALLOCATE(U_b_rr (nf*lu, 1))
  ALLOCATE(PHI_filt (lz, ltheta))
  ALLOCATE(PtimesPHI_filt (ltheta, lz))
  ALLOCATE(Gamma_P (lz, lz))
  
  ALLOCATE(U_filt (lz, 1))
  g_int = 0
  PHI_window = 0
  u_window = 0
  z_window = 0 

  Phi_b_rr = 0
  U_b_rr   = 0
  
  N_filt = 0
  DO ii = 1,lz
     N_filt(ii,nf*lu) = -1
     END DO
     !N_filt = reshape( (/-.1, 0.0 ,0.0,-.1/), shape(N_filt))
     theta_star = 0
     
     ! Initialize the covariance matrix
     CALL identity(P_RLS, ltheta,1.0/W_Rtheta)
     CALL identity(I_lu, lu, 1D0)

end subroutine initialize_RCPE

subroutine RCPE(u_out, u_in, z_in, y_in)
  implicit none
  integer rr
  !integer, intent(in) :: kk
  DOUBLE PRECISION, intent(in), DIMENSION(lu,1) :: u_in
  DOUBLE PRECISION, intent(in), DIMENSION(lz,1) :: z_in
  DOUBLE PRECISION, intent(in), DIMENSION(ly,1) :: y_in
  DOUBLE PRECISION, intent(out), DIMENSION(lu,1) :: u_out
  !DOUBLE PRECISION, intent(out), DIMENSION(ltheta,1) :: theta_out
  DOUBLE PRECISION, DIMENSION(ltheta,1) :: theta_out
  

  ! --------------------------------------------------------------
  ! Regressor
  ! --------------------------------------------------------------
  ! Put data in control and performance buffers
  u_window = cshift(u_window,-1,2)
  z_window = cshift(z_window,-1,2)
  PHI_window = cshift(PHI_window,-lu,1)

  u_window(1:lu,1:1) = u_in !reshape(u_in, shape(u_window(1:lu,1)))
  z_window(1:lz,1:1) = z_in !reshape(z_in, shape(z_window(1:lz,1)))
  
  ! Construct regressor, phi. Size Nc*(lu+lz)+lz
  g_int = g_int + z_in
  !phi(1:Nc*lu,1) = reshape( 1*u_window(:,1:Nc), shape(phi(1:Nc*lu,1)) )
  !phi(Nc*lu+1:Nc*lu+Nc*lz,1)= reshape( 0*1+1*z_window(:,1:Nc), shape(phi(Nc*lu+1:Nc*lu+Nc*lz,1)) )
  phi(Nc*lu+Nc*lz+1:,1)= reshape( 1*g_int, shape(phi(Nc*lu+Nc*lz+1:,1)) )

  ! Construct regressor, PHI. Size lu by lu*(Nc*(lu+lz)+lz)
  CALL KRONECKER_AB(phi_kron, I_lu, lu, lu, transpose(phi),1, Nc*(lu+lz)+lz)
  
  PHI_window(1:lu,:) = phi_kron !reshape(z_in, shape(z_window(1:lz,1)))

  ! --------------------------------------------------------------
  ! Filter
  ! --------------------------------------------------------------
  
  rr = 1
  Phi_b_rr = PHI_window(lu*rr+1:lu*rr+lu*nf,:)
  call vectorize(u_window(:,rr:rr+nf-1 ), lu, nf, U_b_rr)

  PHI_filt = matmul(N_filt, Phi_b_rr)

  write(*,*) matmul(N_filt,U_b_rr)
  U_filt = matmul(N_filt, U_b_rr)
    
  if ( kk> max(Nf+1, Nc+1) )  then
            
     PtimesPhi_filt = matmul(P_RLS, transpose(PHI_filt))
 
     Gamma_P = 1/(lambda + matmul(PHI_filt, PtimesPhi_filt))
  
     P_RLS = P_RLS - matmul(matmul(PtimesPhi_filt,Gamma_P), transpose(PtimesPhi_filt))
     P_RLS = P_RLS/lambda
 
     theta_star = theta_star - matmul(PtimesPhi_filt, matmul(PHI_filt, theta_star)+z_in-U_filt)

     theta_out = theta_star
     u_out = matmul(phi_kron, theta_out)
  else 
        theta_out = 0
        u_out = u_in
  end if
        ! call AppendMatrix2File(transpose(u_in), 1, lu,     'u_window11111111.dat')
        call AppendMatrix2File(transpose(z_in), 1, lz,     'z_window11111111.dat')
        call AppendMatrix2File(transpose(u_out), 1, lu, 'u_RCMR1111111111.dat')
        call AppendMatrix2File(transpose(theta_out), 1, ltheta,     'theta11111111111.dat')
        call writeMatrix2File(phi_kron,lu, lu*Nc*(lu+lz)+lu*lz, 'phi_kron11111111.dat')
        call writeMatrix2File(P_RLS,ltheta, ltheta,       'P_RLS11111111111.dat')
        ! write(*,*) "RCMR output is ", kk, u_out, u_in
        !write(*,*) "Optimized estimator coefficients are ", theta_out
        !write(*,*) "Optimized filter is ", filter_out
        ! write(*,*) kk, u_out, PHI_filt, U_b_rr,theta_star, lambda

end subroutine RCPE




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



! 19 Jan 2017: Following subroutine computes the kronecker product of two matrices 
SUBROUTINE KRONECKER_AB(C, A, Am, An, B, Bm, Bn)
  IMPLICIT NONE

  integer, intent(in) :: Am, An, Bm, Bn
  DOUBLE PRECISION, intent(in), DIMENSION(Am, An) :: A 
  DOUBLE PRECISION, intent(in), DIMENSION(Bm, Bn) :: B 
  DOUBLE PRECISION, intent(out), DIMENSION(Am*Bm,An*Bn) :: C
  integer :: ii, jj 

  DO ii = 1,Am
     DO jj = 1,An
        C(Bm*(ii-1)+1:Bm*ii, Bn*(jj-1)+1:Bn*jj) = A(ii, jj) * B
        END DO
        END DO
END SUBROUTINE KRONECKER_AB




SUBROUTINE vectorize(U, n, m, vecU)
  implicit none

  integer :: ii
  integer, intent(in) :: n, m
  double precision, intent (in), dimension(n,m) :: U
  double precision, intent (out), dimension(n*m,1) :: vecU
  Do ii = 0,m-1
     vecU(ii*n+1:ii*n+n,1) = U(:,ii+1)
     end do
end subroutine vectorize


subroutine writeMatrix(A,n)
  implicit none

  integer :: ii
  integer, intent(in) :: n
  double precision, intent (in), dimension(n,n) :: A
  do ii = 1,n
     write(*,*), A(ii,:)
     end do


end subroutine writeMatrix

subroutine writeMatrix2File(A,n,m, filename)
  implicit none

  integer :: ii
  integer, intent(in) :: n, m
  CHARACTER(20), intent(in)  :: filename
  double precision, intent (in), dimension(n,m) :: A

  !write(*,*), 'writing ' , filename
  open(unit = 8188, file = filename, action = 'write')
  do ii = 1,n
     write(8188,*), A(ii,:)
     end do

     close(8188)


end subroutine writeMatrix2File


subroutine AppendMatrix2File(A,n,m, filename)
  implicit none

  integer :: ii
  integer, intent(in) :: n, m
  CHARACTER(20), intent(in)  :: filename
  double precision, intent (in), dimension(n,m) :: A

  logical :: exist

  inquire(file=filename, exist=exist)
  if (exist) then
     open(unit = 8188, file = filename, status="old", position="append", action = 'write')
     else
        open(unit = 8188, file = filename, status="new", action = 'write')
        end if

        !write(*,*), 'writing ' , filename
        do ii = 1,n
           write(8188,*), A(ii,:)
           end do

           close(8188)


end subroutine AppendMatrix2File

subroutine mainRCPE
  use ModMPI
  use ModGITM
  use ModInputs
  use ModTime
  use EUA_ModMsis00, ONLY: gtd7

  IMPLICIT NONE

  double precision :: diff = -1e32
  double precision :: diff2 = -1e32
  real :: gitmDensity_SS = -1e32
  real :: msisDensity_SS = -1e32
  real :: gitmDensity_AS, msisDensity_AS

  real :: geo_lat, geo_lon, geo_alt, geo_lst
  real, dimension(1:2) :: msis_temp
  real, dimension(1:9) :: msis_dens
  real, dimension(7)  :: ap = 10.0
  logical :: foundLocation = .False.

  double precision :: localVar
  integer :: iError
  integer :: nLines = 0
  integer :: blockCheck = -1

  real :: SSLon, SSLat, ASLon, ASLat
  real :: rLon, rLat, rAlt

  integer :: iLat = -1
  integer :: iLon = -1
  integer :: iAlt = -1
  integer :: iBlock = -1 
  integer :: iiBlock = -1
  integer :: iiProc = -1  

  !call initialize()     !! Initialize RCAC variables
  diff = -1e32
  diff2 = -1e32
  localVar = -1e32
  gitmDensity_SS = -1e32
  msisDensity_SS = -1e32
  gitmDensity_AS = -1e32
  msisDensity_AS = -1e32

  scatter = 0
  if (nLines .eq. 0) then
     call getNumLines(nLines)
     if (iProc .eq. 0) write(*,*) "Found ", nLines, "in syntheticDataFile.txt."
  endif

  if (mod((istep-1)*Dts,Measure_Dts) == 0.0) then
    call readSyntheticData(nLines)
    call printDensity(nLines, diff, diff2, blockCheck)
    localVar = diff
    call MPI_REDUCE(localVar, diff, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, iCommGITM, iError)
    if (diff .eq. -1e32) diff = 0.0 !trying to prevent explosion due to times not matching
  endif

  !When using MSIS
  if (.False.) then
    diff = -1e32
    call get_subsolar(CurrentTime, VernalTime, SSLon, SSLat)
    if (SSLon == 0) then
      SSLon = 1e-2
    endif

    if (SSLat == 0) then
      SSLat = 1e-2
    endif
    
    call LocationIndex(SSLon, SSLat, iiBlock, iLon, iLat, rLon, rLat)                         

    if (iiBlock .eq. 1) then                                                                
      call BlockAltIndex(400*1000.0, iiBlock, iLon, iLat, iAlt, rAlt)                       
      !MSIS rho                                                                               
      geo_lon = SSLon*180.0/pi                                                                
      geo_lat = SSLat*180.0/pi                                                                
      geo_lon = mod(SSLon*180.0/pi + 360.0, 360.0)                                           

      if (geo_lat < -90.0) then                                                    
        geo_lat = -180.0-geo_lat                                                          
        geo_lon = mod(geo_lon+180.0,360.0)                                 
      endif
                                                
      if (geo_lat >  90.0) then                                        
        geo_lat =  180.0-geo_lat                                                          
        geo_lon = mod(geo_lon+180.0,360.0)                                                  
      endif
      
      geo_alt = Altitude_GB(iLon, iLat, iAlt, iiBlock)/1000.0                                 
      geo_lst = mod(utime/3600.0+geo_lon/15.0,24.0)                                        

      CALL GTD7(iJulianDay,utime,geo_alt,geo_lat,geo_lon,geo_lst, &                      
              F107A,F107,AP,48,msis_dens,msis_temp)                                           
      msisDensity_SS = msis_dens(6)                                                       
      gitmDensity_SS = rLon*rLat*Rho(iLon, iLat, iAlt, iiBlock)       + &           
                       (1-rLon)*rLat*Rho(iLon+1, iLat, iAlt, iiBlock) + &        
                       rLon*(1-rLat)*Rho(iLon, iLat+1, iAlt, iiBlock) + &  
                       (1-rLon)*(1-rLat)*Rho(iLon+1, iLat+1, iAlt, iiBlock)               
      diff = msisDensity_SS - gitmDensity_SS  
    endif

    localVar = diff                                                                           
    call MPI_REDUCE(localVar, diff, 1, MPI_DOUBLE_PRECISION, &                             
                  MPI_MAX, 0, iCommGITM, iError)                        
  endif
  !if (lz == 2) then                                                                 
  !  localVar = diff2                                                         
  !  call MPI_REDUCE(localVar, diff2, 1, MPI_DOUBLE_PRECISION, &
  !                  MPI_MAX, 0, iCommGITM, iError)  
  !end if                                    
  !!!!!!!!!!!!
  ! Begin the first RCMR loop                                                 
  if (iProc == 0 .and. mod((istep-1)*Dts,Measure_Dts) == 0.0) then
     if (RCMROutType == "COND") then
       if (lz==1) then
            write(*,*) "Adding", diff, "to zp"
            Z = diff*1.0e10
            call RCPE(mu, mu, Z, Z*0.0)

            RCMR_est(1)  = abs(mu(1,1))
            write(*,*) "Thermal conductivity (iProc = 0):", &
                       rcmr_est0(1), RCMR_est(1)

            if (coefficientToEstimate == "AO2") then
              ThermalConduction_AO2 = RCMR_est(1) + RCMR_est0(1)
            else if (coefficientToEstimate == "AO") then
              ThermalConduction_AO = RCMR_est(1) + RCMR_est0(1)
            else if (coefficientToEstimate == "s") then
              ThermalConduction_s = RCMR_est(1) + RCMR_est0(1)
            endif
            !scattered = RCMR_est(1) + RCMR_est0(1)
            !ThermalConduction_AO2 = 2.0e-4
       endif
     endif
  !else
  !  return
  endif
  

  if (mod((istep-1)*Dts,Measure_Dts) == 0.0) then
    !call MPI_SCATTER(scattered, 1, mpi_double_precision, scatter, 1, &
    !                 mpi_double_precision, 0, iCommGITM, iError) !works
    !if (iProc .eq. 2) then
    !  write(*,*) "Before broadcasting...", ThermalConduction_AO2
    !endif
    !call MPI_BCAST(ThermalConduction_AO2, 1, mpi_double_precision, 0, &                   
    !               iCommGITM, iError)
    !if (iProc .eq. 2) then
    !  write(*,*) "Broadcasted value...", ThermalConduction_AO2
    !endif
    if (coefficientToEstimate == "AO") then
      !ThermalConduction_AO = scatter
      call MPI_BCAST(ThermalConduction_AO, 1, mpi_double_precision, 0, &
                   iCommGITM, iError)
    else if (coefficientToEstimate == "AO2") then
      call MPI_BCAST(ThermalConduction_AO2, 1, mpi_double_precision, 0, &
                   iCommGITM, iError)
    else if (coefficientToEstimate == "s") then
      call MPI_BCAST(ThermalConduction_s, 1, mpi_double_precision, 0, &
                   iCommGITM, iError)
    endif

    kk = kk + 1
  endif

end subroutine mainRCPE


  subroutine getNumLines(nLines)

    implicit none

    integer, intent(out) :: nLines
    integer :: io = 0
    character(len = 21) :: newFileName = "syntheticDataFile.txt"
    
    !BP: Check number of lines in data file for RCMR, then allocated  variables with correct 
    !    size.                                                              
    open (unit = 31, file = newFileName, status = "old")
    do
      read(31, *, iostat=io)
      if (io .ne. 0) then
        exit
      end if
      nLines = nLines + 1
    end do
    close(31)

    if (allocated(satYear) .eqv. .False.) then
      call arraySub(nLines)
    end if
    !!!!!!!!!!!!!!!!!!!!!
    
  end subroutine getNumLines

  subroutine arraySub(nLines)
  use ModGITM, only: nProcs
  INTEGER, intent(in) :: nLines
                                                                       
  allocate(satIndex(nLines), &
             satYear(nLines), &
             satMonth(nLines), &
             satDay(nLines), &
             satHour(nLines), &
             satMinute(nLines), &
             satSecond(nLines), &
             satMeanDensity(nLines), &
             sat2Density(nLines), &
             satLat(nLines), sat2Lat(nLines), &
             satLon(nLines), sat2Lon(nLines), &
             satAlt(nLines), sat2Alt(nLines), &
             time(nLines), date(nLines), &
             scattered(nProcs))

  scattered = 0.0 

  end subroutine arraySub
  

  subroutine readSyntheticData(nLines)
    implicit none
    integer, intent(in) :: nLines
    !integer :: nLines, n, i, lu
    integer :: n, i
    logical :: IsFirstTime = .true.
    character(len = 21) :: newFileName = "syntheticDataFile.txt"

    if (IsFirstTime ) then
      n = nLines
      open (unit = 31, file = newFileName, status = "old")

      ! DATA FORMAT - first line should be                       
      !1833  2002-12-31 00:05:00 1.352269e-12 2.84807e-12 2.1001695e-12 0.00036 0.00056 0.69    
      !write(*,*) "Reading synthetic data..."                                                   
      do i = 1, n
        if (lu == 1) then
          read(31, *) satIndex(i), date(i), time(i), satLon(i), satLat(i), &
                      satAlt(i), satMeanDensity(i)
        else if (lu == 2) then
          read(31, *) satIndex(i), date(i), time(i), satLon(i), satLat(i), &
                      satAlt(i), satMeanDensity(i), sat2Lon(i), sat2Lat(i), &
                      sat2Alt(i), sat2Density(i)
        end if

        read(date(i)(1:4), '(i4)') satYear(i)
        read(date(i)(6:7), '(i2)') satMonth(i)
        read(date(i)(9:10), '(i2)') satDay(i)

        read(time(i)(1:2), '(i2)') satHour(i)
        read(time(i)(4:5), '(i2)') satMinute(i)
        read(time(i)(7:8), '(i2)') satSecond(i)
      end do

      close(31)

      !if (iProc == 0) then                                                           
      !  write(*,*) date(1), date(n), satLon(1), satLon(n)                            
      !end if                           

      IsFirstTime = .false.

    endif

    return

end subroutine readSyntheticData

subroutine printDensity(nLines, diff, diff2, iiBlock)
  use ModTime
  use ModGITM
  !use arrayMod
  use ModMPI
                                         
  real :: LatFind, LonFind, AltFind
  real :: rLon, rLat, rAlt
  real :: gitmDensity
  integer :: iiLat, iiLon, iiAlt
  integer :: loc = -99
  integer :: n, i
  integer, intent(out) :: iiBlock
  !integer, intent(out) :: returnProc, returnProc2
  integer, intent(in) :: nLines
  logical :: exist
  real, intent(out) :: diff, diff2
  real :: percentDifference

  n = nLines
  diff = -1e32
  diff2 = -1e32

  do i = 1, n
    !if ((satYear(i) .eq. iTimeArray(1))   .and. &
    !    (satMonth(i) .eq. iTimeArray(2))  .and. &
    !    (satDay(i) .eq. iTimeArray(3))    .and. &
    !    (satHour(i) .eq. iTimeArray(4))   .and. &
    !    (satMinute(i) .eq. iTimeArray(5)) .and. &
    !    (iTimeArray(6) .eq. satSecond(i))) then

    if ((satYear(i) .eq. iTimeArray(1))   .and. &
        (satMonth(i) .eq. iTimeArray(2))  .and. &
        (satDay(i) .eq. iTimeArray(3))    .and. &
        (satHour(i) .eq. iTimeArray(4))   .and. &
        (satMinute(i) .eq. iTimeArray(5)) .and. &
        (abs(iTimeArray(6) - satSecond(i)) < 3)) then
      loc = i
      exit
    endif
  end do

  if (loc .eq. -99 .and. i .eq. n) then
    write(*,*) "ERROR: Time not found in satellite data"
    return
  endif
  
  ! ------------------------------------------------------                      
  !For first satellite (CHAMP)                                                  
  LonFind = satLon(loc)*pi/180.0
  LatFind = satLat(loc)*pi/180.0
  call LocationIndex(LonFind, LatFind, iiBlock, iiLon, iiLat, rLon, rLat)

  !Finds the satellite location and finds the corresponding mass density        
  !provided by GITM and prints it if it's within 5 seconds of CHAMP        
  if (iiBlock .eq. 1) then
    !if (satHour(loc) .eq. iTimeArray(4) .and. &
    !    satMinute(loc) .eq. iTimeArray(5) .and. &
    !    satSecond(loc) .eq. iTimeArray(6)) then

    if (satHour(loc) .eq. iTimeArray(4) .and. &
        satMinute(loc) .eq. iTimeArray(5) .and. &
        abs(iTimeArray(6) - satSecond(i)) < 3) then
        !iTimeArray(6) .eq. satSecond(loc)) then
      write(*,*) "GITM's time:", iTimeArray(2), "/", iTimeArray(3), "/", &
           iTimeArray(1), iTimeArray(4), ":", iTimeArray(5), ":", iTimeArray(6)
      
      !write(*,*) "Data file time:", satMonth(loc), "/", &                   
      !           satDay(loc), "/", satYear(loc), satHour(loc), ":", &          
      !           satMinute(loc), ":", satSecond(loc)                           

      !write(*,*) "At this time, CHAMP location is:", satLon(loc), &            
      !           satLat(loc), satAlt(loc)                      
      !write(*,*) "iiBlock:", iiBlock                                         
      !write(*,*) "iiLon:", iiLon                              
      !write(*,*) "iiLat:", iiLat                              
      !write(*,*) "rLon:", rLon                               
      !write(*,*) "rLat:", rLat
      
      AltFind = satAlt(loc)*1000.0
      call BlockAltIndex(AltFind, iiBlock, iiLon, iiLat, iiAlt, rAlt)
      !write(*,*) "iiAlt:", iiAlt                                                     
      !write(*,*) "GITM's found location:", Longitude(iiLon, iiBlock)*180/pi, & 
      !           Latitude(iiLat, iiBlock)*180/pi, &                            
      !           Altitude_GB(iiLon, iiLat, iiAlt, iiBlock)                     
      !GITM's density at satellite's location
      
      gitmDensity = rLon*rLat*rAlt*Rho(iiLon, iiLat, iiAlt, iiBlock)+ &
                    (1-rLon)*rLat*rAlt*Rho(iiLon+1, iiLat, iiAlt, iiBlock) + &
                    rLon*(1-rLat)*rAlt*Rho(iiLon, iiLat+1, iiAlt, iiBlock) + &
                    (1-rLon)*(1-rLat)*rAlt*Rho(iiLon+1, iiLat+1, iiAlt, iiBlock) + &
                    rLat*rLon*(1-rAlt)*Rho(iiLon, iiLat, iiAlt+1, iiBlock)+ &
                    (1-rLon)*rLat*(1-rAlt)*Rho(iiLon+1, iiLat, iiAlt+1, iiBlock) + &
                    rLon*(1-rLat)*(1-rAlt)*Rho(iiLon, iiLat+1, iiAlt+1, iiBlock) + &
                    (1-rLon)*(1-rLat)*(1-rAlt)*Rho(iiLon+1, iiLat+1, iiAlt+1, iiBlock)

      !Populate the difference array, z                        
      diff = satMeanDensity(loc) - gitmDensity
      !write(*,*) "Diff in function:", diff                                            
      !inquire(file = "RCMR_locations.txt", exist = exist)                            
      !if (exist) then                                                                
      !  open(33, file = "RCMR_locations.txt", status = "old", position = "append", & 
      !       action = "write")                                                       
      !else                                                        
      !  open(33, file = "RCMR_locations.txt", status = "new", action = "write")    
      !end if                           

      !111 FORMAT(I4, 1X, I2, 1X, I2, 1X, &                                           
      !           I2, 1X, I2, 1X, I2, 1X, &                 
      !           F8.2, 1X, F8.2, 1X, F8.2)                                                    
      !write(33,111) iTimeArray(1), iTimeArray(2), &                                  
      !              iTimeArray(3), iTimeArray(4), &                                  
      !              iTimeArray(5), iTimeArray(6), &                                  
      !              satLon(loc), satLat(loc), satAlt(loc)                            
      !close(33)                                                                       
    endif
  endif

  ! ------------------------------------------------------                                     
  !For second satellite GRACE              
  if (lu == 2) then
    LonFind = sat2Lon(loc)*pi/180.0
    LatFind = sat2Lat(loc)*pi/180.0
    call LocationIndex(LonFind, LatFind, iiBlock, iiLon, iiLat, rLon, rLat)
    if (iiBlock .eq. 1) then
      if (satHour(loc) .eq. iTimeArray(4) .and. &
          satMinute(loc) .eq. iTimeArray(5) .and. &
          iTimeArray(6) .eq. satSecond(loc)) then
        !write(*,*) "At this time, GRACE location is:", sat2Lon(loc), &         
        !           sat2Lat(loc), sat2Alt(loc)                                  
        !write(*,*) "iiBlock:", iiBlock                                         
        !write(*,*) "iiLon:", iiLon                                             
        !write(*,*) "iiLat:", iiLat                                             
        !write(*,*) "rLon:", rLon                                            
        !write(*,*) "rLat:", rLat   

        AltFind = sat2Alt(loc)*1000.0
        call BlockAltIndex(AltFind, iiBlock, iiLon, iiLat, iiAlt, rAlt)
        !write(*,*) "iiAlt:", iiAlt                           
        !write(*,*) "GITM's found location:", Longitude(iiLon, iiBlock)*180/pi,               !           Latitude(iiLat, iiBlock)*180/pi, &                       
        !           Altitude_GB(iiLon, iiLat, iiAlt, iiBlock)                                  
        !GITM's density at satellite's location
        gitmDensity = rLat*rLon*rAlt*Rho(iiLon, iiLat, iiAlt, iiBlock)+ &
                      (1-rLon)*rLat*rAlt*Rho(iiLon+1, iiLat, iiAlt, iiBlock) + &
                      rLon*(1-rLat)*rAlt*Rho(iiLon, iiLat+1, iiAlt, iiBlock) + &
                      (1-rLon)*(1-rLat)*rAlt*Rho(iiLon+1, iiLat+1, iiAlt, iiBlock) + &
                      rLat*rLon*(1-rAlt)*Rho(iiLon, iiLat, iiAlt+1, iiBlock)+ &
                      (1-rLon)*rLat*(1-rAlt)*Rho(iiLon+1, iiLat, iiAlt+1, iiBlock) + &
                      rLon*(1-rLat)*(1-rAlt)*Rho(iiLon, iiLat+1, iiAlt+1, iiBlock) + &
                      (1-rLon)*(1-rLat)*(1-rAlt)*Rho(iiLon+1, iiLat+1, iiAlt+1, iiBlock)

        !write(*,*) "GITM's density:", gitmDensity                                    
	!write(*,*) "Truth data density:", sat2Density(loc)             

	!Populate the difference array, z                               
	diff2 = sat2Density(loc) - gitmDensity
      !returnProc2 = iProc                                                             
      endif
    endif
  endif

end subroutine printDensity





  subroutine clean_mod_rcpe
    if(.not.allocated(satIndex)) RETURN
    deallocate(satIndex)
    deallocate(satYear)
    deallocate(satMonth)
    deallocate(satDay)
    deallocate(satHour)
    deallocate(satMinute)
    deallocate(satSecond)
    deallocate(satMeanDensity)
    deallocate(sat2Density)
    deallocate(satLat)
    deallocate(sat2Lat)
    deallocate(satLon)
    deallocate(sat2Lon)
    deallocate(satAlt)
    deallocate(sat2Alt)
    deallocate(time)
    deallocate(date)
    deallocate(scattered)
    
  end subroutine clean_mod_rcpe
  

END MODULE ModRCPE


