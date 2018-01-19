!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------------
! $Id: RCMR_routines.f90,v 1.7 2017/10/30 14:05:36 ridley Exp $
!
! Author: Asad, UMichigan, 2/2013
!
! Comments: Routines to run RCMR data assimilation
!
! AGB 3/31/13: Added comments, changed variable names to be more descriptive,
!              changed subroutine format for consistency with GITM files,
!              removed test output files, removed unused variables, streamlined
!              MPI calls.
! AGB 10/23/13: Adapted to allow driving of photoelectron heating efficiency
!-----------------------------------------------------------------------------

subroutine run_RCMR

  use ModInputs
  use ModTime
  use ModGITM
  use ModMpi
  use ModRCMR
  use ModSatellites, only: SatCurrentDat, SatAltDat, nRCMRSat, RCMRSat
  use ModEUV, only: sza	
  implicit none

  integer :: status(MPI_STATUS_SIZE), iError, ii_loop
  real :: output_est(1,1), ulimit, llimit, dummy_TEC_calculated, dummy
  double precision :: localVar
  !double precision, dimension(12) :: TEC_calculated = 0     !ANKIT: holds the value of TEC at AA 
  double precision, dimension(:), allocatable :: TEC_calculated     !ANKIT: holds the value of TEC at AA 
  double precision :: TEC_calculated_P1 = 0  !ANKIT: holds the value of TEC at ND
  double precision :: u_sat_level = 0
  allocate(TEC_calculated(lz))

  call MPI_BARRIER(iCommGITM, iError)


  if (RCMROutType == "EDC") then
     !if (iproc == 2) then
     !   !call calc_single_vtec(5, 5, 1, TEC_calculated)
     !   !Ankit5Feb2015 - Output TEC at 12 locations
     !   call calc_single_vtec_interp(1.0* 3.141592653589793/180, 0.0, TEC_calculated(1))
     !   call calc_single_vtec_interp(30.0* 3.141592653589793/180, 0.0, TEC_calculated(2))
     !   call calc_single_vtec_interp(60.0* 3.141592653589793/180, 0.0, TEC_calculated(3))
     !   call calc_single_vtec_interp(90.0* 3.141592653589793/180, 0.0, TEC_calculated(4))
     !   call calc_single_vtec_interp(120.0* 3.141592653589793/180, 0.0, TEC_calculated(5))
     !   call calc_single_vtec_interp(150.0* 3.141592653589793/180, 0.0, TEC_calculated(6))
     !end if
     !if (iproc == 3) then
     !   !call calc_single_vtec(3, 5, 1, TEC_calculated_P1)
     !   call calc_single_vtec_interp(180.0* 3.141592653589793/180, 0.0, TEC_calculated(7))
     !   call calc_single_vtec_interp(210.0* 3.141592653589793/180, 0.0, TEC_calculated(8))
     !   call calc_single_vtec_interp(240.0* 3.141592653589793/180, 0.0, TEC_calculated(9))
     !   call calc_single_vtec_interp(270.0* 3.141592653589793/180, 0.0, TEC_calculated(10))
     !   call calc_single_vtec_interp(300.0* 3.141592653589793/180, 0.0, TEC_calculated(11))
     !   call calc_single_vtec_interp(330.0* 3.141592653589793/180, 0.0, TEC_calculated(12))
     !end if

     if (iproc == TEC_proc) then
        do ii_loop = 1,lz
           call calc_single_vtec_interp(TEC_location(index_lz_tec(ii_loop),1), &
                TEC_location(index_lz_tec(ii_loop) ,2), TEC_calculated(ii_loop))
           !call calc_single_vtec_interp(dummy, dummy, dummy)
           !TEC_calculated(ii_loop) = dummy_TEC_calculated
           !! call calc_single_vtec_interp(TEC_location(index_lz_tec(ii_loop),1), TEC_location(index_lz_tec(ii_loop) ,2), 0.1)
        end do
     endif


  else
     call MPI_BARRIER(iCommGITM, iError)
     localVar = SatAltDat(RCMRSat(1))
     call MPI_REDUCE(localVar, SatAltDat(RCMRSat(1)), 1, MPI_DOUBLE_PRECISION, &
          MPI_MAX, 0, iCommGITM, iError)
     if (iProc==0) then
        do mm=0,(nProcs-1)
           if (Sat_Proc(mm + 1) .ne. 0) then
              SatAltDat(RCMRSat(1)) = Sat_Proc(mm + 1)   !Ankit: what is this doing??
           end if
        end do
     end if
  endif
 
  !ANKIT: Added MPI_Bcast to update on all processors
  !ANKIT:5Feb15 - Added TEC bcasting from 12 locations
  !call MPI_Bcast( TEC_calculated(1:6)   , 6, MPI_DOUBLE_PRECISION,2, iCommGITM, iError)
  !call MPI_Bcast( TEC_calculated(7:12)  , 6, MPI_DOUBLE_PRECISION,3, iCommGITM, iError)
  
  call MPI_Bcast( TEC_calculated  , lz, MPI_DOUBLE_PRECISION, TEC_proc, iCommGITM, iError) !Ankit25Feb2015: Send TEC to all procs

  !call MPI_Bcast( TEC_calculated_P1, 1, MPI_DOUBLE_PRECISION,2, iCommGITM, iError)
  !  write (*,*) iproc, TEC_calculated, TEC_calculated_P1
  !  write(*,*) istep, iproc, "VTEC calculated", TEC_calculated
  
  ! AGB: Initialize the output estimate
  if(RCMROutType == "F107") then
     output_est = f107_est
     llimit     = 70.0
     ulimit     = 400.0
  else if(RCMROutType == "PHOTOELECTRON") then
     output_est = PhotoElectronHeatingEfficiency_est
     llimit     = 0.02
     ulimit     = 0.2
  else if (RCMROutType == "EDC") then           !ANKIT: Eddy diffusion coefficient hard limits
     output_est = EDC_est
     llimit     = 500
     ulimit     = 3000
  else
     write (*,*) "ERROR: unknown RCMR output type", RCMROutType
  end if
  
  
  ! Begin the first RCMR loop
  if (mod((istep-1)*Dts,Measure_Dts) == 0.0 .AND. TRUTH_or_ID==1 .AND. iProc==0) then 
     if (RCMROutType == "EDC") then                                                 !ANKIT 6 Oct 2014
        !zp(1,ustep) = alpha*(TEC_calculated    - TEC_true(iStep,1))                 !!!!!!!!!!!!!!!!!!!!!! FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !zp(2,ustep) = alpha*(TEC_calculated_P1 - TEC_true(iStep,2))
        !ANKIT5Feb2015 - calculate z at 12 locations
        do ii_loop = 1,lz
           !zp(ii_loop,ustep) = alpha*(TEC_calculated(ii_loop) - TEC_true(iStep,ii_loop))
           zp(ii_loop,ustep) = alpha*(TEC_calculated(ii_loop) - TEC_true(iStep,index_lz_tec(ii_loop))) !Ankit25Feb2015: added choice of tec station
        end do
        !write(*,*) ustep, iStep, TEC_calculated, TEC_true(istep,1), zp(1,ustep)
        !write(*,*) ustep, iStep, TEC_calculated_P1, TEC_true(istep,2), zp(2,ustep)
     else
        zp(:,uStep) = SatAltDat(RCMRSat(1)) - SatCurrentDat(RCMRSat(1)) 
     endif
      
     ! Set the averaged error
     if (RCMROutType == "EDC") then          !ANKIT 4 Nov 2014
        !zav(:,ustep) = zp(:,uStep)           ! Remove filter
        !y(1,uStep) = TEC_calculated          ! 1.0e00*zav(1,uStep)     ! Remove scaling
        !y(2,uStep) = TEC_calculated_P1       ! 1.0e00*zav(1,uStep)     ! Remove scaling
        y(:,ustep) = TEC_calculated         !ANKIT5Feb2015
        if (uStep <= 120) then
           zav(:,ustep) = sum(zp(:,1:uStep))/uStep
        else
           zav(:,ustep) = sum(zp(:,uStep-120+1:uStep))/120
        end if
        z(:,uStep) = zav(:,uStep)           !Ankit 25Jan15 - Adds moving average filter to Z
        z(:,uStep) = zp(:,uStep)            !Ankit 25Jan15 - Removes the averaging from Z
     else
        if (uStep <= 120) then
           zav(:,ustep) = sum(zp(:,1:uStep))/uStep
        else
           zav(:,ustep) = sum(zp(:,uStep-120+1:uStep))/120
        end if
        y(1,uStep) = 1.0e12*zav(1,uStep)
        z(:,uStep) = 1.0e12*zav(1,uStep)
     endif

     if (ustep <= l_dim) then
        y_mat(:,l_dim-ustep+1:l_dim) = y(:,1:ustep)
        z_mat(:,l_dim-ustep+1:l_dim) = z(:,1:ustep)
        if (ustep>1) then
           u_mat(:,l_dim-ustep+2:l_dim) = u(:,1:ustep-1)   
        end if
     else if (ustep > l_dim) then
        y_mat = y(:,ustep-(l_dim)+1:ustep)
        u_mat = u(:,ustep-(l_dim):ustep-1)
        z_mat = z(:,ustep-(l_dim)+1:ustep)
     end if

     if (ustep<=Pc) then
        Pcc = ustep-1
        if (Pcc < 0) then
           Pcc=0
        end if
     else
        Pcc=Pc
     end if
     lu2=lu*lu

     if (ustep>C_on) then
        control_on = 1
     end if

     if (RCMROutType == "EDC") then          ! ANKIT 6 Nov 2014
        !CALL RCMR_onestep(u_out, ustep, EDC_est, zp(:,ustep), y(:,ustep))
        CALL RCMR_onestep(u_out, ustep, EDC_est, z(:,ustep), y(:,ustep), output_est)
        ! u_out is currently 1 by 1. Need to make it lu by 1
        ! EDC_est is u_in, that is the last value of u. Is 1 by 1
        ! zp(:,ustep) is z_in. Should be lz by 1. Is 1 by 1
        ! y_in = 0
     else
        CALL RCMR_Function(llimit, ulimit, Nc, ustep, s_dhat, lz, lu, lu2, ly, &
             l_dim, Pc, Pcc, control_on, dbuffer, C_on, dhat, eta, lambda, usum, &
             T,  R2, y_mat, z_mat, u_mat, P1, u_out, theta1, UB)
     endif

     ! write(*,*) istep, ustep, u_out(1,1)     

     ! Enforce the realistic physical limitations
     if (u_out(1,1) < llimit) then
        u_out(1,1) = llimit
     end if
     if (u_out(1,1) >= ulimit) then
        u_out(1,1) = ulimit 
     end if

     up(1,ustep) = u_out(1,1)

     ! Find the averaged error
     
     !! 4 November 2014: Add rate saturation
     u_sat_level = 5 ! rate saturation 
     !!if (RCMROutType == "EDC") then          ! ANKIT 4 Nov 2014
     if (.false.) then                               ! ANKIT 16 Jan 2014, Removing the saturation, and adding filtering
        if ((abs(up(1,uStep) - up(1,uStep-1)) > u_sat_level) .and. (ustep > C_on+2))  then
           if (up(1,uStep) > up(1,uStep-1)) then
              up(1,uStep) = up(1,uStep-1) + u_sat_level
           else
              up(1,uStep) = up(1,uStep-1) - u_sat_level
           endif
        endif
        u(:,ustep) = up(:,uStep)             ! Remove filter
     else
        if (uStep <= uFiltLength) then
           u(:,ustep) = sum(up(:,1:uStep))/uStep
        else
           u(:,ustep) = sum(up(:,uStep-uFiltLength+1:uStep))/uFiltLength
        end if
     endif

     ! write(*,*) istep, ustep, up(1,ustep), u(1,ustep), output_est






     ! Save the output if GITM is initialized
     ! AGB Question: would it not be better to ask this before calculating
     !               a new output estimate?

     if (ustep > C_on) then
        output_est = u(1, ustep)
     else if (ustep <= C_on) then
        u(1,ustep) = output_est(1,1)
     end if

     if(RCMROutType == 'F107') then
        f107_est  = output_est(1,1)
        f107a_est = f107_est
     else if(RCMROutType == "PHOTOELECTRON") then
        PhotoElectronHeatingEfficiency_est = output_est(1,1)
     elseif(RCMROutType == "EDC") then !ANKIT 6 Oct 2014
        EDC_est = output_est
        !write (*,*) "EDC estimate is ", EDC_est, ustep
     end if

     ustep = ustep + 1
  end if

  ! END of the first RCMR loop

  ! Send the F10.7 Estimate to all the different processors
  ! AGB Question: why use MPI_SCATTER and not MPI_Bcast?
  if (iProc==0) then
     scattered(1) = output_est(1,1)
  end if

  call MPI_SCATTER(scattered, 1, mpi_double_precision, scatter, 1, &
       mpi_double_precision, 0, iCommGITM, iError)

  if(RCMROutType == 'F107') then
     f107_est  = scatter    
     f107a_est = f107_est
  else if(RCMROutType == 'PHOTOELECTRON') then
     PhotoElectronHeatingEfficiency_est = scatter
  else if(RCMROutType == 'EDC') then !ANKIT 6 Oct 2014
     EDC_est = scatter
     if (iProc == TEC_proc) call write_TEC_data
  end if

end subroutine run_RCMR

subroutine RCMR_Function(llimit, ulimit, Nc, ustep, s_dhat, lz, lu, lu2, ly, &
     l_dim, Pc, Pcc, control_on, dbuffer, C_on, dhat, eta, lambda, usum, T, &
     R2, y_mat, z_mat, u_mat, P1, u_out, theta1, UB)

  use ModBlasLapack, only: LAPACK_getrf, LAPACK_getrs
  implicit none

  ! *
  !  Input variables
  ! *

  real, intent(in) :: llimit, ulimit
  integer, intent(in) :: Nc, ustep, s_dhat, lu, ly, lz, l_dim, Pc, Pcc
  integer, intent(in) :: control_on, dbuffer, C_on
  integer, dimension(1,s_dhat), intent(in) :: dhat

  double precision :: eta
  double precision, dimension(1,1), intent(in) :: lambda, usum
  double precision, dimension(s_dhat*lz, lu), intent(in) :: T
  double precision, dimension(s_dhat*lz,s_dhat*lz), intent(in) :: R2
  double precision, dimension(ly, l_dim), intent(in) :: y_mat ! y_data
  double precision, dimension(lz, l_dim), intent(in) :: z_mat ! z_data
  double precision, dimension(lu, l_dim), intent(in) :: u_mat ! u_data

  ! *
  !  Input and output variables
  ! *

  double precision, dimension(lu*Nc+ly*Nc,lu*Nc+ly*Nc) :: P1

  ! AGB MOVED FROM INPUT

  integer :: iii, lb_z, errorhandler, lu2
  double precision, dimension(1,1) :: xi_RLS, r_RLS
  double precision, dimension(lu*Nc+ly*Nc,1) :: pi_RLS, k_RLS
  double precision, dimension(lu, 1) :: umod

  ! *
  !  Work variables
  ! *

  integer, dimension(lu) :: PIVOTARRAY
  integer, dimension(1,s_dhat) :: dhatfl ! used with flipped dhat

  double precision intr_a
  double precision, dimension(lu2) :: WORK
  double precision, dimension(1,1) :: intr_b
  double precision, dimension(s_dhat*lz,1) :: zmod ! mat dim r=lz, c=size(dhat, control_on)
  double precision, dimension(lu,1) :: Us
  double precision, dimension(lz,1) :: zr
  double precision, dimension(lu,lu) :: To_invert, Idty, Inverted !FA=lu, SA=lu
  double precision, dimension(ly,Nc) :: py
  double precision, dimension(ly*(Nc),1) :: pym
  double precision, dimension(lu, Nc) :: pu
  double precision, dimension(lu*Nc,1) :: pum
  double precision, dimension(lu*Nc+ly*Nc,1) :: H1
  double precision, dimension(1,lu*Nc+ly*Nc) :: Hm
  double precision, dimension(1,lu*Nc+ly*Nc) :: K1

  ! *
  !  Output variables
  ! *

  double precision, dimension(1,1), intent(out) :: UB
  double precision, dimension(lu,1), intent(out) :: u_out
  double precision, dimension(lu,Nc*ly+Nc*lu), intent(out) :: theta1

  ! Identity matrix
  iDty     = 0.0;
  Inverted = 0.0;
  do iii=1,lu
     Idty(iii,iii) = 1.0
     Inverted(iii,iii) = 1.0
  end do

  ! Initialization
  zmod = 0.0

  ! *
  !  In this section, we calculate Uhat
  ! *

  if (ustep > maxval(dhat)) then
     zmod(1:lz,1) = z_mat(:,l_dim)
     dhatfl(1,:)  = dhat(1,s_dhat:1:-1)

     if (s_dhat>1) then
        make_zmod: do iii=1,(s_dhat-1)
           lb_z = dhat(1,s_dhat)
           zmod(1+(iii*lz):1+(iii*lz)+lz-1,1) = z_mat(:, &
                l_dim-( dhat(1,s_dhat)-dhatfl(1,iii+1))) 
        end do make_zmod
     end if

     umod(:,1) = u_mat(:,l_dim-dhat(1,s_dhat)+1);
     zr(:,1)   = z_mat(:,size(z_mat,2));

     To_Invert = matmul(matmul(transpose(T),R2),T) &
          + eta*matmul(matmul(transpose(zr),zr),Idty)

     ! invert this
     call LAPACK_getrf(lu, lu, To_Invert, lu, PIVOTARRAY, errorHandler)
     CALL LAPACK_getrs('n', lu, lu, To_Invert, lu, PIVOTARRAY, Inverted, lu, &
          errorHandler)

     Us = matmul(matmul(matmul(Inverted,transpose(T)),R2), &
          (-zmod+matmul(T,umod)))

     UB(1,1) = Us(1,1)
     if (Us(1,1)>=ulimit) then
        Us(1,1) = ulimit
     else if(Us(1,1)<=llimit) then
        Us(1,1) = llimit
     end if

     if (ustep > dbuffer) then
        !Strictly proper
        py = y_mat(:,l_dim-Nc-dhat(1,s_dhat):l_dim-dhat(1,s_dhat)-1)
        !Exactly proper
        !py = y_mat(:,l_dim-Nc-dhat(1,s_dhat)+1:l_dim-dhat(1,s_dhat))

        pym = reshape(py, shape(pym))
        pu  = u_mat(:,(l_dim+1)-Nc-dhat(1,s_dhat):(l_dim+1)-dhat(1,s_dhat)-1)

        if(minval(pu) == 0.0) then
           pu = 0.0
        end if

        pum = reshape(pu, shape(pum))

        H1(1:ly*(Nc),1) = pym(:,1)
        H1(ly*(Nc)+1:ly*(Nc)+lu*Nc,1) = pum(:,1)
        Hm = matmul(transpose(H1), P1)

        pi_RLS = matmul(P1,H1)
        r_RLS = 1/(lambda(1,1)+matmul(transpose(H1),pi_RLS))
        k_RLS = r_RLS(1,1)*pi_RLS
        xi_RLS = Us-matmul(theta1, H1)
        theta1 = theta1 + matmul(xi_RLS,transpose(k_RLS))
        P1 = (P1 - matmul(k_RLS,transpose(pi_RLS))) 

        py = y_mat(:,l_dim-Nc:l_dim-1) !Strictly proper
        !py = y_mat(:,l_dim-Nc+1:l_dim) !Exactly proper
        pym = reshape(py, shape(pym))

        pu = u_mat(:,l_dim-Nc+1:l_dim)
        pum = reshape(pu, shape(pum))
	
        H1(1:ly*(Nc),1) = pym(:,1)
        H1(ly*(Nc)+1:ly*(Nc)+lu*Nc,1) = pum(:,1)

        if (control_on == 1) then
           u_out =  matmul(theta1,H1)
        else
           u_out = 0.0
        end if
     else
        u_out = 0.0
        theta1 = theta1
     end if	
  else
     u_out = 0.0
     theta1 = theta1
  end if

  return
end subroutine RCMR_Function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! NEW ONE STEP RCAC 
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RCMR_onestep(u_out1, kk, u_in, z_in, y_in, ESTIMATE)
  use ModRCMR 
  use ModBlasLapack, only: LAPACK_getrf, LAPACK_getrs
  implicit none 
  integer iii, jjj
  integer :: errorhandler
  DOUBLE PRECISION :: det
  DOUBLE PRECISION, DIMENSION(lz,1) :: z_hat
  integer, intent(in) :: kk
  DOUBLE PRECISION, intent(in), DIMENSION(lu,1) :: u_in
  DOUBLE PRECISION, intent(in), DIMENSION(lz,1) :: z_in
  DOUBLE PRECISION, intent(in), DIMENSION(ly,1) :: y_in
  DOUBLE PRECISION, intent(in), DIMENSION(lu,1) :: ESTIMATE
  DOUBLE PRECISION, intent(out), DIMENSION(lu,1) :: u_out1

  logical :: exist
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!REGRESSOR
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Build regressor here
  u_h(:,2:NC)     = u_h(:,1:NC-1)
  z_h(:,2:NC+1)   = z_h(:,1:NC)
  y_h(:,2:NC+1)   = y_h(:,1:NC)
  DO iii = 1,lu
     u_h(iii,1)        = u_in(iii,1)
  END DO
  DO iii = 1,lz
     z_h(iii,1)        = z_in(iii,1)
  END DO
  DO iii = 1,ly
     y_h(iii,1)        = y_in(iii,1)
  END DO
  
  ! Stack up u's and performance here
  !IF (kk > Nc) THEN
  IF (kk > 0) THEN  ! Initialized the regressor before controller
     Uphi = u_h(:,1:NC)
     !! To have SP controller, indices run from 2 to Nc+1. 
     !! P controller will have indices running from 1 to Nc.
     if (RegZ == 1) THEN
        Vphi = z_h(:,2:NC+1)  
     ELSE
        Vphi = y_h(:,2:NC+1)
     ENDIF
     UphiVec = reshape(Uphi, shape(UphiVec))
     VphiVec = reshape(Vphi, shape(VphiVec))
     
     phi_reg(1:lu*NC,:) = UphiVec
     if (RegZ == 1) THEN
        phi_reg(lu*NC+1:(lz+lu)*NC,:) = VphiVec
     ELSE
        phi_reg(lu*NC+1:(ly+lu)*NC,:) = VphiVec
     ENDIF
     !! call kronecker(PHI, phi_reg)
     DO iii = 1,lu
        DO jjj = 1,lphi
           PHI(iii, jjj+(lu-1)*lphi) = phi_reg(jjj,1)
        ENDDO
     ENDDO
  ENDIF

!  write(*,*) "Phi regressor is ", phi_reg
!  write(*,*) "PHI ", PHI

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !! FILTERS 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Ubar(:,2:nf_u)             = Ubar(:,1:nf_u-1)
  Zbar(:,2:nf_z+1)           = Zbar(:,1:nf_z)
  PHIbar(lu+1:lu*(nf_u+1),:) = PHIbar(1:lu*nf_Nu,:)   !! Fixed lu
  DO iii = 1,lu
     Ubar(iii,1)        = u_in(iii,1)
  END DO
  DO iii = 1,lz
     Zbar(iii,1)        = z_in(iii,1)
  END DO
  DO jjj = 1,lu
     DO iii = 1,ltheta
        PHIbar(jjj,iii) = PHI(jjj,iii)
     END DO
  END DO

!  z_filtered      = reshape(matmul(Nz, reshape(Zbar, (/lz*(nf_z+1)/))) - matmul(Dz, reshape(Zfbar, (/lz*nf_Dz/))), (/lz,1/))
!  u_filtered      = reshape(matmul(Nu, reshape(Ubar, (/lu*nf_Nu/))) - matmul(Du, reshape(Ufbar, (/lz*nf_Du/))), (/lz,1/))
!  phi_filtered    = matmul(Nphi, PHIbar) - matmul(Dphi, PHIfbar)

  z_filtered      = reshape(matmul(Nz, reshape(Zbar, (/lz*(nf_z+1)/))) - matmul(Dz, reshape(Zfbar, (/lz*nf_z/))), (/lz,1/))
  u_filtered      = reshape(matmul(Nu, reshape(Ubar, (/lu*nf_u/)))     - matmul(Du, reshape(Ufbar, (/lz*nf_u/))), (/lz,1/))
  phi_filtered    = matmul(Nphi, PHIbar) - matmul(Dphi, PHIfbar)




  !!print*, z_filtered

  Ufbar(:,2:nf_u)     = Ufbar(:,1:nf_u-1)
  Zfbar(:,2:nf_z)     = Zfbar(:,1:nf_z-1)
  PHIfbar(lz+1:lz*nf_u,:) = PHIfbar(1:lz*nf_u-lz,:)   !! fix this, this starts from lz

  DO iii = 1,lz
     Ufbar(iii,1)        = u_filtered(iii,1)
  END DO
  DO iii = 1,lz
     Zfbar(iii,1)        = z_filtered(iii,1)
  END DO
  DO jjj = 1,lz
     DO iii = 1,ltheta
        PHIfbar(jjj,iii) = PHI_filtered(jjj,iii)
     END DO
  END DO
  
  z_hat = z_filtered + matmul(PHI_filtered, reshape(theta_h(:,1), (/ltheta,1/))) - u_filtered

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !! THE controller
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  IF ( (kk > Nc) .and. (kk > C_on) ) THEN         !! This has to be fixe
  IF (kk > C_on) THEN
     IF (W_Ru == 0.0) THEN
        XX = PHI_filtered
        zp1 = z_filtered - u_filtered
        Rp = Rz
        call identity(Rbarinv, lv, 1/(W_Rz+W_Ruf))
     ELSE
        XX(1:lz, :) = PHI_filtered
        XX(lz+1:lz+lu, :) = PHI   !! Fix this, this is supposed to be pHI
        zp1(1:lz, :) = z_filtered - u_filtered
        zp1(lz+1:lz+lu, :) = 0
        Rp(1:lz, 1:lz) = Rz
        Rp(lz+1:lz+lu, lz+1:lz+lu) = Ru
        DO iii = 1,lz
           Rbarinv(iii,iii) = 1/(W_Rz+W_Ruf)
        ENDDO
        DO iii = lz+1,lz+lu
           Rbarinv(iii,iii) = 1/W_Ru
        ENDDO
     ENDIF
     
     Ginv = (lambda1*Rbarinv  + matmul(matmul(XX, PP),transpose(XX)))
     tau  = 0.0
     Do iii = 1,lv
        tau(iii,iii) = 1.0
     enddo

     call LAPACK_getrf(lv, lv, Ginv, lv ,IPIV, INFO)
     CALL LAPACK_getrs('n', lv, lv, Ginv, lv, IPIV, Tau, lv, INFO)
     

     !     det = Ginv(1,1)*Ginv(2,2) - Ginv(1,2)*Ginv(2,1)
     !     tau(1,1) = Ginv(2,2)/det
     !     tau(1,2) = -Ginv(1,2)/det
     !     tau(2,1) = -Ginv(2,1)/det
     !     tau(2,2) = Ginv(1,1)/det
     !     Tau     = 1/Ginv !WORK


     PpXxpTauInv = matmul(matmul(PP, transpose(XX)), Tau)
     XxTheta     = reshape(matmul(XX, theta_h(:,1)), (/size(XX,1),1/))
     Rzp         = matmul(matmul(Rbarinv, Rp), zp1)
     !    PRdelta     = matmul(PP, Rdelta)
     !    dTheta      = reshape(theta_h(:,2) - theta_h(:,1), (/ltheta,1/))
     thetaout    = reshape(theta_h(:,1), (/ltheta,1/)) - &
          matmul(PpXxpTauInv, (Xxtheta+Rzp))
          !    matmul( matmul(PRdelta, (matmul(transpose(XX), transpose(PpXxpTauInv))-Iltheta)), dTheta)
     
     PP          = (PP - matmul(matmul(matmul(PP, transpose(XX)), Tau), matmul(XX, PP)))/lambda1

     if (.false.) then
        write(*,*) "z_hat is ", z_hat
        write(*,*) "Zf is ", z_filtered
        write(*,*) "Uf is ",  u_filtered
        write(*,*) "Phif is ", phi_filtered
        
        write(*,*) "XX \n", XX
        write(*,*) "Zp1 \n", zp1
        write(*,*) "Rp \n", Rp
        write(*,*) "Rbarinv \n", Rbarinv
        write(*,*) "Ginv \n", Ginv
        write(*,*) "Tau \n", Tau
        write(*,*) "PpXxpTauInv \n", PpXxpTauInv
        write(*,*) "Xxtheta \n", XxTheta
        write(*,*) "Rzp \n", Rzp
        write(*,*) "thetaout \n", thetaout
        write(*,*) "PP \n", PP
     endif

        theta_h(:,2) = theta_h(:,1)
     DO iii = 1,ltheta
        theta_h(iii,1)        = thetaout(iii,1)
     END DO
     u_out1 = matmul(PHI, thetaout)

     !     write(121) kk, transpose(u_in), transpose(z_in), transpose(y_in), transpose(u_out1), transpose(thetaout)
     !if (mod(kk, 1) == 0) then
     !   write(579, "(300F12.6)") thetaout  ! ltheta = 300
     !endif

     !if (mod(kk, 10) == 0) then  ! ANKIT: this is a problem, I want to write at the last step :(
     !   rewind(580)
     !   do iii = 1,ltheta
     !      write(580, "(300F16.9)") PP(iii, :)  ! ltheta = 300
     !   enddo
     !endif
     
  else
     u_out1 = ESTIMATE
  ENDIF
  

  !Ankit 5May2015: The following code writes theta to theta.dat
  IF (kk > C_on) THEN  
     inquire(file="theta.dat", exist=exist)
     if (exist) then
        open(12345, file="theta.dat", status="old", position="append", action="write")
     else
        open(12345, file="theta.dat", status="new", action="write")
     end if
     
     
     do iii = 1,ltheta
        write(12345, *) thetaout(iii,1)
     end do
  
     close(12345)
  endif

end subroutine RCMR_onestep

