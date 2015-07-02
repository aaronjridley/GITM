!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------------
! $Id: RCMR_routines.f90,v 1.5 2013/10/24 18:53:28 agburr Exp $
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

  integer :: status(MPI_STATUS_SIZE), iError
  real :: output_est, ulimit, llimit
  double precision :: localVar

  call MPI_BARRIER(iCommGITM, iError)
  localVar = SatAltDat(RCMRSat(1))
  call MPI_REDUCE(localVar, SatAltDat(RCMRSat(1)), 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, 0, iCommGITM, iError)

  if (iProc==0) then
     do mm=0,(nProcs-1)
        if (Sat_Proc(mm + 1) .ne. 0) then
           SatAltDat(RCMRSat(1)) = Sat_Proc(mm + 1)
        end if
     end do
  end if

  ! AGB: Initialize the output estimate
  if(RCMROutType == "F107") then
     output_est = f107_est
     llimit     = 70.0
     ulimit     = 400.0
  else if(RCMROutType == "PHOTOELECTRON") then
     output_est = PhotoElectronHeatingEfficiency_est
     llimit     = 0.02
     ulimit     = 0.2
  else
     write (*,*) "ERROR: unknown RCMR output type", RCMROutType
  end if

  ! Begin the first RCMR loop
  if (mod((istep-1)*Dts,Measure_Dts) == 0.0 .AND. TRUTH_or_ID==1 .AND. &
       iProc==0) then
        
     zp(:,uStep) = SatAltDat(RCMRSat(1)) - SatCurrentDat(RCMRSat(1))

     ! Set the averaged error
     if (uStep <= 120) then
        zav(:,ustep) = sum(zp(:,1:uStep))/uStep
     else
        zav(:,ustep) = sum(zp(:,uStep-120+1:uStep))/120
     end if
        
     y(1,uStep) = 1.0e12*zav(1,uStep)
     z(:,uStep) = 1.0e12*zav(1,uStep)

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

     CALL RCMR_Function(llimit, ulimit, Nc, ustep, s_dhat, lz, lu, lu2, ly, &
          l_dim, Pc, Pcc, control_on, dbuffer, C_on, dhat, eta, lambda, usum, &
          T,  R2, y_mat, z_mat, u_mat, P1, u_out, theta1, UB)

     ! Enforce the realistic physical limitations
     if (u_out(1,1) < llimit) then
        u_out(1,1) = llimit
     end if
     if (u_out(1,1) >= ulimit) then
        u_out(1,1) = ulimit 
     end if

     up(1,ustep) = u_out(1,1)

     ! Find the averaged error

     if (uStep <= 90) then
        u(:,ustep) = sum(up(:,1:uStep))/uStep
     else
        u(:,ustep) = sum(up(:,uStep-90+1:uStep))/90
     end if

     ! Save the output if GITM is initialized
     ! AGB Question: would it not be better to ask this before calculating
     !               a new output estimate?

     if (ustep > C_on) then
        output_est = u(1, ustep)
     else if (ustep <= C_on) then
        u(:,ustep) = output_est
     end if

     if(RCMROutType == 'F107') then
        f107_est  = output_est
        f107a_est = f107_est
     else if(RCMROutType == "PHOTOELECTRON") then
        PhotoElectronHeatingEfficiency_est = output_est
     end if

     ustep = ustep + 1
  end if

  ! END of the first RCMR loop

  ! Send the F10.7 Estimate to all the different processors
  ! AGB Question: why use MPI_SCATTER and not MPI_Bcast?
  if (iProc==0) then
     scattered(:) = output_est
  end if

  call MPI_SCATTER(scattered, 1, mpi_double_precision, scatter, 1, &
       mpi_double_precision, 0, iCommGITM, iError)

  if(RCMROutType == 'F107') then
     f107_est  = scatter    
     f107a_est = f107_est
  else if(RCMROutType == "PHOTOELECTRON") then
     PhotoElectronHeatingEfficiency_est = scatter
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
