!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModB0Matvec

  use ModPotentialField, ONLY:  iComm, iProc, nProc, iProcTheta, iProcPhi, &
       nTheta, nPhi, &
       tmpXPhi0_II,tmpXPhipi_II, &
       nThetaLgr,nThetaSml,nPhiLgr,nPhiSml, &
       nProcThetaLgr,nProcThetaSml,nProcPhiLgr,nProcPhiSml, &
       nR, nTheta, nPhi, Radius_I, SinTheta_I, &
       dRadiusNode_I, dTheta_I, dCosTheta_I, dThetaNode_I, dPhiNode_I, &
       Br_II, set_boundary, &
       UseCosTheta, RadiusNode_I, Theta_I, SinThetaNode_I, dCosThetaNode_I, &
       iRTest, iThetaTest, iPhiTest, ThetaNode_I, Phi_I, PhiNode_I

  use ModMPI
  use ModUtilities, ONLY: flush_unit
  use ModIoUnit, ONLY: STDOUT_

  implicit none

contains

  !============================================================================

  subroutine matvec(x_C, y_C, n)

    use ModPotentialField, ONLY: B0_DF, UsePreconditioner, UseHypre, &
         nR, nTheta, d_I, e_I, e1_I, e2_I, f_I, f1_I, f2_I
    use ModLinearSolver, ONLY: Lhepta, Uhepta
    use ModHypre,        ONLY: hypre_preconditioner

    integer, intent(in) :: n
    real, intent(in)    :: x_C(n)
    real, intent(out)   :: y_C(n)
    !--------------------------------------------------------------------------

    ! Calculate y = laplace x in two steps
    call get_gradient(x_C, B0_DF)
    call get_divergence(B0_DF, y_C)

    if(UsePreconditioner)then
       if(UseHypre)then
          call hypre_preconditioner(n, y_C)
       else
          ! Preconditioning: y'= U^{-1}.L^{-1}.y
          call Lhepta(        n, 1, nR, nR*nTheta, y_C, d_I, e_I, e1_I, e2_I)
          call Uhepta(.true., n, 1, nR, nR*nTheta, y_C,      f_I, f1_I, f2_I)
       end if
    end if

  end subroutine matvec

  !============================================================================

  subroutine get_gradient(x_C, Grad_DG)

    real, intent(in):: x_C(nR,nTheta,nPhi)
    real, intent(out):: Grad_DG(3,nR+1,nTheta+1,nPhi+1)

    real, allocatable, save :: x_G(:,:,:)

    integer:: iR, iTheta, iPhi

    ! real:: r, GradExact_D(3)

    !--------------------------------------------------------------------------
    if(.not.allocated(x_G))then
       allocate(x_G(0:nR+1,0:nTheta+1,0:nPhi+1))
       ! Initialize so that corners are all set
       x_G = 0.0
    end if

    call set_boundary(x_C, x_G)
    
    ! This initialization is only for the corners
    Grad_DG = 0.0

    do iPhi = 1, nPhi
       do iTheta = 1, nTheta
          do iR = 1, nR+1
             Grad_DG(1,iR,iTheta,iPhi) = &
                  (x_G(iR,iTheta,iPhi) - x_G(iR-1,iTheta,iPhi)) &
                  / dRadiusNode_I(iR)
          end do
       end do
    end do

    if(UseCosTheta)then
       do iPhi = 1, nPhi
          do iTheta = 1, nTheta+1
             do iR = 1, nR
                Grad_DG(2,iR,iTheta,iPhi) = &
                     (x_G(iR,iTheta,iPhi) - x_G(iR,iTheta-1,iPhi)) &
                     *SinThetaNode_I(iTheta) &
                     / (Radius_I(iR)*dCosThetaNode_I(iTheta))
             end do
          end do
       end do
    else
       do iPhi = 1, nPhi
          do iTheta = 1, nTheta+1
             do iR = 1, nR
                Grad_DG(2,iR,iTheta,iPhi) = &
                     (x_G(iR,iTheta,iPhi) - x_G(iR,iTheta-1,iPhi)) &
                     / (Radius_I(iR)*dThetaNode_I(iTheta))
             end do
          end do
       end do
    end if

    do iPhi = 1, nPhi+1
       do iTheta = 1, nTheta
          do iR = 1, nR
             Grad_DG(3,iR,iTheta,iPhi) = &
                  (x_G(iR,iTheta,iPhi) - x_G(iR,iTheta,iPhi-1)) &
                  / (Radius_I(iR) &
                  *max(1e-10,SinTheta_I(iTheta))*dPhiNode_I(iPhi))
          end do
       end do
    end do

    ! Calculate discretization error for the l=m=1 harmonics
    !iR = iRTest; iPhi = iPhiTest; iTheta = iThetaTest
    !
    !r = Radius_I(iR)
    !GradExact_D  = (/ &
    !     (1+2*rMax**3/RadiusNode_I(iR)**3)/(1+2*rMax**3) &
    !     *sin(Theta_I(iTheta))*cos(Phi_I(iPhi)), &
    !     (r-rMax**3/r**2)/(1+2*rMax**3)/r &
    !     *cos(ThetaNode_I(iTheta))*cos(Phi_I(iPhi)), &
    !     -(r-rMax**3/r**2)/(1+2*rMax**3)/r*sin(PhiNode_I(iPhi)) /)
    !
    !write(*,*) 'magnetogram at test cell=', Br_II(iTheta,iPhi)
    !do iDim = 1, 3
    !   write(*,*) 'Grad, Exact, Error=', &
    !        Grad_DG(iDim,iR,iTheta,iPhi), GradExact_D(iDim), &
    !        Grad_DG(iDim,iR,iTheta,iPhi) - GradExact_D(iDim)
    !end do

  end subroutine get_gradient

  !============================================================================

  subroutine get_divergence(b_DG, DivB_C)

    use ModPotentialField, ONLY: nR, Radius_I, dRadius_I, &
         dPhi_I, SinTheta_I, dCosTheta_I, RadiusNode_I, &
         SinThetaNode_I, nTheta, nPhi

    real, intent(in) :: b_DG(3,nR+1,nTheta+1,nPhi+1)
    real, intent(out):: DivB_C(nR,nTheta,nPhi)

    ! real:: r, DivExact_D(3), Div_D(3)
    integer:: iR, iTheta, iPhi
    !--------------------------------------------------------------------------
    do iPhi = 1, nPhi
       do iTheta = 1, nTheta
          do iR = 1, nR
             DivB_C(iR,iTheta,iPhi) = &
                  ( RadiusNode_I(iR+1)**2*b_DG(1,iR+1,iTheta,iPhi)   &
                  - RadiusNode_I(iR)**2  *b_DG(1,iR  ,iTheta,iPhi) ) &
                  / (Radius_I(iR)**2 *dRadius_I(iR)) &
                  + &
                  ( SinThetaNode_I(iTheta+1)*b_DG(2,iR,iTheta+1,iPhi)   &
                  - SinThetaNode_I(iTheta)  *b_DG(2,iR,iTheta  ,iPhi) ) &
                  / (Radius_I(iR)*dCosTheta_I(iTheta)) &
                  + &
                  ( b_DG(3,iR,iTheta,iPhi+1) &
                  - b_DG(3,iR,iTheta,iPhi) ) &
                  / (Radius_I(iR)*max(1e-10,SinTheta_I(iTheta))*dPhi_I(iPhi))
          end do
       end do
    end do

    ! Calculate discretization error for the l=m=1 harmonics
    !iR = iRTest; iPhi = iPhiTest; iTheta = iThetaTest
    !r = Radius_I(iR)
    !
    !Div_D(1) = ( RadiusNode_I(iR+1)**2*b_DG(1,iR+1,iTheta,iPhi)   &
    !     - RadiusNode_I(iR)**2  *b_DG(1,iR  ,iTheta,iPhi) ) &
    !     / (Radius_I(iR)**2 *dRadius_I(iR))
    !
    !Div_D(2) = ( SinThetaNode_I(iTheta+1)*b_DG(2,iR,iTheta+1,iPhi)   &
    !     - SinThetaNode_I(iTheta)  *b_DG(2,iR,iTheta  ,iPhi) ) &
    !     / (Radius_I(iR)*dCosTheta_I(iTheta))
    !
    !Div_D(3) = ( b_DG(3,iR,iTheta,iPhi+1) - b_DG(3,iR,iTheta,iPhi) ) &
    !     / (Radius_I(iR)*SinTheta_I(iTheta)*dPhi_I(iPhi))
    !
    !DivExact_D = &
    !     (/ 2*SinTheta_I(iTheta), &
    !     (1-SinTheta_I(iTheta)**2)/SinTheta_I(iTheta), &
    !     - 1/SinTheta_I(iTheta) /)
    !
    !DivExact_D = DivExact_D &
    !     *(r-rMax**3/r**2)/(1+2*rMax**3)/r**2*cos(Phi_I(iPhi))
    !
    !do iDim = 1, 3
    !   write(*,*) 'Div_D, Exact, Error=', Div_D(iDim), DivExact_D(iDim), &
    !        Div_D(iDim) - DivExact_D(iDim)
    !end do
    !   
    !write(*,*)'testlaplace=', DivB_C(iR,iTheta,iPhi)
    !write(*,*)'location   =', maxloc(abs(DivB_C))
    !write(*,*)'max laplace=', maxval(abs(DivB_C))
    !write(*,*)'avg laplace=', sum(abs(DivB_C))/(nR*nThetaAll*nPhiAll)
    !
    !stop

  end subroutine get_divergence

end module ModB0Matvec

!==============================================================================

program potential_field

  ! Solve 3D potential field with given Br at inner boundary,
  ! radial field at outer boundary.

  use ModPotentialField
  use ModB0Matvec, ONLY: get_gradient, get_divergence, matvec
  use ModLinearSolver, ONLY: bicgstab, prehepta
  use ModHypre,        ONLY: hypre_initialize, hypre_solver
  use ModPlotFile, ONLY: save_plot_file
  use ModUtilities, ONLY: flush_unit
  use ModIoUnit, ONLY: STDOUT_
  use ModMpi

  implicit none

  integer :: nIter=10000
  real    :: r, DivBMax, DivBMaxAll
  integer :: n, i, iError, iR, iPhi, iTheta
  !--------------------------------------------------------------------------

  call MPI_init(iError)
  call MPI_comm_rank(iComm,iProc,iError)
  call MPI_comm_size(iComm,nProc,iError)

  call read_fdips_param

  if(DoReadMagnetogram .and. iProc == 0) call read_magnetogram

  call MPI_bcast(UseCosTheta, 1, MPI_LOGICAL, 0, iComm, iError)
  call MPI_bcast(nThetaAll, 1, MPI_INTEGER, 0, iComm, iError)
  call MPI_bcast(nPhiAll,   1, MPI_INTEGER, 0, iComm, iError)
  if (.not. allocated(Br_II)) allocate(Br_II(nThetaAll,nPhiAll))

  call MPI_bcast(Br_II, nThetaAll*nPhiAll, MPI_REAL, 0,  iComm, iError)

  if(UseTiming) TimeStart = mpi_wtime()

  call init_potential_field

  if(.not.DoReadMagnetogram)then
     allocate(Br_II(nThetaAll,nPhiAll))
     do iPhi = 1, nPhiAll; do iTheta = 1, nThetaAll; 
        ! magnetogram proportional to the l=m=n harmonics
        n = 1 ! or 2
        Br_II(iTheta,iPhi) = sin(Theta_I(iTheta))**n *cos(n*Phi_I(iPhi))

        ! Exact solution
        do iR = 1, nR
           r = Radius_I(iR)
           Potential_C(iR,iTheta,iPhi) = Br_II(iTheta,iPhi) &
                * (r**n - rMax**(2*n+1)/r**(n+1)) &
                / (n    + (n+1)*rMax**(2*n+1))
        end do
     end do; end do

     write(*,*)'rTest    =',Radius_I(iRTest)
     write(*,*)'PhiTest  =',Phi_I(iPhiTest)
     write(*,*)'ThetaTest=',Theta_I(iThetaTest)
     write(*,*)'BrTest   =',Br_II(iThetaTest,iPhiTest)
     write(*,*)'PotTest  =',Potential_C(iRTest,iThetaTest,iPhiTest)

  end if

  n = nR*nTheta*nPhi
  ! write(*,*) 'iProc, n', iProc, n

  if(UsePreconditioner .or. UseHypre)then

     allocate(d_I(n), e_I(n), f_I(n), e1_I(n), f1_I(n), e2_I(n), f2_I(n))

     i = 0
     do iPhi = 1, nPhi; do iTheta = 1, nTheta; do iR = 1, nR

        i = i + 1
        e_I(i)  = RadiusNode_I(iR)**2 &
             /(Radius_I(iR)**2 * dRadiusNode_I(iR) * dRadius_I(iR))

        f_I(i)  = RadiusNode_I(iR+1)**2 &
             /(Radius_I(iR)**2 * dRadiusNode_I(iR+1) * dRadius_I(iR))

        if(UseCosTheta)then
           e1_I(i) = SinThetaNode_I(iTheta)**2 / &
                (Radius_I(iR)**2 * dCosThetaNode_I(iTheta)  *dCosTheta_I(iTheta))

           f1_I(i) = SinThetaNode_I(iTheta+1)**2 /&
                (Radius_I(iR)**2 * dCosThetaNode_I(iTheta+1)*dCosTheta_I(iTheta))
        else
           e1_I(i) = SinThetaNode_I(iTheta) / &
                (Radius_I(iR)**2 * dThetaNode_I(iTheta)  *dCosTheta_I(iTheta))

           f1_I(i) = SinThetaNode_I(iTheta+1) /&
                (Radius_I(iR)**2 * dThetaNode_I(iTheta+1)*dCosTheta_I(iTheta))
        end if

        !e1_I(i) = 0.0

        !f1_I(i) = 0.0

        e2_I(i) = 1/(Radius_I(iR)**2 * SinTheta_I(iTheta)**2 &
             * dPhiNode_I(iPhi) * dPhi_I(iPhi))

        !e2_I(i) = 0.0

        f2_I(i) = 1/(Radius_I(iR)**2 * SinTheta_I(iTheta)**2 &
             * dPhiNode_I(iPhi+1) * dPhi_I(iPhi))

        !f2_I(i) = 0.0

        d_I(i)  = -(e_I(i) + f_I(i) + e1_I(i) + f1_I(i) + e2_I(i) + f2_I(i))

        if(iR     == 1)      d_I(i)  = d_I(i) + e_I(i) ! inner BC
        if(iR     == 1)      e_I(i)  = 0.0
        if(iR     == nR)     d_I(i)  = d_I(i) - f_I(i) ! outer BC
        if(iR     == nR)     f_I(i)  = 0.0

        if (iProcTheta ==0) then
           if(iTheta == 1)      e1_I(i) = 0.0
        end if
        if (iProcTheta == nProcTheta-1) then
           if(iTheta == nTheta) f1_I(i) = 0.0
        end if

        if(.not.UseHypre)then
           if (iProcPhi == 0) then
              if(iPhi   == 1)   e2_I(i) = 0.0
           end if
           if (iProcPhi == nProcPhi-1) then
              if(iPhi   == nPhi)f2_I(i) = 0.0
           end if
        end if
     end do; end do; end do

     ! A -> LU
     if(.not.UseHypre)call prehepta(n, 1, nR, nR*nTheta, PrecondParam, &
          d_I, e_I, f_I, e1_I, f1_I, e2_I, f2_I)

  end if

  if(UseHypre)call hypre_initialize

  UseBr = .true.
  call matvec(Potential_C, Rhs_C, n)
  Rhs_C = -Rhs_C

  UseBr = .false.

  call flush_unit(STDOUT_)
  call mpi_barrier(iComm, iError)

  if(UseHypre .and. .not.UsePreconditioner)then
     call hypre_solver
  else
     call bicgstab(matvec, Rhs_C, Potential_C, .false., n, &
          Tolerance, 'rel', nIter, iError, DoTest=DoTestMe, iCommIn=iComm)
  end if

  UseBr = .true.

  if(UseTiming) TimeEnd = mpi_wtime()
  if(iProc == 0) &
       write(*,*)'nIter, Tolerance, iError=', nIter, Tolerance, iError

  ! report maximum divb
  call get_gradient(Potential_C, B0_DF)
  call get_divergence(B0_DF, DivB_C)
  DivbMax = maxval(abs(DivB_C))
  if(nProc > 1)then
     call MPI_reduce(DivBMax, DivBMaxAll, 1, MPI_REAL, MPI_MAX, 0, &
          iComm, iError)
     if(iProc==0) DivBMax = DivBMaxAll
  end if
  if(iProc ==0)then
     write(*,*) 'max(abs(divb)) = ', DivBMax
     write(*,*) 'nProcTheta, nProcPhi=', nProcTheta, nProcPhi
  end if
  if(UseTiming) write(*,*) 'running time=', TimeEnd - TimeStart

  if(DoSavePotential)then
     allocate(PlotVar_VC(6,nR,nTheta,nPhi))
     PlotVar_VC = 0.0
     PlotVar_VC(1,:,:,:) = Potential_C
     PlotVar_VC(2,:,:,:) = &
          0.5*(B0_DF(1,1:nR,1:nTheta,1:nPhi) + &
          B0_DF(1,2:nR+1,1:nTheta,1:nPhi))
     PlotVar_VC(3,:,:,:) = &
          0.5*(B0_DF(2,1:nR,1:nTheta,1:nPhi) + &
          B0_DF(2,1:nR,2:nTheta+1,1:nPhi))
     PlotVar_VC(4,:,:,:) = &
          0.5*(B0_DF(3,1:nR,1:nTheta,1:nPhi) + &
          B0_DF(3,1:nR,1:nTheta,2:nPhi+1))
     PlotVar_VC(5,:,:,:) = DivB_C
     PlotVar_VC(6,:,:,:) = Rhs_C

     ! Note the fake processor index to be used by redistribute.pl
     write(NameFile,'(a,2i2.2,a,i3.3,a)') &
          trim(NameFilePotential)//'_np01', nProcTheta, nProcPhi,'_', &
          iProcTheta + iProcPhi*nProcTheta, '.out'

     ! Save divb, potential and RHS for testing purposes
     call save_plot_file(NameFile, TypeFileIn=TypeFilePotential, &
          StringHeaderIn='potential field', &
          NameVarIn='r theta phi pot br btheta bphi divb rhs', &
          Coord1In_I=Radius_I(1:nR), &
          Coord2In_I=Theta_I(1:nTheta), &
          Coord3In_I=Phi_I(1:nPhi), &
          VarIn_VIII=PlotVar_VC)

     deallocate(PlotVar_VC)
  end if

  call save_potential_field

  call MPI_FINALIZE(iError)

end program potential_field
!==============================================================================
subroutine CON_stop(String)

  character(len=*), intent(in):: String

  write(*,*) 'ERROR:', String
  stop

end subroutine CON_stop
!==============================================================================
