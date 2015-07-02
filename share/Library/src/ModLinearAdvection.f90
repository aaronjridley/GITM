!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModLinearAdvection

  implicit none
  real:: BetaLim=2.0
  public :: test_linear_advection
contains
  !===========advance_lin_advection==========================================!
  !DESCRIPTION: the procedure integrates the linear advection equation,      !
  !using a one-stage second order scheme                                     !
  !==========================================================================!
  subroutine advance_lin_advection_plus(&
       CFLIn_I,         &
       nX,            &
       nGCLeft,       &
       nGCRight,      &
       FInOut_I,      &
       BetaIn, UseConservativeBC, IsNegativeEnergy )      


    !--------------------------------------------------------------------------!
    !                                                                          !
    integer,intent(in):: nX           !Number of meshes                        !
    real,intent(in),dimension(nX)::CFLIn_I      !Time step                     !
    integer,intent(in):: nGCLeft      !The solution in the ghost cells is not  !
    integer,intent(in):: nGCRight     !advanced in time but used as the        ! 
    !                                 boundary condition, if any               !
    !                                                                          !
    !In: sol. to be advanced; Out: advanced sol.                               !
    real,dimension(1-nGCLeft : nX+nGCRight),intent(inout)::  FInOut_I          !
    real, optional::BetaIn

    !If this parameter IS PRESENT, the flux via the boundaries is nullified
    logical, optional, intent(in):: UseConservativeBC
    logical, optional, intent(out):: IsNegativeEnergy

    !--------------------------------------------------------------------------!
    integer:: iX,iStep,nStep
    real,dimension(1-max(nGCLeft,2):nX+max(nGCRight,2))::F_I
    real,dimension(0:nX+1):: FSemiintUp_I,FSemiintDown_I
    real,dimension(0:nX+1)::CFL_I
    character(LEN=*),parameter :: NameSub = 'advance_lin_advection_plus'
    !--------------------------------------------------------------------------!
    if(present(BetaIn))then
       BetaLim = BetaIn
    else
       BetaLim = 2.0
    end if

    if(present(IsNegativeEnergy))IsNegativeEnergy = .false.

    nStep=1+int(maxval(CFLIn_I)); CFL_I(1:nX) = CFLIn_I/real(nStep)
    CFL_I(0) = CFL_I(1); CFL_I(nX+1) = CFL_I(nX)
    F_I(1-nGCLeft : nX+nGCRight) = FInOut_I(1-nGCLeft : nX+nGCRight)

    !Check for positivity
    if(any(F_I(1:nX)<0.0))then
       write(*,*)'Before advection F_I < 0 in '//NameSub
       if(present(IsNegativeEnergy))then
          IsNegativeEnergy = .true.
          return
       else
          write(*,*)'F_I:', F_I
          call CON_stop('Error in '//NameSub )
       end if
    end if

    !One stage second order upwind scheme

    do iStep=1,nStep

       !Boundary condition at the left boundary
       if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft )
       !Boundary condition at the right boundary
       if(nGCRight<2)F_I(nX+1+nGCRight:nX+2     ) = F_I(nX+nGCRight)

       do iX=0,nX

          ! f_(i+1/2):
          FSemiintUp_I(iX) = F_I(iX)+&
               0.50*(1.0 - CFL_I(iX))*df_lim(F_I(iX-1:iX+1))
       end do
       ! f_(i-1/2): 
       FSemiintDown_I(1:nX) = FSemiintUp_I(0:nX-1)

       !Ensure the conservation of sum(F_I(1:nX), if needed
       if(present(UseConservativeBC))then
          FSemiintDown_I(1) = 0.0
          FSemiintUp_I( nX) = 0.0
       end if

       !\
       ! Update the solution from f^(n) to f^(n+1):
       !/
       F_I(1:nX) = F_I(1:nX)+CFL_I(1:nX)*&
            (FSemiintDown_I(1:nX)-FSemiintUp_I(1:nX))
    end do

    FInOut_I(1-nGCLeft:nX+nGCRight) = F_I(1-nGCLeft:nX+nGCRight)
    if(any(FInOut_I(1:nX)<0.0))then
       write(*,*)'After advection F_I <0 in '//NameSub
       write(*,*)'F_I:',F_I
       if(present(IsNegativeEnergy))then
          IsNegativeEnergy = .true.
          return
       else
          call CON_stop('Error in '//NameSub )
       end if
    end if
  end subroutine advance_lin_advection_plus

  !===========advance_lin_advection==========================================!
  !DESCRIPTION: the procedure integrates the log-advection equation, in the  !
  !conservative or non-conservative formulation, at a logarithmic grid, using!
  !a one-state second order scheme                                           !
  !==========================================================================!
  subroutine advance_lin_advection_minus(&
       CFLIn_I,         &
       nX,            &
       nGCLeft,       &
       nGCRight,      &
       FInOut_I,      &
       BetaIn, UseConservativeBC, IsNegativeEnergy)      
    !--------------------------------------------------------------------------!
    !
    integer,intent(in):: nX           !Number of meshes                        !
    real,intent(in),dimension(nX)::CFLIn_I      !Time step   
    integer,intent(in):: nGCLeft      !The solution in the ghost cells is not  !
    integer,intent(in):: nGCRight     !advanced in time but used as the        ! 
    !boundary condition, if any              !
    real,dimension(1-nGCLeft:nX+nGCRight),intent(inout):: &
         FInOut_I                     !In: sol. to be advanced; Out: advanced  !
    real, optional :: BetaIn

    !If this parameter IS PRESENT, the flux via the boundaries is nullified
    logical, optional, intent(in):: UseConservativeBC

    logical, optional, intent(out):: IsNegativeEnergy 

    !sol. index 0 is to set boundary         !
    !condition at the injection energy.      !
    !--------------------------------------------------------------------------!
    integer:: iX,iStep,nStep
    real,dimension(1-max(nGCLeft,2):nX+max(nGCRight,2))::F_I
    real,dimension(0:nX+1):: FSemiintUp_I,FSemiintDown_I
    character(LEN=*),parameter :: NameSub = 'advance_lin_advection_minus'
    !--------------------------------------------------------------------------!


    real,dimension(0:nX+1)::CFL_I
    !--------------------------------------------------------------------------!
    if(present(BetaIn))then
       BetaLim = BetaIn
    else
       BetaLim = 2.0
    end if
    if(present(IsNegativeEnergy))IsNegativeEnergy = .false.

    nStep=1+int(maxval(CFLIn_I)); CFL_I(1:nX) = CFLIn_I/real(nStep)
    CFL_I(0) = CFL_I(1); CFL_I(nX+1) = CFL_I(nX)


    F_I(1-nGCLeft : nX+nGCRight) = FInOut_I(1-nGCLeft:nX+nGCRight)

    !Check for positivity
    if(any(F_I(1:nX)<0.0))then
       write(*,*)'Before advection F_I <0 in '//NameSub

       if(present(IsNegativeEnergy))then
          IsNegativeEnergy = .true.
          return
       else
          write(*,*) 'F_I:', F_I
          call CON_stop('Error in '//NameSub )
       end if
    end if

    !One stage second order upwind scheme

    do iStep=1,nStep
       !Boundary condition at the left boundary
       if(nGCLeft  < 2) F_I(            -1:0-nGCLeft) = F_I( 1 - nGCLeft )
       !Boundary condition at the right boundary
       if(nGCRight < 2) F_I(nX +1 + nGCRight:nX+2   ) = F_I(nX + nGCRight)
       do iX=1,nX+1
          ! f_(i-1/2):
          FSemiintDown_I(iX) = F_I(iX)&
               -0.50 * (1.0 - CFL_I(iX)) *  df_lim(F_I(iX-1:iX+1))
       end do
       ! f_(i+1/2):
       FSemiintUp_I(1:nX) = FSemiintDown_I(2:nX+1)

       !Ensure the conservation of sum(F_I(1:nX), if needed
       if(present(UseConservativeBC))then
          FSemiintDown_I(1) = 0.0
          FSemiintUp_I( nX) = 0.0
       end if
       !\
       ! Update the solution from f^(n) to f^(n+1):
       !/
       F_I(1:nX) = F_I(1:nX) - CFL_I(1:nX) * (FSemiintDown_I(1:nX) - FSemiintUp_I(1:nX))
    end do


    FInOut_I(1:nX)=F_I(1:nX)

    if(any(FInOut_I(1:nX)<0.0))then
       write(*,*)'After advection F_I <0 in '//NameSub
       write(*,*)'F_I:', F_I
       if(present(IsNegativeEnergy))then
          IsNegativeEnergy = .true.
          return
       else
          call CON_stop('Error in '//NameSub )
       end if
    end if
    !------------------------------------ DONE --------------------------------!
  end subroutine advance_lin_advection_minus
  !============================================================================!
  subroutine test_linear_advection
    ! Added Nov. 2009 by R. Oran

    use ModIoUnit, ONLY: UNITTMP_

    integer,parameter                          :: nGCLeft=1, nGCRight=1     
    integer,parameter                          :: nCell= 80,nStep=40
    integer                                    :: iCell, iStep
    real,dimension(nCell),parameter            :: CFL = 0.9
    real,dimension(1-nGCLeft : nCell+nGCRight) :: Fplus_I = 0.0,&
         Fminus_I=0.0 ! solution vectors
    logical                                    :: IsNegativeEnergy
    ! character(len=4)                         :: NameStage
    character(len=40)                          :: FileNameTec
    ! ------------------------------------------------------
    ! Initial condition - create a uniform spectrum

    Fplus_I  = 1.0 ! moving right
    Fminus_I = 1.0 ! moving left
    ! Start looping over time steps
    do iStep =1,nStep
       call advance_lin_advection_plus(CFL,nCell,nGCLeft,nGCRight,Fplus_I,&
            2.0,.true.,IsNegativeEnergy)
       call advance_lin_advection_minus(CFL,nCell,nGCLeft,nGCRight,Fminus_I,& 
            2.0,.true.,IsNegativeEnergy)

    end do

    ! write sln to file
    ! Deafault output file name
    FileNameTec = 'test_linear_advection.tmp'

    ! use this filename when an output file for each iteration is desired
    ! write(NameStage,'(i4.4)') iStep
    ! FileNameTec = 'LinAdvOut/Linear_advection_n_'//trim(NameStage)//'.dat'

    open(UNITTMP_,file=FileNameTec,form='formatted', access='sequential',&
         status='replace')
    write(UNITTMP_, '(a)') 'Title: Test Linear Advection'
    write(UNITTMP_, '(a)') 'Variables = "Cell","I+","I-"'
    write(UNITTMP_, '(a,i3.3,a,i5.5,a)') 'Zone I= ',nCell,' F=point'
    do iCell =1,nCell
       write(UNITTMP_,'(i3.3,e14.6,e14.6)') iCell,Fplus_I(iCell), Fminus_I(iCell)
    end do
    close(UNITTMP_)

  end subroutine test_linear_advection
  !====================SUPERBEE LIMITER =======================================!
  real function df_lim(F_I)
    real,dimension(0:2),intent(in)::F_I
    !--------------------------------------------------------------------------!
    integer,parameter:: i=1
    real:: dF1,dF2
    !--------------------------------------------------------------------------!
    dF1 = F_I(i+1) - F_I(i)
    dF2 = F_I(i) - F_I(i-1)

    !df_lim=0 if dF1*dF2<0, sign(dF1) otherwise:

    df_lim = sign(0.50,dF1) + sign(0.50,dF2)
    dF1 = abs(dF1)
    dF2 = abs(dF2)
    df_lim = df_lim * min(max(dF1,dF2), BetaLim*dF1, BetaLim*dF2)
    !---------------------------------- DONE ----------------------------------!
  end function df_lim
  !============================================================================!
end module ModLinearAdvection
