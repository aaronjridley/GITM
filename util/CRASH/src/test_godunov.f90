!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=====================Test for the Godunov scheme with the eos==========!!!!
!
! Thermodynamical variables and other notations
!         \rho, Rho - the mass density
!         {\cal E}, E - internal energy of the unit of mass
!         e, i - electron, ion
!        V, vol - volume or volumetric
!        \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
!             derivative at constant volume
!        T_{e,i}, Te,Ti - electron and ion temperature
!        iMaterial - integer variable, a signature of the material:
!        iMaterial=0 - xenon
!        iMaterial=1 - beryllium  
program test_Godunov
  use CRASH_ModAtomicMass
  use CRASH_ModEos
  implicit none
  integer,parameter::nX=1000,nDim=1,nVar=3,iMaterial=0
  real::Cons_VC(nVar,1:nX),Prime_VG(nVar,0:nX+1),Flux_VF(nVar,1:nX+1)
  real::CMax_F(1:nX+1),Gamma_G(0:nX+1),InternalEnergy_G(0:nX+1)
  real::GammaMax
  real,parameter::cLength=5.0e-3    !5 mm
  real::Dt,Time
  real,parameter::DX=cLength/nX
  integer::iX
  real,parameter::RhoInit=100.0 !Approximately 30 times the normal density
  real,parameter::UInit = 1.50*RhoInit*(8.31e+3/122.2)*3.0e3 !3000 K 
  real,parameter::UPiston = 3.0e+4 !30 km/s
  real,parameter::CFL=0.8
  real,parameter::TimeOut=5.0e-7 !500 ns
  !^CFG COPYRIGHT UM
  !=====================Equation Of State (EOS)===========================!!!!
  !
  ! Thermodynamical variables and other notations
  !         \rho, Rho - the mass density
  !         {\cal E}, E - internal energy of the unit of mass
  !         e, i - electron, ion
  !        V, vol - volume or volumetric
  !        \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
  !             derivative at constant volume
  !        T_{e,i}, Te,Ti - electron and ion temperature
  !        iMaterial - integer variable, a signature of the material:
  !        iMaterial=0 - xenon
  !        iMaterial=1 - beryllium         

  UsePreviousTe = .true.

  Dt=0.0; Time =0.0 
  !Initial condition, in the frame of reference comoving with a piston
  do iX=1,nX
     Cons_VC(1,iX)=RhoInit
     Cons_VC(2,iX)=-UPiston*RhoInit
     Cons_VC(3,iX)= UInit +0.50*RhoInit*UPiston**2
  end do

  open(24,file='Godunov_scheme_for_Xe',status='replace')
  write(24,'(a)') 'X [mm] Rho [100 kg/m^3]   u [10,000 m/s]  P [MBar]  Gamma'
  do 
     write(*,*)'iStep,StartTime=',Time
     !Get primitives:
     do iX=1,nX
        Prime_VG(1,iX)=Cons_VC(1,iX)
        Prime_VG(2,iX)=Cons_VC(2,iX)/Prime_VG(1,iX)
        InternalEnergy_G(iX)=Cons_VC(3,iX)-0.50*Prime_VG(1,iX)*Prime_VG(2,iX)**2
        call eos(iMaterial  ,               & !Input: sort of material
             Prime_VG(1,iX),                & !Input mass density,[kg/m^3] 
             ETotalIn=InternalEnergy_G(iX), & !Input energy density[J/m^3]
             PTotalOut=Prime_VG(3,iX),      & !Output,pressure [Pa]
             GammaOut=Gamma_G(iX)           ) !Output,polytropic index  
     end do
     if(Time>TimeOut)then
        do iX=1,nX
           write(24,'(5(E14.6,2X))')DX*(iX-0.5)/1.0e-3,&
                Prime_VG(1,iX)/RhoInit,&
                Prime_VG(2,iX)/1.0e4,&
                Prime_VG(3,iX)/1.0e11,&
                Gamma_G(iX)
        end do
        exit
     end if
     !Fix ends
     Prime_VG(:,0)   =Prime_VG(:, 1); Gamma_G(0)   =Gamma_G(1)
     Prime_VG(:,nX+1)=Prime_VG(:,nX); Gamma_G(nX+1)=Gamma_G(nX)

     InternalEnergy_G(0)   =InternalEnergy_G(1)
     InternalEnergy_G( nX+1)   =InternalEnergy_G(nX) 
     !Reflection at the left end
     Prime_VG(2,0) = - Prime_VG(2,0)


     !Get fluxes
     do iX=1,nX+1
        GammaMax=5.0/3.0
        call get_godunov_flux(1,(/1.0/),3,Prime_VG(:,iX-1),Prime_VG(:,iX),&
             Flux_VF(:,iX),CMax_F(iX),GammaMax,GammaMax,&
             (InternalEnergy_G(iX-1)-Prime_VG(3,iX-1)/(GammaMax-1.0))/Prime_VG(1,iX-1),&
             (InternalEnergy_G(iX)-Prime_VG(3,iX)/(GammaMax-1.0))/Prime_VG(1,iX))
     end do
     Cons_VC(:,1:nX)=Cons_VC(:,1:nX)+&
          (Dt/DX)*(Flux_VF(:,1:nX)-Flux_VF(:,2:nX+1))
     Dt=min(CFL*DX/maxval(CMax_F(1:nX+1)),1.001*(TimeOut-Time))

     Time = Time + Dt
  end do
  close(24)
contains
  !=====================================================================!
  !Calculates the Godunov flux with the exact Riemann solver. The
  !interface is the most generic (nDim may be 1,2,3), the state vector
  !of nVar .ge. 2+nDim may inslude an arbitrary number of extra state 
  !variable, for all of the the continuity equation being implied. At 
  !the same time, the interface is the least convenient: it requires to
  !pass nDim and nVar and the unity vector normal to the face which is 
  !(/1,0,0/), (/0,1,0/), (0,0,1/) for 3D cartesian grid while calculate 
  !the fluxes along x,y,z axis respectively. 
  subroutine get_godunov_flux(&
       nDim,Normal_D, & !Unity vector normal to face 
       nVar,LeftState_V,RightState_V,&!Input states
       Flux_V, CMax, &  !Output flux and max perturbation speed
       GammaL,GammaR,&  !Optional Gamma, if needed
       Energy0L,Energy0R, & !Optional Energy0, if needed
       DoTest)
    use ModExactRS
    
    integer,intent(in)::nDim,nVar
    real,intent(in),dimension(nDim)::Normal_D !Unity vector normal to the face
    real,intent(in),dimension(nVar)::LeftState_V,RightState_V
    real,intent(out),dimension(nVar)::Flux_V
    real,intent(out)::CMax
    real,intent(in),optional::GammaL,GammaR !Effective gamma, left and right
    real,intent(in),optional::Energy0L,Energy0R
    logical,optional::DoTest
    real::Rho, Un, P, StateStar_V(nVar)
    real::RhoSide,UnSide,GammaSideM1Inv,E0=0.0
    integer::iVar
    real,parameter::GammaHere=5.0/3.0
    !----------------------------------------
    !Take scalars

    RhoL=LeftState_V(1)   ;RhoR=RightState_V(1)
    PL  =LeftState_V(nVar);PR  =RightState_V(nVar)
    UnL =sum( LeftState_V(2:1+nDim) * Normal_D)
    UnR =sum(RightState_V(2:1+nDim) * Normal_D)

    !Take the parameters at the Contact Discontinuity (CD)

    call pu_star(GammaL,GammaR)

    if(UnStar>0.0)then
       !The CD is to the right from the face
       !The Left gas passes through the face
       RhoSide = RhoL ; UnSide = UnL
       StateStar_V=LeftState_V
       if(present(GammaL))then
          GammaSideM1Inv=1.0/(GammaL-1.0)
       else
          GammaSideM1Inv=1.0/(GammaHere-1.0)
       end if
       if(present(Energy0L))then
          E0=Energy0L
       else
          E0=0.0
       end if
    else
       !The CD is to the left from the face
       !The Right gas passes through the face
       RhoSide = RhoR ; UnSide = UnR
       StateStar_V=RightState_V
       if(present(GammaR))then
          GammaSideM1Inv=1.0/(GammaR-1.0)
       else
          GammaSideM1Inv=1.0/(GammaHere-1.0)
       end if
       if(present(Energy0R))then
          E0=Energy0R
       else
          E0=0.0
       end if
    end if

    !Take the parameters at the face

    call sample(0.0,Rho, Un, P, GammaL,GammaR)
    StateStar_V(1)=Rho
    StateStar_V(2:1+nDim) = StateStar_V(2:1+nDim)+(Un-UnSide)*Normal_D
    StateStar_V(nVar) = P
    do iVar=2+nDim,nVar-1
       StateStar_V(iVar) = StateStar_V(iVar)*(Rho/RhoSide)
    end do

    !Calculate flux (1) take conservative variable

    StateStar_V(nVar) = (1.0 +  GammaSideM1Inv) *  StateStar_V(nVar) &
         + 0.5 * Rho *sum(StateStar_V(2:1+nDim)**2) + Rho*E0
    StateStar_V(2:1+nDim) = StateStar_V(2:1+nDim) * Rho

    ! (2) take advective part of the flux

    Flux_V = StateStar_V * Un 

    ! (3) add the pressure tensor

    Flux_V(2:1+nDim) = Flux_V(2:1+nDim) + P * Normal_D

    CMax = max( WR, -WL)
    if(present(DoTest))then
       if(DoTest)then
          write(*,'(a)')'Test for get_godunov_flux:'
          write(*,*)'Left state:', LeftState_V
          write(*,*)'Right state:',RightState_V
          write(*,*)'Star state:',StateStar_V
          write(*,*)'PStar,UStar=',PStar,UnStar
          write(*,*)'Flux:',Flux_V
          write(*,*)'CMax',CMax
          if(present(GammaL))write(*,*)'GammaL,GammaR=',GammaL,GammaR
       end if
    end if
  end subroutine get_godunov_flux
  !==================================================================!
end program test_Godunov
!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test


