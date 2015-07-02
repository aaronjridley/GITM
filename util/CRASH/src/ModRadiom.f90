!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!------- radiom.f90
! \
module CRASH_M_RADIOM
  ! /
  ! The module solves the main equations of the non-LTE
  ! model
  !
  use CRASH_M_projE
  !Projects the radiation energy densities in user-defined groups 
  !onto a refined grid  
  implicit none
  PRIVATE !Except


  !  Flags :
  logical,save :: setXubar=.true.
  !  Options :
  logical,save :: jhg=.false.,sumMode=.false.,projMode=.true. &
       ,doNew=.true.
  real,save :: Auger=1,Auger_1oB=20d0,Auger_B=0,Gaunt=0.2
  ! work arrays for "corrUbar" :
  integer,parameter :: mxN=1+20*5
  real :: QJ_1,dQJ ,corrUs(mxN),QJs(mxN) 

  real,parameter :: aSaha=6.02e21, b_CovR=1.34e13*0.2 ,one=1.0
  !\
  ! For Indirect Equation of state either ERad or ERad/B(Te) may be
  ! used as inputs.
  !/
  logical:: UseERadInput = .false.

  public:: caltz0,prep_projE,PrepCorrUBar,printversion, UseERadInput

  !-------
contains
   !-------
  subroutine printVersion
    implicit none
    character(*),parameter :: author='M.Busquet' &
         ,wher='ARTEP,inc' &
         ,version='CRASH-2012.006.a' &
         ,library='RADIOM' &
         ,dateModified='2012-06-19 10:58' &
         ,compiled='2012.07.31 08:59:15' &
         ,collab='M.Klapisch & I.Sokolov' &
         ,credits='D.Fyfe & J.H.Gardner'
    !--------------
    write(*,'(a)')'...................................'
    write(*,'(a)')'library=',library
    write(*,'(a)')'version=',version
    write(*,'(a)')'modified on ',dateModified
    write(*,'(a)')'compiled on ',compiled
    write(*,'(a)')'author=',author
    write(*,'(a)')' at ',wher
    write(*,'(a)')' in collaboration with : ',collab
    write(*,'(a)')' credits also to : ',credits
    write(*,'(a)')'...................................'
  end subroutine printVersion
  
 !-------
  subroutine calTz0(Te,Ne, Tz, ERad_I, RhoDTzDRho)
    !-
    !-     same as "calTZ" but w/o hnug(0:ng) & ubar1
    !
    !
    use CRASH_M_projE,only : nbIn
    use ModCOnst, ONLY: cEvToK
    use CRASH_ModMultiGroup,ONLY:get_planck_g_from_temperature
    real,intent(IN) :: Te,ne
    real,dimension(:),intent(IN) :: ERad_I
    real,intent(OUT) :: Tz
    real,intent(out), optional::RhoDTzDRho
    real,parameter :: roVar=0.015d0,roVarP=1+roVar
    real:: lgVar=0.0, QJ1, Ne1, ubar1, Tz1
    real :: at32s,ubar,QJ
    !Array of EOverB
    real:: EoB(nbIn), TeK=0, PlanckSi=0
    integer:: iBin
    if(nbIn.le.0) then
       write(*,*)'-P- prep_projE not done before "caltz0"'
       call CON_stop('-P- calTZ: prep_projE not done')
    end if
   
    EoB = ERad_I
    if(UseERadInput)then
       TeK = Te*cEvToK
       do iBin=1,nBin
          call get_planck_g_from_temperature(iBin, TeK, EgSI=PlanckSi)
          if(0.2*EoB(iBin)>0)then
             EoB(iBin)=EoB(iBin)/max(0.2*EoB(iBin),PlanckSi)
          else
             EoB(iBin)=0
          end if
       end do
    end if
    !!!$
    EoB = 0.0
    !!!$
    at32s=aSaha*te*sqrt(te)/Ne
    tz=te
    if(present(RhoDTzDRho))RhoDTzDRho = 0.0
    if(at32s.le.one) return
    ! this is already done in "setErad", should pass by arg ?
    call xubar0(te,ne,EoB ,ubar)	! no more a function
    QJ=log(at32s)
    tz=te*ubar/QJ
    !if(tz<0.9*Te.or.tz>1.1*Te)then
    !   write(*,*)'Te=',Te,'  Tz=',Tz, 'Ne=',Ne
    !   call CON_stop('Te/=Tz')
    !end if
    if(present(RhoDTzDRho))then
       if(lgVar==0.0) lgVar=log(roVarP)			
       ne1=ne*roVarP				
       call xubar0(te,ne1,EoB,ubar1)		
       QJ1=QJ-lgVar	! =log(at32s)	
       tz1=te*ubar1/QJ1				
       RhoDTzDRho = (tz1-tz)/roVar		
    end if
  end subroutine calTz0
  !===================

  subroutine xubar0(Te_in,Ne_in ,EoB,ubar)
    ! 
    ! same as "xubar" w/o check of Eground
    !
    use CRASH_M_projE
    use CRASH_M_expTab
    real :: ubar
    integer :: ngr
    real,intent(IN) :: Te_in,Ne_in
    real,intent(IN) :: EoB(nbIn)

    real,parameter :: A_LTE_limit=7.38905609	! =exp(2)
    real,parameter :: smallR=1d-200,two=2

    integer :: nfgp,nOut,ig
    real :: ts,at32s,QJ,betapm ,bu3
    real :: du,duu,s,rm,frad
    logical :: gotoLTE
    real :: u,ey		! for 'exp_tab'  u -> ey=exp(-u), ex_u=1-ey

    logical,save :: lastWasPrep=.false.
    ngr=nbIn
    nfgp=ngr+1
    ts=sqrt(Te_in)/Ne_in
    at32s= aSaha*Te_in*ts
    if(at32s.le.A_LTE_limit) then
       ubar=log(at32s)	! thus Tz/Te will be  1
       return
    end if
  
    QJ=log(at32s)
    betapm= (b_CovR/Gaunt)*ts*Te_in**3        
    !	     --------
    call projSP_E(Te_in,EoB,nOut,gotoLTE)	
    !	     --------
    if(gotoLTE) then
       ubar=QJ ! =log(at32s)		! thus Tz/Te will be  1
       return
    end if
    ! perform the "non overflow" alogrithm, u=including  jhg correction
    ubar=0
    s=one
    if(nOut.lt.2) call CON_stop('-E- xubar0 : nbOut<2')
    duu=Uout(2)-Uout(1)
    do ig=2,nOut+1
       u  = (Uout(ig)+Uout(ig-1))/two
       du = (Uout(ig)-Uout(ig-1))
       bu3=betapm*u**3
       include 'exp_tab.h'
       frad=at32s*ey 					&
            *(Auger+Auger_1oB/betapm 			&
            +Auger_B*betapm 				&
            +bu3*SPout(ig-1)) 				&
            /(Auger +Auger_1oB/betapm 			&
            +Auger_B*betapm				&
            +bu3					)
       rm=frad/s*(du/duu)
       rm=max(rm,smallR)
       duu=du
       ubar=(ubar+u*rm)/(one+rm)
       s=one+one/rm		! sumMode
    end do
    call CorrUbar(ubar,QJ)
  end subroutine xubar0
  !=======================

  subroutine CorrUbar(ubar,QJ)

    real, intent(IN) :: QJ
    real :: ubar,d
    real, parameter :: QJ_limit=2.1
    integer :: iq
    !	
    if(setXubar) return
    ! 
    iq=int((QJ-QJ_1)/dQJ) +1
    if(iq <= 0 .or. QJ <= QJ_limit) then
       ubar=QJ
    elseif(iq.ge.mxN) then
       ubar=ubar*corrUs(mxN)
    else
       d=(QJ-QJs(iq))/dQJ
       d=(one-d)*corrUs(iq)+d*corrUs(iq+1)
       ubar=ubar*d
    end if
    ! 
  end subroutine CorrUbar
  !-------
  subroutine prepCorrUbar
    !
    !  tabulate correction to UBAR (see corrUbar)
    !
    use CRASH_M_projE,only : Umin,Umax	,nbIn,Efirst,Elast
    integer,parameter :: mxu=1+10*5 
    real :: QJ_2 !log(6e4)
    real :: ubar,ne,te,at32s,qj,r,eg1 
    real,dimension(0:mxN) :: Ugs
    real,dimension(0:nBIn) :: EovB
    integer :: iq
    ! 
    EovB = 1.0    ! 0 & 1 give same result
    te=0.1
    r=(Umax/Umin)**(one/mxu)
    eg1=Umin*te
    do iq=1,mxN
       Ugs(iq-1)=eg1
       eg1=eg1*r
    end do
    Ugs(mxN)=eg1
    setXubar=.true.     !Disable the call for CorrUBar
    QJ_1=2.
    QJ_2=log(6.d4)
    dQJ=(QJ_2-QJ_1)/(mxN-1)
    qj=QJ_1
    do iq=1,mxN
       at32s=exp(qj)
       ne=te*sqrt(te)*(aSaha/at32s)
       call xubar0(Te,Ne,EovB(0:nbIn),ubar)
       corrUs(iq)=qj/ubar
       QJs(iq)=qj
       qj=qj+dQJ
    end do
    ! 
    setXubar=.false.	

  end subroutine prepCorrUbar
  !================
  ! \
end module CRASH_M_RADIOM
 ! /
