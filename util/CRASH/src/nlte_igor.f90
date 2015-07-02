!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! ***  nlte.f90  ***
!
!
!  code will use :
!
! initialisation:
!	call exp_tab8()		! prepare tabulated exponentials
!		call prepCorrUbar()	! or it will called automatically
!	call prep_projE(hnuGr(0:nbGr),nbGR)	! using group. definition of the run
!
!  runtime:
!       call set_kbr
!	call nLTE_EOS(...
!	call correctKnuJnu(...	! not yet implemented
!
! ToDo :
!	- check Auger
!	- use Sum(E/B...) for FB thresholds ?
!	- grey processing 
!                          AT32S * exp(-UM/gama)
!                  1 =  -----------------------  == AT32S * exp(-Um*Te/Tz/gama)
!                           1 + BT72S * UM**3
!	- analytical fit in transparent approximation
!
!-----------
! \ 
module CRASH_M_NLTE
  ! /
  !This module solves the equation for the effective energy.
  !E_{eff} = E_{tot}^{NLTE} - (3/2) k_B x rho x Z^*(T_Z) * (T_e-T_z)
  !where E_{tot} is the input parameter for the "inverse" equation of state
  !(from which the temperatures should be solved) and E_{eff} is the tabulated
  !LTE equation of state with (rho, T_Z) being the input parameters for the 
  !latter.
  use CRASH_M_RADIOM,only : caltz0
  use CRASH_M_projE,only : mxOut
  use CRASH_M_localProperties
  implicit none

  logical :: useLTE=.false.
  logical :: useEElog=.false.
  logical :: useZbrent=.true.
  logical :: dbg=.false.

  !  NOTE :  use either    (useEElog=.F.,useZbrent=.F.)  or  (useEElog=.T.,useZbrent=.F.)
  !
  !  external call :   
  !		call LTE_EOS_inv
  !		call LTE_EOS_dir
  !
  !  user should fill these before calling  LTE_EOS, 
  ! Te 	: eV
  ! Ne	: el/cm3
  ! Natom : atom/cm3
  ! Etot	: internal energy
  ! Ptot	: total pressure
  !
  ! KbR_E, kBr_P :   C_subV, C_subP , including ro if required,
  !  i.e.    Etot=1.5 *(Zbar+1) * kbR_E * T
  !	   Ptot=     (Zbar+1) * kbR_P * T
  ! Efloor : lower bound of total Energy in the table.
  !	   when  Etot.LE.Efloor, LTE will be turned on
  !
  ! \
  !  CAUTION !
  ! /
  ! Ecold,Tcold : lower bound of total Energy in the table for the actual density
  !	during the process of NLTE< it may come that the input value EEeff turns less than Ecold
  !      Tcold will be returned
  !
  !  useLTE : .T. to save CPU (ex. for H at high temperature): forces LTE
  !
  integer,save :: ng_rad=-1
  real,save,dimension(mxout) :: EoBIn
  real,save,dimension(0:mxout) :: hnu_rad
  ! parameter for LTE_EOS_dir  convergence on 'Zbar'
  real,save :: epsD=1d-5
  integer,save :: niterMax=20

contains
  !-------
  subroutine NLTE_EOS(Natom,RO_in, Te_in,Ee_in,Et_in,Pe_in,Pt_in &
       , estim_Zbar,estim_Tz,estim_Te  &
       ,Zbar_out,Tz_out  ,Te_out,Ee_out,Et_out,Pe_out,Pt_out,Cv_out &
       ,RhoDTzDRho)
    !
    !  to account for  Et=Ee+Ei w/ Te=Ti : use  Ee_* or Et_*
    !
    !  direct/inverse EOS is used depending of presence of  E*_in or Te_in
    ! 
    !   CALL SETERAD(Erad,Brad,hnuG,nbG)  has been set before.
    !
    use CRASH_M_localProperties,only : ro,Ni,zSmall
    implicit none
    real,optional,intent(IN) :: Natom,ro_in 	! atomic density in   atom/cm3 
    real,optional,intent(IN) :: Te_in,Ee_in,Et_in,Pe_in,Pt_in 	! one and only one
    real,optional,intent(IN) :: estim_Zbar,estim_Tz,estim_Te
    real,optional,intent(OUT) :: Zbar_out,Tz_out
    real,optional,intent(OUT) :: Te_out,Ee_out,Et_out,Pe_out,Pt_out,Cv_out
    real,optional,intent(OUT) :: RhoDTzDRho
    real :: Ee,Pe,Ne,D 
    real :: Te,Tz,Cv,zBar,eCheck, eDiff
    integer :: nIter
    real,parameter :: x_0=0,x_3o2=1.5d0,x_1=1d0,x_1o2=0.5d0

    ! variables for NLTE EOS_inv

    real,parameter :: tol_EE=1d-4	! convergence param. for zbrentEE (linear function)
    logical :: inverse,direct,oneT,twoT


    if(present(Natom)) then
       !	 if(NIgiven .or. ROgiven) goto 104
       Ni=Natom
       if(present(RO_in)) goto 102
       ro=(Ni/avogadro)*atoMass
    else
       if(.not.present(RO_in))goto 102
       ro=RO_in
       Ni=avogadro*(ro/atoMass)
    end if
    !
    ! check argument list consistency
    ! and check  1T (Te=Ti)  .or. 2T ?
    !
140 twoT=present(Ee_in).or. present(Ee_out) .or. present(Pe_in) .or. present(Pe_out)
    oneT=present(Et_in) .or. present(Pt_in) .or. present(Et_out) .or. present(Pt_out)
    if(twoT)then
       if(oneT) goto 103
       zion=0
       direct=present(Te_in) .or. present(Ee_out)
       inverse=present(Ee_in) .or. present(Pe_in) .or. present(Te_out) 
    else
       if(.not.oneT) goto 103
       zion=1 
       direct=present(Te_in) .or. present(Et_out)
       inverse=present(Et_in) .or. present(Pt_in) .or. present(Te_out)
    end if
    if(inverse.and.direct) goto 101
    if(.not.inverse.and..not.direct) goto 101

    if(ng_rad.lt.0) then
       write(*,*)'-P- setErad has to be called before,'
       write(*,*)'    whenever Erad or Brad are modified'
       stop '-P- setErad has to be called before NLTE_OPAC'
    end if

    !--
    !- direct EOS -
    if(direct) then
       !-
       te=Te_in
       call LTE_EOS_dir(te,Ee,Pe,Zbar,Cv)
       if(useLTE.or.zBar<2*zSmall) then
          !- waiting for using the real direct EOS
          tz=te
          if(present(RhoDTzDRho))RhoDTzDRho = 0.0
          go to 200
       else
          
          d=100.
	  nIter=0
          DIR:	  do while ( abs(d).gt.epsD .and. niter.lt.niterMax)
             Ne=Zbar*Ni
             call calTz0(Te,Ne, Tz, EoBIn(1:ng_rad))

             call LTE_EOS_dir(tz,Ee,Pe,Zbar,Cv)	! ro : in module M_localProperties
             if(zBar<=0.0)EXIT DIR	
             d=(Ne-Ni*zbar)/(Ne+Ni*zbar)
             niter=niter+1
	  end do DIR
          if(nIter==nIterMax)call CON_stop('No convergence in NLTE_EOS_dir')
          call correctEOS(zbar , Te, Tz ,EE=Ee, Pe=Pe, Cv=Cv)
       end if   !if non-LTE
       !--
    else	! if(present(TE_in))
       !--
       !- inverse EOS -
       if(present(Ee_in) ) then
          ee=Ee_in
       else
          ee=Et_in
       end if

       !Caclulate the cold equation of state.
       ! call getEcold(ro,Ecold,Tcold)		! 110827


       if(useLTE) then		! 110827
          !- 
          call LTE_EOS_inv(te,ee,pe,Zbar,Cv)! ro, {Erad,Brad} : in module M_localProperties
          Ne=Zbar*Natom
          tz=te
          if(present(RhoDTzDRho))RhoDTzDRho = 0.0
          goto 200
       else
          !- inverse EOS, non LTE -
          !  should use bracket(Tzdiff,tz1,tz2) + zbrent(tzDif,tz1,tz2)  <<<<<<<
          !   	ro= ; Erad= ; Brad= 
          !   T	e= ; Ne=  : update at each iteration ? use Cv_loc ?

          !- waiting for using the real inverse EOS
          call LTE_EOS_inv(tz,ee,pe,Zbar,Cv)	! zion flag must be dealt with 

          Ne=Zbar*Ni
          if(Zbar.lt.2*zSmall) then
             ! for low Z, LTE is a sensible approximation
             te=tz
             if(present(RhoDTzDRho))RhoDTzDRho = 0.0
             goto 200
          end if
          !
          te=tz
          call calTZ0(Te,Ne, Tz, EoBIn(1:ng_rad))

          if(ee-x_3o2*kbR_E*(zbar+zion)*(Te-Tz).lt.Efloor) then
             ! for low EE, LTE is a sensible approximation
             te=tz
             if(present(RhoDTzDRho))RhoDTzDRho = 0.0
             goto 200
          end if
          d=100.
          eDiff=0.0
	  nIter=0
          INV:	  do while ( (abs(d).gt.epsD.or.abs(eDiff)>tol_EE*eE)&
               .and. niter.lt.niterMax)
             Ne=Zbar*Ni
             call calTz0(Te,Ne, Tz, EoBIn(1:ng_rad))

             call LTE_EOS_dir(tz,eCheck,Pe,Zbar,Cv)! ro : in module M_localProperties
             if(zBar<0)call CON_stop('Negative zBar')
             call correctEOS(zbar , Te, Tz ,EE=eCheck, Pe=Pe, Cv=Cv)	
             d=(Ne-Ni*zbar)/(Ne+Ni*zbar)
             eDiff=eE - eCheck
             !Recalculate Te
             Te = Te + eDiff/Cv
             nIter=nIter+1
	  end do INV
          if(nIter==nIterMax)then
             write(*,*)'Ne=',Ne
             write(*,*)'Ni=',Ni
             write(*,*)'zBar=',zBar
             write(*,*)'Ee=',ee
             write(*,*)'eCheck=',eCheck
             write(*,*)'Cv=',Cv
             write(*,*)'Tz=',Tz
             write(*,*)'Te=',Te
             call CON_stop('No Convergence in EOS_NLTE_Inv')
          end if
       end if	! useLTE
       !--
    end if 	! if(present(TE_in))
    if(present(RhoDTzDRho))then
       Ne=Zbar*Ni
       call calTz0(Te,Ne, Tz, EoBIn(1:ng_rad),RhoDTzDRho)
    end if
200 if(present(zbar_out)) zbar_out=zbar
    if(present(Tz_out)) Tz_out=Tz
    if(present(Te_out)) Te_out=te
    if(present(Ee_out)) Ee_out=ee
    if(present(Et_out)) Et_out=ee
    if(present(Pe_out)) Pe_out=pe
    if(present(Pt_out)) Pt_out=pe
    if(present(Cv_out)) Cv_out=Cv

    return
    !
103 write(*,*)'-P- LTE_EOS require  using either  Ee/Pe .xor. Et/Pt'
    stop '-P- LTE_EOS.103'
102 write(*,*)'-P- LTE_EOS require  Natom .xor. RO_in, or prior use to "set_RO_NI"'
    stop '-P- LTE_EOS.102'
101 write(*,*)'-P- LTE_EOS require  TE_in .xor. EE_in/ET_in (direct EOS .xor. inverse EOS)'
    stop '-P- LTE_EOS.101'
104 write(*,*)'-W- LTE_EOS : use either optional argument (Ni or Ro) or prior call to Set_NI_RO !! '
    !	stop '-P- LTE_EOS.104'
    goto 140
    !
  end subroutine NLTE_EOS
  !-------
  subroutine correctEOS(zbar , Te,Tz &
       ,Ee,Pe,Cv)
    ! Te,Tz, kBro_E,kBro_P,Etot,Ptot : thorugh module variables
    use CRASH_M_localProperties,only : Zsmall
    implicit none
    real,intent(IN) :: zbar,Te,Tz
    real,optional,intent(INOUT) :: Ee
    real,optional,intent(INOUT) :: Pe
    real,optional,intent(INOUT) :: Cv
    real :: zdt,zp
    real,parameter :: x_3o2=1.5d0
    ! kbRO is the proportionnality factor in the code units (generally = Boltzmann cst * density)
    !   so  3/2* kbRo *Te  is the kinetic (translational) energy at Te
    if(1.01*zBar<zSmall)write(*,*)'Te,Tz,zBar=',Te,Tz,zBar
    !if(Te<0.01.or.Te>5e3.or.Tz<0.01.or.Tz>5e3.or.Tz<0.01*Te.or. Te<0.01*Tz.or.zBar<0)then
    !   write(*,*)'Correct EOS:Tz=',Tz,'  Te=',Te,'  zBar=',zBar
    !   call CON_stop('Stop')
    !end if
    !  small Zbar would yield diverging the "EEeff equ." (see correctEOS + EEdiff)
    !no	zp=(zbar+zion)
    zp=(max(zbar,Zsmall) + zion)
    zdt=zp * ( Te-Tz)
    if(present(Ee))then
       if(Ee*1.01<=0)call CON_stop('negative Ee in correct EOS')
       Ee=Ee + x_3o2 * kBr_E * zdt
    end if
    if(present(Pe)) then
       if(Pe*1.01<=0) call CON_stop('negative Pe in correct EOS')
       Pe=Pe +         kBr_P * zdt
    end if
    if(present(Cv)) then
       if(1.01*Cv<=0)call CON_stop('negative Cv in correct EOS')
       Cv=Cv*(Tz/Te)+ x_3o2 * kBr_E * zp*(1-Tz/Te)
    end if
  end subroutine correctEOS
  !==========================

  subroutine setErad(eg_o_bg,hnug,ng)
    !Subroutine re-assigns 
    use CRASH_M_projE,only : prep_projE
    implicit none
    integer,intent(IN) :: ng                 !The dimension of array EOverB
    real,dimension(0:ng),intent(in) :: hnug  !The group boundaries, in eV
    real,dimension(*),intent(in) :: eg_o_bg !Input array of EOverB
    integer :: i

    if(ng.ne.ng_rad.and.ng>0) then
       call prep_projE(hnug,ng)
       ng_rad = ng
       hnu_rad(0:ng_rad)=hnug(0:ng)
    elseif(ng==0)then
       ng_rad = ng
       write(*,*)
       write(*,'(a)')'-W- RADIOM activated with NG=0'
       write(*,*)
    end if
    if(ng>0)then
       EoBIn(1:ng_Rad)=eg_o_bg(1:ng_Rad)
    else
       EoBIn(1)=eg_o_bg(1)
    end if
  end subroutine setErad

 
  ! \
end module CRASH_M_NLTE
! /
!-------



