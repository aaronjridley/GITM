!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-  testNLTE.f90 -  in progress
!	testing NLTE_EOS, using zBrent for  te,tz finding given ro,Etot & EOS=ZTF_Perrot
!      compile with :  radiom.f90, eos_material.f90

!- testZTF.f90 - checked
!      testing  eos_material (ZTF_EOS_dir, ZTF_EOS_inv & ZTF*s
!	compule with : eos_material.f90, testZTF.f90

!- testradiom.f90 - checked
!	testing RADIOM routines : CALTZ, CALTE, 
!  f90 -q -O -z 2 -en -Z899 -Z124 -Z938 -Z1553 -o tn -g testNLTE.f90 radiom.o
!	

program testNLTE

  use CRASH_M_EOS,only : set_kBr
  use CRASH_M_expTab,only: exp_tab8
  use CRASH_M_localProperties,only : kBr_E,kBr_P ,zion,atoMass &
       ,ERGperEV,DYNEperEV,avogadro
  use M_RADIOM, only : prep_projE,prepCorrUbar &
       ,printVersion  !xx ,calTz0
  use CRASH_M_NLTE,only : nlte_EOS ,useLTE ,ng_rad,EoB &!, set_RO_NI &
       ,useEElog,useZbrent


  implicit none
  real :: Ni_l,Te_l,ro_l,Ee_l,Tz_l &
       ,zbar_l, pe_l, Cv_l  ,tz_new 	! ,Ne
  real :: eeff,d,dl
  real,external :: EEdiff,EEdiff_ln

  real ::  ro_in , ee_in, te_in  &
       ,LTE_ee_outD ,LTE_zb_outD ,LTE_pe_outD ,LTE_cv_outD &
       ,LTE_te_outI ,LTE_zb_outI ,LTE_pe_outI ,LTE_cv_outI &
       ,NLTE_ee_outD,NLTE_zb_outD,NLTE_tz_outD,NLTE_cv_outD &
       ,NLTE_te_outI,NLTE_zb_outI,NLTE_tz_outI,NLTE_cv_outI 

  integer :: ir,it
  real,dimension(0:45),parameter :: gr45=(/ &
       1., 10., 48.75, 50.5, 58.5, 67.5 ,68.5, 79.25, 91., 98.5, 105. &
       , 127.5 ,150., 167.33, 184.67, 202., 219.33, 236.67 ,254., 267., 287., 300. &
       , 308., 309.5 ,320., 335., 345., 361., 364., 369. ,376., 378., 403.5, 435.5 &
       , 436.5, 446.2 ,474., 480., 540., 600., 757., 950. ,1250., 2500. &
       , 12000., 100000. &
       /)
  real,dimension(301) :: ones =1.d0

  !
  call printVersion()
  !\
  ! Tabulates the exponential function.
  !/
  call exp_tab8()

  !\
  ! Coefficients for transforming from the user defined grid to
  ! the refined logrithmic-uniform internal fixed grid
  !/ 
  call prep_projE(gr45,45)
  ng_rad=45

  !\
  ! Initialize and calculate some internal arrays
  !/
  call prepCorrUbar()
  !
  !=============  Set material + Method  (for the test)  ==================
  call setMethod()


  ! \
  !  Write 'fort.45' file for checking  identity of  EOS_direct x EOS_inverse :
  ! /
  write(45,'(10a)') '  ro_in  ee_in te_in ' &
       ,' LTE_ee_outD LTE_zb_outD LTE_pe_outD LTE_cv_outD' &
       ,' LTE_te_outI LTE_zb_outI LTE_pe_outI LTE_cv_outI' &
       ,' NLTE_ee_outD NLTE_zb_outD NLTE_tz_outD NLTE_cv_outD' &
       ,' NLTE_te_outI NLTE_zb_outI NLTE_tz_outI NLTE_cv_outI' 

  EoB(1:ng_rad)=0

  goto 1000		! skip this branch, except to visualize testFunction

  !==================================================
  !  vizualize the  values of the test function vs EEeff for one (ro,EE)
  ! 
1001 write(46,*) '  ro_l  Ni_l  Eeff  ee_l   Te_l  Tz_l  diff  diff_ln'
  RO_IN=3e-7
  ro_l=RO_IN
  Ni_l=ro_l*avogadro/atoMass
  call set_kBr(Natom=Ni_l)
  ! call set_RO_NI(atoMass=atoMass,RO_in=ro_l)
  ee_l=1e-15	! until 1e-8
  Eeff=1.4188E-08
  do ir=-12,-8
     Eeff=1.4188 *10d0**ir
     write(46,'(a,1p,e13.4)') 'Eeff=',Eeff
     do it=-150,-50
        ee_l=10.d0**(it*0.1d0)
        d=EEdiff(Eeff,ee_l, Te_l,Tz_l)
        dl=EEdiff_ln(log(Eeff),ee_l, Te_l,Tz_l)
        write(46,*) ro_l,Ni_l,Eeff,ee_l, Te_l,Tz_l,d,dl
     end do
  end do

1000 continue
  !==================================================
  !
  ! loop on  density  and temperature, testing RADIOM on all sets of conditions

  LOOPro:   do ir= -8,2	! -4,-4	!
     write(45,*) 
     RO_IN=3*10d0**ir
     ro_l=RO_IN
     Ni_l=ro_l*avogadro/atoMass
     write(*,*)'---- NI=',Ni_l
     ! density dependant , will be passed through MODULE:
     call set_kBr(Natom=Ni_l)

     LOOPte:	   do  it= 0,20	! 0,11	! 4,4	!
        TE_IN=10d0**(it*0.25)
        Te_l=TE_IN
        print*
        !\\\\\\\\\\\\
        !============ Direct EOS, LTE, elec. onlyt (i.e. code uses Te + Ti )
        !////////////
        write(*,*)'- direct EOS , LTE -'
        useLTE=.true.
        print*
        write(*,*)'LTE: ni,te=',Ni_l,Te_l
        Ee_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! reset value, just to check
        RO_in=ro_l
        EE_in=ee_l
        call NLTE_EOS(Natom=Ni_l, Te_in=Te_l, Ee_out=Ee_l, Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pe_out=pe_l, Cv_out=Cv_l)
        LTE_ee_outD=ee_l
        LTE_pe_outD=pe_l
        LTE_cv_outD=cv_l
        LTE_zb_outD=zbar_l
        write(*,*)' -> Ee,Tz=',Ee_l,Tz_l
        write(*,*)' -> Zbar=',zbar_l,' Pe=',Pe_l,' Cv=',Cv_l

        !\\\\\\\\\\\\
        !============ Inverse EOS, LTE
        !////////////
        print*
        write(*,*)'- inverse EOS , LTE-'
        write(*,*)'LTE: ni,Ee=',Ni_l,Ee_l
        te_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! reset value, just to check
        call NLTE_EOS(Natom=Ni_l, EE_IN=Ee_l, TE_out=te_l, Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pe_out=pe_l, Cv_out=Cv_l)
        LTE_te_outI=te_l
        LTE_pe_outI=pe_l
        LTE_cv_outI=cv_l
        LTE_zb_outI=zbar_l
        write(*,*)' -> Te,Tz=',Te_l,Tz_l
        write(*,*)' -> Zbar=',zbar_l,' Pe=',pe_l,' Cv=',Cv_l


        Te_l=TE_IN
        ro_l=RO_IN
        Ni_l=ro_l*avogadro/atoMass
        ! density dependant , will be passed thorugh MODULE:
        call set_kBr(Natom=Ni_l)
        !\\\\\\\\\\\\\
        !============= Direct EOS, LTE  w/  Te=Ti ( this Cv=(Zbar+1)*...
        !/////////////
        print*
        write(*,*)'- direct EOS , LTE -'
        write(*,*)'LTE: ni,(T=te=ti)=',Ni_l,Te_l
        Ee_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! reset value, just to check
        call NLTE_EOS(Natom=Ni_l, Te_in=Te_l, Et_out=Ee_l, Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pt_out=pe_l, Cv_out=Cv_l)
        write(*,*)' -> Et,Tz=',Ee_l,Tz_l
        write(*,*)' -> Zbar=',zbar_l,' Pt=',Pe_l,' Cv=',Cv_l,' TE_in=',TE_IN

        !\\\\\\\\\\\\\
        !============= inverse EOS, LTE  w/  Te=Ti ( this Cv=(Zbar+1)*...
        !/////////////
        print*
        write(*,*)'- inverse EOS , LTE-'
        write(*,*)'LTE: ni,Et=',Ni_l,Ee_l
        te_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! reset value, just to check
        call NLTE_EOS(Natom=Ni_l, Et_IN=Ee_l, TE_out=te_l, Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pt_out=pe_l, Cv_out=Cv_l)
        write(*,*)' -> Te,Tz=',Te_l,Tz_l
        write(*,*)' -> Zbar=',zbar_l,' Pt=',pe_l,' Cv=',Cv_l



        !\\\\\\\\\\\\\
        !============= direct EOS, Non-LTE, elec. only (i.e. code uses Te + Ti )
        !/////////////
166     print*

        ! As the first step, define  Erad, Brad arrays
        ! Where Erad are the group-integrated radiation energy, in arbitrary units,
        ! and Brad are the integrals of the Planckian spectrum, in the same units ,
        ! with the temperature of radiation being equal to the electron temperature,
        ! at the beginning of the time step: even if the EOS is applied not in the
        ! starting time for this timestep.
        ! So that  Brad=Erad in  CTE  conditions. 
        ! It is important not to add some bias because of time or location beeing different.

        
        
        ! note : LTE would be recovered with Erad(:)=Brad(:)

        ro_l=RO_IN
        Ni_l=ro_l*avogadro/atoMass
        ! density dependant , will be passed thorugh MODULE:
        call set_kBr(Natom=Ni_l)
        Te_l=TE_IN
        write(*,*)'- ---------- -'
        print*
        write(*,*)'- direct EOS , nonLTE -'
        write(*,*)'nonLTE: ni,te=',Ni_l,Te_l
        useLTE=.false.
        Ee_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! reset value, just to check

        call NLTE_EOS(Natom=Ni_l,Te_in=Te_l,Ee_out=Ee_l,Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pe_out=pe_l, Cv_out=Cv_l)
        NLTE_ee_outD=ee_l
        NLTE_tz_outD=Tz_l
        NLTE_cv_outD=cv_l
        NLTE_zb_outD=zbar_l
        write(*,*)' -> Ee,Te,Tz=',Ee_l,Te_L,Tz_l	
        write(*,*)' -> Zbar=',zbar_l,' Pe=',pe_l,' Cv=',Cv_l
        !\\\\\\\\\\\\\
        !============= Inverse EOS, Non-LTE , elec. only (i.e. code uses Te + Ti )
        !/////////////
        ro_l=RO_IN
        Ni_l=ro_l*avogadro/atoMass
        ! density dependant , will be passed thorugh MODULE:
        call set_kBr(Natom=Ni_l)
        Te_l=TE_IN
        print*
        write(*,*)'- inverse EOS , nonLTE (not completed) -'
        print*
        write(*,*)'nonLTE: ni,Ee=',Ni_l,Ee_l
        te_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1		! reset value, just to check
        call NLTE_EOS(Natom=Ni_l, Ee_in=Ee_l, Te_out=Te_l, Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pe_out=pe_l, Cv_out=Cv_l)
        write(*,*)'Te_out,Tz_out,zbar=',te_l,tz_l,zbar_l,' TE_in=',TE_IN
        NLTE_te_outI=te_l
        NLTE_tz_outI=tz_l
        NLTE_cv_outI=cv_l
        NLTE_zb_outI=zbar_l
        write(*,*)' -> EE,Te,Tz=',EE_l,Te_l,Tz_l 
        write(*,*)' -> Zbar=',zbar_l,' Pe=',pe_l,' Cv=',Cv_l

	goto 296

 !\\\\\\\\\\\\\
 !============= direct EOS, Non-LTE , el+ion.  (i.e. code uses Te=Ti )
 !/////////////
        ro_l=RO_IN
        Ni_l=ro_l*avogadro/atoMass
        ! density dependant , will be passed thorugh MODULE:
        call set_kBr(Natom=Ni_l)
        Te_l=TE_IN
        write(*,*)'- ---------- -'
        print*
        write(*,*)'- direct EOS , nonLTE -'
        write(*,*)'nonLTE: ni,te=Ti=',Ni_l,Te_l
        useLTE=.false.
        Ee_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! reset value, just to check
        call NLTE_EOS(Natom=Ni_l,Te_in=Te_l,Et_out=Ee_l,Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pt_out=pe_l, Cv_out=Cv_l)
        NLTE_ee_outD=ee_l
        NLTE_tz_outD=Tz_l
        NLTE_cv_outD=cv_l
        NLTE_zb_outD=zbar_l
        eeff=ee_l-1.5d0*kBr_E*(Te_l-Tz_l)*(zbar_l+zion)
        write(*,*)' -> Et,Te,Tz=',ee_l,Te_L,Tz_l,' Tz_new=',Tz_new
        write(*,*)' -> Zbar=',zbar_l,' Pt=',pe_l,' Cv=',Cv_l
        !\\\\\\\\\\\\\
        !============= Inverse EOS, Non-LTE , el+ion.  (i.e. code uses Te=Ti )
        !/////////////
        ro_l=RO_IN
        Ni_l=ro_l*avogadro/atoMass
        ! density dependant , will be passed thorugh MODULE:
        call set_kBr(Natom=Ni_l)
        Te_l=TE_IN
        print*
        write(*,*)'- inverse EOS , nonLTE (not completed) -'
        print*
        write(*,*)'nonLTE: ni,Et=',Ni_l,ee_l
        te_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1		! reset value, just to check
        call NLTE_EOS(Natom=Ni_l, Et_in=ee_l, Te_out=Te_l, Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pt_out=pe_l, Cv_out=Cv_l)
        NLTE_te_outI=te_l
        NLTE_tz_outI=tz_l
        NLTE_cv_outI=cv_l
        NLTE_zb_outI=zbar_l
        write(*,*)' -> Et,Te,Tz=',ee_l,Te_l,Tz_l,' Tz_new=',Tz_new
        write(*,*)' -> Zbar=',zbar_l,' Pt=',pe_l,' Cv=',Cv_l


        !\
        !  write results to check 
        !/
296     write(45,*) ro_in , ee_in, te_in  &
             ,LTE_ee_outD ,LTE_zb_outD ,LTE_pe_outD ,LTE_cv_outD &
             ,LTE_te_outI ,LTE_zb_outI ,LTE_pe_outI ,LTE_cv_outI &
             ,NLTE_ee_outD,NLTE_zb_outD,NLTE_tz_outD,NLTE_cv_outD &
             ,NLTE_te_outI,NLTE_zb_outI,NLTE_tz_outI,NLTE_cv_outI 

     end do LOOPte
  end do LOOPro

  !
  call CON_Stop('- Normal End -')
contains
  subroutine setMethod
    !ONLY used in the test
    use CRASH_M_EOS, ONLY: setOptions
    
    implicit none
    integer ib,ie
    
    write(0,*)'input 2 values : 1 for useZbrent, 1 for useEElog'
    read(*,*) ib,ie
    call setOptions(ib.gt.0,ie.gt.0,.false.)
    
    write(0,*)'calls setZTF(au)'
    
    ! initialize Thomas_Fermi EOS (analytical fit)
    call setZTF('au')
    !===================================================
    !Check the consitency of the input parameters
    call verify()
    return
  end subroutine setMethod
end program testNLTE


subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  write(*,*)StringError
  stop
end subroutine CON_stop
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test
