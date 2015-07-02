!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=====================Above is for ZTF
!HERE CAL means CACculated and tabulated EOS in contrast with
!ZTF which is the analytical fit for the Thomas-Fermi model
module CRASH_M_EOS
  !The CRASH meterial number
  SAVE
  integer:: iMaterial = -1
  logical:: UseCrashEos = .false.
contains
  !------
  subroutine set_kbR(ro,Natom)
    use CRASH_M_localProperties,only : kBr_E,kBr_P ,ERGperEV,DYNEperEV &
         ,avogadro,EVperK,Atomass,kB_ztf_E,kB_ztf_P	! ,ro,NI
    use ModConst,ONLY: cBoltzmann,cEV,cKToEv,cEVToK, cAtomicMass
    use CRASH_ModEos,ONLY:cAtomicMassCRASH_I
    implicit none
    real,optional,intent(IN) :: ro,Natom
    if(useCrashEos)then
       if(present(ro) ) then
          call CON_stop('set_kbr for CRASH requires nAtomic as as input')
          kBr_E= ro/(cAtomicMass * cAtomicMassCRASH_I(iMaterial))& !Natomic
               * cBoltzmann * cEVToK 
          kBr_P= kBr_E
       else
          kBr_E= Natom * (1.0e7*cBoltzmann) * cEVToK  !erg/(cm3*eV) 
          kBr_P= Natom * (1.0e7*cBoltzmann) * cEVToK 
       end if
    else
       if(present(ro) ) then
          kBr_E= ro * kB_ztf_E
          kBr_P= ro * kB_ztf_P
       else
          kBr_E= Atomass * (Natom/avogadro) * kB_ztf_E
          kBr_P= Atomass * (Natom/avogadro) * kB_ztf_P
       end if
    end if
    return
  end subroutine set_kBr
  !=======================
  subroutine setOptions(brent,EElog,caleos)
    use CRASH_M_NLTE,only : useZbrent,useEElog
    implicit none
    logical,optional,intent(IN) :: brent,EElog,caleos
    
    if(present(brent)) useZbrent=brent
    if(present(EElog)) useEElog=EElog
    ! write(0,*) '- - useEElog,useZbrent=',useEElog,useZbrent
    return
  end subroutine setOptions
end module CRASH_M_EOS
!===============
!  to develop and test "radiom.f90"  we use this simple EOS based on
!  Thomas-Fermi  fit  by  F.Perrot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! \
!!   PROBABLY IT SHOULD BE BETTER TO PASS "ro" through argument all along the hierachy of calls
!! /
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!===below are the interface routines 
!------
subroutine verify()	! check transmission of Z,A,...
  use CRASH_M_localProperties
  implicit none
  write(*,*)'verify: Z,A=',atonum,atomass,' roS=',roSolid
  if(atonum.gt.0  .and. atomass.gt.0) return
  stop '-P-  atoNum or atoMass not set'
end subroutine verify
!==================================

subroutine LTE_EOS_dir(te,Etot,Ptot,Zbar,Cv)	! ro : in module M_localProperties
  use CRASH_M_ZTF,only : ZTF_EOS_dir,setCALEOS  ! 
  use CRASH_M_localProperties,only : atoNum,ro,zion !ModLocalProperties.f90
  use CRASH_ModEos, ONLY: eos   
  use CRASH_M_EOS,ONLY:iMaterial,useCALEOS=>UseCrashEos
  use ModConst
  implicit none
  real,intent(IN) :: te   !in eV
  real,intent(OUT) :: Etot,Ptot,Zbar  !erg/cm3, dyne/cm2, 1
  real, optional, intent(OUT) :: Cv ! erg/(cm3*eV)
  real,external :: zbrent_EE_inv

  if(useCALEOS) then
      if(zion==0) then !The electron contribution only
        call eos(iMaterial=iMaterial,&
                 Rho=ro*1.0e3,  &   !ro is in g/cm3
                 TeIn=Te*ceVToK,&   !Te is in eV
                 eElectronOut=ETot,&
                 pElectronOut=PTot,&
                 zAverageOut=zBar,&
                 CvElectronOut=Cv)
     else
        call eos(iMaterial=iMaterial,&
                 Rho=ro*1.0e3,  &   !ro is in g/cm3
                 TeIn=Te*ceVToK,&   !Te is in eV
                 eTotalOut=ETot,&
                 pTotalOut=PTot,&
                 zAverageOut=zBar,&
                 CvTotalOut=Cv)
     end if 
     !\
     !CONVERT
     !/
     ETot = ETot * 10.0   !J/m3 = 10^7 erg/10^6cm3 = 10 erg/cm3
     PTot = PTot * 10.0   !J/m3 = 10^7 erg/10^6cm3 = 10 erg/cm3
     Cv   = Cv   * 10.0   !J/m3 = 10^7 erg/10^6cm3 = 10 erg/cm3
     Cv = Cv * cEvToK     !erg/cm3K = erg/cm^3eV *(eV/K)
  else
     call ZTF_EOS_dir(te,Etot,Ptot,Zbar,Cv)
  end if
  !Check positivity
 
  if(pTot<0.0)then
     write(*,*)'Negative pressure in direct LTE_EOS =',ptot,' dyne/cm2'
     write(*,*)'TeIn=',Te,' eV'
     
     write(*,*)'eTot=',eTot,' erg/cm3'
     
  
     write(*,*)'zBar=',zBar
  
     
     call CON_stop('Unphysical output parameter(s) from dir LTE EOS')
  end if
end subroutine LTE_EOS_dir	! ro : in module M_localProperties

!------
subroutine LTE_EOS_inv(te,Etot,Ptot,Zbar,Cv)	! ro : in module M_localProperties
  use CRASH_M_localProperties,only : ro, zion
  use CRASH_M_ZTF,only : ZTF_EOS_inv,setCALEOS
  use CRASH_ModEos, ONLY: eos
  use CRASH_M_EOS,ONLY:iMaterial,useCALEOS=>UseCrashEos
  use ModConst, ONLY:cKToEv,cEVToK
  implicit none
  real,intent(IN) :: Etot
  real,intent(OUT) :: te,Ptot,Zbar
  real,optional,intent(OUT) :: Cv

  if(useCALEOS) then
     if(zion==0) then !The electron contribution only
        call eos(iMaterial=iMaterial,&
                 Rho=ro*1.0e3,         &  !ro is in g/cm3
                 eElectronIn=ETot*0.10,&  !eTot is in erg/cm3
                 TeOut=Te,&
                 pElectronOut=PTot,&
                 zAverageOut=zBar,&
                 CvElectronOut=Cv)
     else
        call eos(iMaterial=iMaterial,&
                 Rho=ro*1.0e3,         &  !ro is in g/cm3
                 eTotalIn=ETot*0.10,   &  !eTot is in erg/cm3
                 TeOut=Te,&
                 pTotalOut=PTot,&
                 zAverageOut=zBar,&
                 CvTotalOut=Cv)
     end if
     !\
     !CONVERT
     !/
     Te   = Te * cKToeV   
     PTot = PTot * 10.0   !J/m3 = 10^7 erg/10^6cm3 = 10 erg/cm3
     Cv   = Cv   * 10.0   !J/m3 = 10^7 erg/10^6cm3 = 10 erg/cm3
     Cv   = Cv * cEvToK   !erg/cm3K = erg/cm^3eV *(eV/K)
  else
     call ZTF_EOS_inv(te,Etot,Ptot,Zbar,Cv)
  end if
  !Check positivity
 
  if(pTot<0.0)then
     write(*,*)'Negative pressure in inverse LTE_EOS =',ptot,' dyne/cm2'
     write(*,*)'TeIn=',Te,' eV'
     
     write(*,*)'eTot=',eTot,' erg/cm3'
     
  
     write(*,*)'zBar=',zBar
  
     
     call CON_stop('Unphysical output parameter(s) from inv LTE EOS')
  end if

end subroutine LTE_EOS_inv	! ro : in module M_localProperties

!------
subroutine getEcold(ro,EEcold,TEcold)
  use CRASH_M_ZTF,only : getEcoldZTF,setCALEOS
  use CRASH_M_localProperties,only : atoNum
  use CRASH_M_Eos,            ONLY : useCALEOS=>UseCrashEos

  implicit none
  real,intent(IN) :: ro
  real,intent(OUT) :: EEcold,TEcold
  real :: Zcold,CVcold,Pcold

  if(useCALEOS) then
     !At the time being, Ecold=0 and for
     !any positive energy density the LTE regime
     !is not applied. In the future, if we know
     !some good estimate, for Ecold(rho) below which
     !the LTE approximation can be employed, we can save
     !the computation time at low energy densities (or
     !temperatures)
     EEcold=0
     TEcold=0
  else
     call getEcoldZTF(atoNum,ro,EEcold,TEcold)
  end if
end subroutine getEcold
!============================================

