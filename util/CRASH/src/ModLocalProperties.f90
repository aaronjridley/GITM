!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
 ! *** mod_local.f90 ***
 !------
module CRASH_M_localProperties
  implicit none

  real,save :: zion = 1 ! 0 if code use 2 T (Te + Ti)
  real :: atoNum  !Atomic number
  real :: Atomass !Atomic mass
  real :: roSolid !Parameter of the EOS function

  !  Eflor, TEfloor may vary with the EOS in use. 
  !  Are to be set by the calling program
  real :: Efloor=1d-5   !Minimal energy density
  real :: TEfloor=1./40.!Minimal electron temperature

  real :: Pcold=0	!Parameter of the EOS function

  !\
  ! density dependant , can be passed through MODULE:
  !/
  real :: kBr_E,kBr_P  ! conversion factor T -> E [erg/cm3], or P [dyne], includes "ro"
  !  kBr_E, kbR_P  are dependant on the units used.  kinetic contribution to 
  !  energy and pressure,
  !  they are generally proportional to the actual bulk density  ro
  !  are   3/2*kBr_E * Zp * T  and kBr_P * Zp * T  ,  
  !  zp =Zbar if only elec. contrib,  zp=Zbar+1   with  elec+ion contib. (w/ Te=Ti)


  real :: ro 	       ! bulk density [g/cm3] shared by all routines, and unchanged
  real :: Ni           ! ionic density [cm-3]

  real,parameter :: kB_ztf_E=1.8591817e-8
  real,parameter :: kB_ztf_P=kb_ztf_E
  real,parameter :: ErgPEReV=1.60219e-12
  real,parameter :: DYNEperEV=ErgPEReV 
  real,parameter :: avogadro=6.02e23
  real,parameter :: EVperK = 11604.1
	
  ! flags set by CALL setRoNi
  ! logical,save :: roGiven=.false.,niGiven=.false.

  real,parameter :: Zsmall=0.050

 !  Zsmall is a lower bound to be used in the  Eef eq. (see correctEOS and EEdiff)

end module CRASH_M_localProperties
!======================
subroutine set_ZA(Z,A)
  use CRASH_M_localProperties,only : atoNum,atoMass
  implicit none
  real,intent(in) :: Z,A
  atoNum=Z
  atoMass=A
end subroutine set_ZA
