!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModPowerLawEos
  use ModConst,ONLY: cRadiation
  implicit none
  !This modelu implements the EOS of the following kind:
  !C_V = cESimT4*T^3
  !E/V = cESimT4*T^4/4
  !p   = cESimT4*T^4/12
  !By default the constant is defined as for a pair plasma at T>1 MeV:
  
  real :: cESimT4 = 7.0 * cRadiation 
  !(both electron and positrons have the energey density of 
  !(7/8) of the black radiation
contains
  !=================================
  subroutine read_cesimt4
    use ModReadParam
    call read_var('cESimT4',cESimT4)
  end subroutine read_cesimt4
  !=================================
  subroutine eos_esimt4( TeIn,ETotalIn,pTotalIn,&
       TeOut,ETotalOut,pTotalOut,GammaOut,CVTotalOut)

    real,    optional, intent(in)  :: TeIn          ! temperature SI[K]
    real,    optional, intent(in)  :: ETotalIn      ! internal energy density [J/m^3]
    real,    optional, intent(in)  :: pTotalIn      ! pressure, SI [Pa]

    real,    optional, intent(out) :: TeOut         ! temperature SI[K]
    real,    optional, intent(out) :: pTotalOut     ! pressure, SI [Pa]
    real,    optional, intent(out) :: ETotalOut     ! internal energy density [J/m^3]
    real,    optional, intent(out) :: GammaOut      ! polytropic index
    real,    optional, intent(out) :: CVTotalOut    ! Specific heat per the unit of volume,[J/(K*m^3)]

    real:: T, E, p, CV
    !--------------------

    if(present(TeIn))then

       T = TeIn

    elseif(present(ETotalIn))then

       T = sqrt(sqrt(ETotalIn*4.0/cESimT4))

    elseif(present(pTotalIn))then

       T = sqrt(sqrt(pTotalIn*12.0/cESimT4))
    else

       call CON_stop('In eos for E \sim T^4 none of Te|ETotal|pTotal is among input parameters')

    end if

    if(present(TeOut))      TeOut     = T
    if(present(ETotalOut))  ETotalOut = 0.250*cESimT4*T**4
    if(present(PTotalOut))  pTotalOut = (cESimT4/12.0)*T**4
    if(present(GammaOut))   GammaOut  = 4.0/3
    if(present(CVTotalOut)) CVTotalOut = cESimT4*T**3
  end subroutine eos_esimt4
end module CRASH_ModPowerLawEos
