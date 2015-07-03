!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module Mod_GLOW

  implicit none
  
  integer, parameter  :: NMAJ = 3       !Number of major species
  integer, parameter  :: NEX = 20       !Number of ionized/excited spc
  integer, parameter  :: NW = 20        !Number of airglow emission wavelengths
  integer, parameter  :: NC = 10        !Number of component production terms
  integer, parameter  :: ALTSMAX=400      !for each emission  
  integer, parameter  :: NST = 6        !Number of states produced by photo-
                                          !ionization/dissociation
  integer, parameter  :: NEI = 10       !Number of states produced by e- impact 
  integer, parameter  :: LMAX = 123     !Number of wavelength ints for solar flux
  integer, parameter  :: iO2PP = 1
  integer, parameter  :: iO2DP = 2
  integer, parameter  :: iO4SP = 3
  integer, parameter  :: iN2P  = 4
  integer, parameter  :: iNP   = 5
  integer, parameter  :: iO2P  = 6
  integer, parameter  :: iNOP  = 7
  integer, parameter  :: i_Glow = 15

  integer, allocatable, dimension(:) :: IIMAXX

  real, parameter     :: PI = 3.141592653589793
  integer :: IDate,ipc,ildb,JMAX,NBINS,GFIRST,ETFIRST,EPFIRST
  real    :: UT, GLAT, GLONG,ISCALE,JLOCAL,KCHEM, &
       F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC, SZA, &
       DIP, EFRAC, IERR

  real, allocatable, dimension(:) :: ZZ, ZO, ZN2, ZO2, ZNO,PHITOP
  real, allocatable, dimension(:) :: ZNS, ZND, ZRHO, ZE, ZTN
  real, allocatable, dimension(:) :: ZTI, ZTE, Z, EHEAT
  real, allocatable, dimension(:) :: TEZ, ECALC,XNO,ENER,DEL
  real, allocatable, dimension(:,:) :: ZXDEN,ZMAJ,ZCOL,PESPEC,SESPEC
  real, allocatable, dimension(:,:) :: SION,UFLX,DFLX,ZETA,PIA
  real, allocatable, dimension(:,:) :: SIGS,PE,PIN,PHONO,OUTF
  real, allocatable, dimension(:,:,:) :: SIGA,SEC,SIGEX,SIGIX
  real, allocatable, dimension(:,:,:) :: PHOTOI, PHOTOD,AGLW,ZCETA
 
  real ::  WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX), &
       VCB(NW),waves(LMAX),wavel(LMAX),Interp_Fac(LMAX)

  real :: D(8), T(2), SW(25), &
        OARR(30), TPI(NMAJ)  
          
  data WaveS/                                                            &
       0.50, 1.00, 2.00, 4.00, 8.00, 18.0, 23.0, 32.00,  44.0, 60.0, 70.0, 80.0, 90.0, &
       100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0,    &
       190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0,    &
       280.0, 290.0, 300.0, 310.0, 320.0, 330.0,  340.0, 350.0,    &
       360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0, 430.0, 440.0,    &
       450.0, 460.0, 470.0, 480.0, 490.0, 500.0, 510.0, 520.0, 530.0,    &
       540.0, 550.0, 560.0, 570.0, 580.0,                                &
       590.0, 600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0,    &
       680.0, 690.0, 700.0, 710.0, 720.0, 730.0, 740.0, 750.0,    &
       760.0, 770.0, 780.0, 790.0, 800.0, 810.0, 820.0, 830.0, 840.0,    &
       850.0, 860.0, 870.0, 880.0, 890.0, 900.0, 910.0, 920.0, 930.0,    &
       940.0, 950.0, 960.0, 970.0, 980.0,                                &
       990.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0, 1100.0,    &
       1150.0, 1210.0, 1220.0, 1250.0, 1300.0, 1350.0, 1400.0, 1450.0,   &
       1500.0, 1550.0, 1600.0, 1650.0, 1700.0/

data WaveL/                                                            &
       1.00, 2.00, 4.00, 8.00, 18.0, 23.0, 32.00,                        &
       44.0, 60.0, 70.0, 80.0, 90.0,                                     &
       100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0,    &
       190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0,    &
       280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0,    &
       360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0, 430.0, 440.0,    &
       450.0, 460.0, 470.0, 480.0, 490.0, 500.0, 510.0, 520.0, 530.0,    &
       540.0, 550.0, 560.0, 570.0, 580.0,                                &
       590.0, 600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0,    &
       680.0, 690.0, 700.0, 710.0, 720.0, 730.0, 740.0, 750.0,    &
       760.0, 770.0, 780.0, 790.0, 800.0, 810.0, 820.0, 830.0, 840.0,    &
       850.0, 860.0, 870.0, 880.0, 890.0, 900.0, 910.0, 920.0, 930.0,    &
       940.0, 950.0, 960.0, 970.0, 980.0,                                &
       990.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0, 1100.0,    &
       1150.0, 1210.0, 1220.0, 1250.0, 1300.0, 1350.0, 1400.0, 1450.0,   &
       1500.0, 1550.0, 1600.0, 1650.0, 1700.0, 1750.0/

 DATA SW/25*1./

contains
  !===========================================================
  subroutine init_mod_glow

   
    if(allocated(PHITOP)) deallocate(PHITOP)
    if(allocated(IIMAXX)) deallocate(IIMAXX)
    if(allocated(ZMAJ))   deallocate(ZMAJ)
    if(allocated(ZCOL))   deallocate(ZCOL)
    if(allocated(PESPEC)) deallocate(PESPEC)
    if(allocated(SESPEC)) deallocate(SESPEC)
    if(allocated(SION))   deallocate(SION)
    if(allocated(PIA))    deallocate(PIA)
    if(allocated(UFLX))   deallocate(UFLX)
    if(allocated(DFLX))   deallocate(DFLX)
    if(allocated(AGLW))   deallocate(AGLW)
    if(allocated(ZETA))   deallocate(ZETA)
    if(allocated(EHEAT))  deallocate(EHEAT)
    if(allocated(TEZ))    deallocate(TEZ)
    if(allocated(ECALC))  deallocate(ECALC)
    if(allocated(ZCETA))  deallocate(ZCETA)
    if(allocated(SIGS))   deallocate(SIGS)
    if(allocated(PE))     deallocate(PE)
    if(allocated(PIN))    deallocate(PIN)
    if(allocated(SIGA))   deallocate(SIGA)
    if(allocated(SEC))    deallocate(SEC)
    if(allocated(SIGEX))  deallocate(SIGEX)
    if(allocated(SIGIX))  deallocate(SIGIX)
    if(allocated(PHOTOI)) deallocate(PHOTOI)
    if(allocated(PHOTOD)) deallocate(PHOTOD)
    if(allocated(XNO))    deallocate(XNO)
    if(allocated(OUTF))   deallocate(OUTF)
    if(allocated(ENER))   deallocate(ENER)
    if(allocated(DEL))    deallocate(DEL)
    if(allocated(PHONO))  deallocate(PHONO)
    if(allocated(ZTN))    deallocate(ZTN)
    if(allocated(ZTE))    deallocate(ZTE)
    if(allocated(ZTI))    deallocate(ZTI)
    if(allocated(ZRHO))   deallocate(ZRHO)
    if(allocated(ZO))     deallocate(ZO)
    if(allocated(ZO2))    deallocate(ZO2)
    if(allocated(ZN2))    deallocate(ZN2)
    if(allocated(ZNO))    deallocate(ZNO)
    if(allocated(ZNS))    deallocate(ZNS)
    if(allocated(ZND))    deallocate(ZND)
    if(allocated(ZXDEN))  deallocate(ZXDEN)
    if(allocated(ZE))     deallocate(ZE)

 
    allocate(PHITOP(NBINS))
    allocate(IIMAXX(NBINS))
    allocate(ZMAJ(NMAJ,JMAX))
    allocate(ZCOL(NMAJ,JMAX))
    allocate(PESPEC(NBINS,JMAX))
    allocate(SESPEC(NBINS,JMAX))
    allocate(SION(NMAJ,JMAX))
    allocate(PIA(NMAJ,JMAX))
    allocate(UFLX(NBINS,JMAX))
    allocate(DFLX(NBINS,JMAX))
    allocate(AGLW(NEI,NMAJ,JMAX))
    allocate(ZETA(NW,JMAX))
    allocate(EHEAT(JMAX))
    allocate(TEZ(JMAX))
    allocate(ECALC(JMAX))
    allocate(ZCETA(NC,NW,JMAX))
    allocate(SIGS(NMAJ,NBINS))
    allocate(PE(NMAJ,NBINS))
    allocate(PIN(NMAJ,NBINS))
    allocate(SIGA(NMAJ,NBINS,NBINS))
    allocate(SEC(NMAJ,NBINS,NBINS))
    allocate(SIGEX(NEI,NMAJ,NBINS))
    allocate(SIGIX(NEI,NMAJ,NBINS))
    allocate(PHOTOI(NST,NMAJ,JMAX))
    allocate(PHOTOD(NST,NMAJ,JMAX))
    allocate(XNO(JMAX))
    allocate(OUTF(11,JMAX))
    allocate(ENER(NBINS))
    allocate(DEL(NBINS))
    allocate(PHONO(NST,JMAX))
    allocate(ZTN(JMAX))
    allocate(ZTE(JMAX))
    allocate(ZTI(JMAX))
    allocate(ZRHO(JMAX))
    allocate(ZO(JMAX))
    allocate(ZO2(JMAX))
    allocate(ZN2(JMAX))
    allocate(ZNO(JMAX))
    allocate(ZNS(JMAX))
    allocate(ZND(JMAX))
    allocate(ZXDEN(NEX,JMAX))
    allocate(ZE(JMAX))
    
   
end subroutine init_mod_glow
 
end module Mod_GLOW


   
