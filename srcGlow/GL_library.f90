!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------
subroutine GL_settime(itime,utime)
 
  use Mod_Glow
  implicit none
  integer, intent(in) :: itime(7)
  real, intent(in)    :: utime
  integer             :: year,month,day,yyyddd
  year = itime(1)
  month = itime(2)
  day = itime(3)
  ut = utime
  
  call ymd_to_ydoy(year,month,day,IDate)
  
end subroutine GL_settime

!------------------------------------------

subroutine GL_setloc(alt,n,lat,lon,npb,ierror)

  use Mod_Glow
  implicit none
  integer, intent(in)  :: n,npb
  real, intent(in)     :: alt(n),lat,lon
  integer,intent(out)  :: iError
  
  if (allocated(ZZ)) deallocate(ZZ)
  allocate(ZZ(n),stat=iError)
  ZZ = alt * 100.
  JMAX = n
  NBINS = npb
  GLAT = lat*180/PI
  GLONG = lon*180/PI

  if (JMAX.gt.ALTSMAX) then
     iError=1 
  else 
     iError = 0
  endif
end subroutine GL_setloc
  
!-----------------------------------------------------

subroutine GL_init
  use Mod_Glow
  implicit none
  
  call init_mod_glow
  PIA(:,:) = 0
  
end subroutine GL_init

!------------------------------------------------------

subroutine GL_interp_flux(GITM_wave_l,GITM_wave_s,GITM_flux,nGitmWaves)
  use Mod_Glow
  implicit none
  
  integer, intent(IN) :: nGitmWaves
  real, intent(IN)    :: GITM_wave_l(nGitmWaves), GITM_wave_s(nGitmWaves)
  real, intent(IN)    :: GITM_flux(nGitmWaves)

  integer          :: i,j, istart, ifinish, incL, incS, countL,countS
  integer          :: countdiff,icount,incU,incD
  real             :: S_Factor, L_Factor
  real, dimension(nGitmWaves)  :: dwaveL,dwaveS
  real :: temp1, temp2, temp3
  integer :: l

  !Interpolate GITM flux to solar flux grid used by GLOW
  ! This is complicated because GITM uses individual lines, ie,
  ! wavelow = wavehigh for some bins.  GLOW does not like this
  ! so it is necassary to not use those particular bins in the 
  ! interpolation.  Their energy is added to the bin after the 
  ! interp coefficients are determined.

  !Note: Bins are arranged from high wavelength to short wavelength



  do l = 1,lmax
     wave1(l) = waveS(lmax-(l-1))
     wave2(l) = waveL(lmax-(l-1))
  enddo
  
  SFLUX = 0.0

  do i = 1, LMAX

     if (wave2(i) .ge. GITM_wave_l(1)) then 
        istart = i
     endif

     if (WAVE1(i) .ge. GITM_wave_s(nGitmWaves)) then 
        ifinish = i-1
     endif

  enddo

  do i = istart, ifinish
     
     dwaveL = GITM_wave_l - wave2(i)
     dwaveS = GITM_wave_s - wave1(i)

     where (dwaveL .lt. 0) dwaveL = 1000
     where (dwaveS .gt. 0) dwaveS = -1000
     countL = 1
     countS = 1
     do j = 1, nGitmWaves  
        if (dwaveL(j) .lt. dwaveL(countL)) countL = j
     enddo
     do j = 1, nGitmWaves  
        if (dwaveS(j) .gt. -900) then
           countS = j
           goto 300
        endif
     enddo

300  incL = 0
     incS = 0
     incU = 0
     incD = 0

100  if (GITM_wave_l(countL+incL) .eq. GITM_wave_s(countL+incL)) then
        if (dwaveL(countL+incU) .ge. dwaveL(countL-incD) .or. &
           GITM_wave_l(countL+incU+1) .lt. wave2(i)) then
           incD = incD - 1
           incL = incD
        else 
           incU = incU + 1
           incL = incU
        endif
        goto 100
     endif
    
     incU = 0
     incD = 0

200  if (GITM_wave_s(countS+incS) .eq. GITM_wave_l(countS+incS)) then
        if (dwaveS(countS+incU) .le. dwaveS(countS-incD).or. &
           GITM_wave_s(countS+incU+1) .gt. WAVE1(i)) then
           incD = incD - 1
           incS = incD
        else 
           incU = incU + 1
           incS = incU
        endif
        goto 200
     endif
     
     countS = countS+incS
     countL = countL+incL
     countdiff = countS - countL
     if(countdiff .eq. 0) then
        
        L_factor = (GITM_wave_l(countL) - wave2(i))/ &
             (GITM_wave_l(countL) - GITM_wave_s(countL))
        S_factor = (GITM_wave_l(countL) - WAVE1(i))/ &
             (GITM_wave_l(countL) - GITM_wave_s(countL))
        SFLUX(i) = GITM_flux(countL)*(S_Factor-L_factor) 
        
     else if (countdiff .ge. 1) then
        L_factor = (GITM_wave_l(countL) - wave2(i))/ &
             (GITM_wave_l(countL) - GITM_wave_s(countL))
        S_factor = (GITM_wave_l(countS) - WAVE1(i))/ &
             (GITM_wave_l(countS) - GITM_wave_s(countS))
        SFLUX(i) = GITM_flux(countL)*(1-L_factor) + &
             GITM_flux(countS) *S_Factor         
        
        do icount = 1, countdiff-1
           if (GITM_wave_l(countL+icount) .ne. &
                GITM_wave_s(countL+icount)) then
              SFLUX(i) = SFLUX(i) + GITM_flux(countL+icount)
           endif
        enddo
     endif
  enddo
  
  !Deal with beginning and end of GLOW spectrum if not already
  if (istart .ne. 1) then
     L_factor = (GITM_wave_l(1) - wave2(1))/ &
          (GITM_wave_l(1) - GITM_wave_s(1))
     S_factor = (GITM_wave_s(1) - WAVE1(1))/ &
          (GITM_wave_s(1) - 0)
     SFLUX(i) = GITM_flux(1)*(S_Factor+L_factor)
  endif

  if (ifinish .ne. LMAX) then
     L_factor = (wave2(LMAX)-GITM_wave_l(nGitmWaves))/ &
          (GITM_wave_l(nGitmWaves) - GITM_wave_s(nGitmWaves))
     S_factor = (GITM_wave_l(nGitmWaves) - WAVE1(LMAX))/ &
          (GITM_wave_l(nGitmWaves) - GITM_wave_s(nGitmWaves))
     SFLUX(i) = GITM_flux(nGitmWaves)*(S_Factor+L_factor)
  endif
  !Now take care of individual lines in GITM
  do j = 1, nGitmWaves 
     if (GITM_wave_l(j) .eq. GITM_wave_s(j)) then
        do i = 1, LMAX 
           if (WAVE1(i) .le. GITM_wave_l(j) .and. &
                wave2(i) .ge. GITM_wave_l(j)) &
                SFLUX(i) = SFLUX(i) + GITM_flux(j)
          
        enddo
     endif
  enddo


  !Convert to cm
  SFLUX = SFLUX/(10000.0)
  

end subroutine GL_interp_flux
!------------------------------------------------------  

subroutine GL_setND(rho,o,n2,o2,no,ns,nd)

  use Mod_Glow
  implicit none
  real, intent(in), dimension(JMAX) :: rho,o,n2,o2,no,ns,nd

  ZRHO = rho
  ZO  = o
  ZO2 = o2
  ZN2 = n2
  ZNO = no
  ZNS = ns
  ZND = nd

end subroutine GL_setND

!-------------------------------------------------------

subroutine GL_settemp(tn,ti,te)

  use Mod_GLOW
  implicit none
  real, intent(in) :: tn(JMAX),ti(JMAX),te(JMAX)
  
  ZTN = tn
  ZTE = te
  ZTI = ti

end subroutine GL_settemp

!-------------------------------------------------------

subroutine GL_setID(op, o2p, nop, n2p, np, e)
  
  use Mod_GLOW
  implicit none
  real, intent(in), dimension(JMAX) :: op,o2p,nop,e,n2p,np 
  integer :: j
  
  ZXDEN(:,:) = 0
  ZXDEN(iO4SP,:) = op
  ZXDEN(iO2P,:)  = o2p
  ZXDEN(iNOP,:)  = nop
  ZXDEN(iN2P,:)  = n2p
  ZXDEN(iNP,:)   = np

  ZE = e

end subroutine GL_setID

!-------------------------------------------------------

subroutine GL_getvals(ipin,emisrate,upeflx,downeflx,photoErate,ebins)
  
  use Mod_GLOW
  implicit none
  integer, intent(in) :: ipin
  real, intent(out), dimension(NW,JMAX   ) :: emisrate
  real, intent(out), dimension(NBINS,JMAX) :: upeflx,downeflx,photoErate
  real, intent(out), dimension(NBINS)      :: ebins

  real :: STL,XMLAT,SZAC,DIPC,GLATC,GLONGC,SZAD,DIPD,SZACD,DIPCD,AD
  real :: GLATS, GLONGS, AP, TOTPI,TOTSI, XMLONG
  integer :: ICONJ,i, j, NS, iw,N

  ipc = ipin
  STL = (UT/240.+GLONG) / 15.
  IF (STL .LT. 0.) STL = STL + 24.
  IF (STL .GT. 24.) STL = STL - 24.


  !!!!!KCHEM should be set to 1 since this means that GLOW doesn't do any unnessacary 
  ! chemistry calculations.  However, the O2P density from GITM is funny just below 
  ! F region alts, so we use KCHEM = 3 to make glow values better.

  ISCALE = 2
  XUVFAC = 3.
  JLOCAL = 0
  KCHEM = 1
  ICONJ = 0
  HLYBR = 0.
  FEXVIR = 0.
  HLYA = 0.
  PHITOP(:) = 0.0
  XMLONG = 0.0
  XMLAT = 0.0

  GFIRST = 1
  ETFIRST = 1
  EPFIRST = 1
  ZETA(:,:) = 0
 
  call GLOW
   

  GLATS = GLAT
  GLONGS = GLONG

  CALL GEOMAG_glow(0,GLONG,GLAT,XMLONG,XMLAT)

  CALL GEOMAG_glow(1,GLONG,GLAT,XMLONG,-XMLAT)


!Downward flux at conjugate point = upward flux at location
  do N = 1, NBINS
     PHITOP(N) = UFLX(N,JMAX)
  enddo

!Call Glow for conjugate point (same neutral atmosphere)
  
  If (ICONJ .eq. 1) call GLOW

!Upward flux at conjugate point = downward flux at location
  do N = 1, NBINS
     PHITOP(N) = UFLX(N,JMAX)
  enddo

!Call Glow again at location
  SZAC = SZA
  DIPC = DIP
  GLATC = GLAT
  GLONGC = GLONG
  GLAT = GLATS
  GLONG = GLONGS

  IF (ICONJ .eq. 1) call GLOW

 SZAD = SZA * 180. / PI
      DIPD = DIP * 180. / PI
      SZACD = SZAC * 180. / PI
      DIPCD = DIPC * 180. / PI
!      write (6,444) IDATE, UT, GLAT, GLONG, F107, F107A, AP
!  444 FORMAT (' Date=',i5,' UT=',f6.0,' Lat=',f5.1,' Lon=',f6.1, &
!             ' F107=',f4.0,' F107A=',f4.0,' Ap=',f4.0)
!      WRITE (6,445) SZAD, STL, DIPD, GLATC, GLONGC, SZACD, EFRAC, IERR
!  445 FORMAT (' SZA=',F5.1,' LST=',F5.2,' Dip=',F5.1,' Clat=',F5.1, &
!        ' Clon=',F6.1,' CSZA=',F5.1,' Ec=',F6.3,' Ie=',I1)
!
! write (6,690)
!  690 format ('  Z     Photoion   EIion    Ecalc     O+(2P)    ',&
!           'O+(2D)    O+(4S)     N+         N2+       O2+       NO+')
      do 750 j=1,JMAX
        do 700 i=1,nmaj
          tpi(i) = 0.
          do 700 ns=1,nst
            tpi(i) = tpi(i) + photoi(ns,i,j)
  700     continue
        totpi = tpi(1) + tpi(2) + tpi(3) + phono(1,j)
        totsi = sion(1,j) + sion(2,j) + sion(3,j)
!        write (6,730) z(j),totpi,totsi,ecalc(j),(zxden(i,j),i=1,7)
!  730   format (1x, 0p, f5.1, 1p, 10e10.2)
  750 continue

 !write (6,780)
 ! 780 format ('   z     3371   4278   5200   5577   6300',&
 !          '   7320  10400   3466   7774   8446')
 !     write (6,790) (z(j), (zeta(iw,j),iw=1,10), j=1,jmax)
 ! 790 format (1x, f5.1, 10f7.1)
 !     write (6,795)  (vcb(iw),iw=1,10)
 ! 795 format (' VCB:',11f7.0)

emisrate = ZETA

upeflx = UFLX
downeflx = DFLX
ebins = ENER

photoErate = PESPEC

!do j=1,JMAX
!   write(*,*) ZXDEN(iO4SP,J), ZZ(J)   
!   write(*,*) ZETA(5,J), ZZ(j),j
!enddo
! do j=1,LMAX
!    write(*,*) SFLUX(j)
! end do




end subroutine GL_getvals

!------------------------------------------------------

subroutine GL_clean

  use Mod_Glow
  implicit none
  

  if (allocated(ZZ)) deallocate(ZZ)
  if (allocated(ZRHO)) deallocate(ZRHO)
  if (allocated(ZO)) deallocate(ZO)
  if (allocated(ZO2)) deallocate(ZO2)
  if (allocated(ZN2)) deallocate(ZN2)
  if (allocated(ZNO)) deallocate(ZNO)
  if (allocated(ZNS)) deallocate(ZNS)
  if (allocated(ZND)) deallocate(ZND)
  if (allocated(ZTN)) deallocate(ZTN)
  if (allocated(ZTE)) deallocate(ZTE)
  if (allocated(ZTI)) deallocate(ZTI)
  if (allocated(ZXDEN)) deallocate(ZXDEN)
  if (allocated(ZE)) deallocate(ZE)
  if (allocated(Z)) deallocate(Z)
  if(allocated(IIMAXX)) deallocate(IIMAXX)
  if(allocated(ZMAJ))   deallocate(ZMAJ)
  if(allocated(ZCOL))   deallocate(ZCOL)
  if(allocated(PESPEC)) deallocate(PESPEC)
  if(allocated(SESPEC)) deallocate(SESPEC)
  if(allocated(SION))   deallocate(SION)
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
  if(allocated(PHITOP)) deallocate(PHITOP)
end subroutine GL_clean


!-------------------------------------------------------


subroutine ymd_to_ydoy(year,month,day,yyyyddd)

  implicit none
  integer, intent(in) :: year,month,day
  integer,intent(out) :: yyyyddd


  integer, dimension(1:12) :: dayofmon
  integer :: doy, imonth

  dayofmon(1) = 31
  dayofmon(2) = 28
  dayofmon(3) = 31
  dayofmon(4) = 30
  dayofmon(5) = 31
  dayofmon(6) = 30
  dayofmon(7) = 31
  dayofmon(8) = 31
  dayofmon(9) = 30
  dayofmon(10) = 31
  dayofmon(11) = 30
  dayofmon(12) = 31
  
  if (mod(year,4).eq.0) dayofmon(2) = dayofmon(2) + 1

  doy = 0
  do imonth = 1, month-1
     doy = doy + dayofmon(imonth)
  enddo
  
  doy = doy + day
 
  yyyyddd = (year * 1000) + doy

end subroutine ymd_to_ydoy
