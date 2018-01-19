subroutine calc_electron_temperature(iBlock,eHeatingp,iHeatingp,eHeatingm,iHeatingm,iHeating, lame, lami)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModConstants
  use ModTime
  use ModInputs
  use ModUserGITM
  
  implicit none

  integer, intent(in) :: iBlock
  real(kind=8), intent(in), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: eHeatingp, iHeatingp, eHeatingm, &
       iHeatingm, iHeating, lame, lami
 
  real, dimension(nLons,nLats,0:nAlts+1) :: tn, te, ti, etemp, itemp, nn, ni, ne
  real, dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: lam_e, lam_i, sinI2, sinI
  real(kind=8), dimension(nLons,nLats,nAlts) :: eConduction, iConduction

  real(kind=8), dimension(nLons,nLats,nAlts) :: eThermo=0.0, eHeatadv=0.0, eAdiab=0.0

  real, dimension(nLons,nLats,0:nAlts+1) :: alts 
  real :: tipct = 1.1

! calc thermoelectric currents                                                                                                              
  real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2) :: PartialJPara
  real, dimension(-1:nLons+2,-1:nLats+2) :: JParaAlt

!!! For Electron Heat Flux
  real, dimension(nLons,nLats) :: eflux, mlats                    ! W/m2 

!!!! Electron heat flux from Yue

  real :: Dst,Kp,bb1,kk1,P1,P2,bb2,kk2
  real :: x1,y1,z1,aa,bb,cc,temp,mlon,mlat
  real :: geo_lat,geo_lon
  integer, dimension(7) :: iTime
  integer, external :: jday
  integer :: DoY

! for Tridiagnal solver
  real(kind=8), dimension(nAlts) :: a, b, c, d, u, testm     !
  real(kind=8) :: zu, zl, lamu, laml, nne, tte, nni, tti, neu, nel, uiu, uil
  real(kind=8) :: xcoef, hcoef, fcoef, ilam, m
  
  real, dimension(nLons,nLats,0:nAlts+1) :: uiup

  integer :: iLon, iLat, iAlt

  tn = Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock)*TempUnit(1:nLons,1:nLats,0:nAlts+1)
  nn = NDensity(1:nLons,1:nLats,0:nAlts+1,iBlock)
  ne = IDensityS(1:nLons,1:nLats,0:nAlts+1,ie_,iBlock)
  ni = IDensityS(1:nLons,1:nLats,0:nAlts+1,ie_,iBlock)
  uiup = IVelocity(1:nLons,1:nLats,0:nAlts+1,iUp_,iBlock)  
  alts = Altitude_GB(1:nLons,1:nLats,0:nAlts+1,iBlock)  

  
  do iLon = -1, nLons+2
    do iLat = -1, nLats+2
      do iAlt = -1, nAlts+2
        eTemperature(iLon,iLat,iAlt,iBlock) = &
            0.5*(1000.0 - &
                 Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt))*&
             (1.0 + tanh( (Altitude_GB(iLon,iLat,iAlt,iBlock) - 1050.0e+03)/100.0e+03)) + &
                 Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt)

!         iTemperature(iLon,iLat,iAlt,iBlock) = &
!            0.5*(eTemperature(iLon,iLat,iAlt,iBlock) + &
!                  Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt))

        iTemperature(iLon,iLat,iAlt,iBlock) = &
            0.5*(eTemperature(iLon,iLat,iAlt,iBlock) - &
                 Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt))*&
             (1.0 + tanh( (Altitude_GB(iLon,iLat,iAlt,iBlock) - 1900.0e+03)/200.0e+03)) + &
                 Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt)
      enddo !iAlt = -1, nAlts+2
    enddo !iLat = -1, nLats+2
  enddo !iLon = -1, nLons+2

!!!! Upper Boundary Conditions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!![Reference: 1. J.T. Hastings and R.G. Roble, Planet. Space Sci., Vol.25, pp.209, 1977.  Equation(31)
!!!!!!!           2. M.W. Liemohn et al, J. Geohys. Res., Vol.105, pp.27,767, 2000   Plate 1.            ]




!  Kp=2.0
!  Dst=-Kp*20.0
!  
!  call time_real_to_int(CurrentTime, iTime)
!  DoY = jday(iTime(1), &
!       iTime(2), &
!       iTime(3)) 
!  
!  bb1=2365-508*exp(-((DoY-jday(1998,6,12))/180.)**2)
!  kk1=7.8+6.2*sin(DoY*pi/365.)
!  P1=1200.+ 500.*(1+cos(DoY*4*pi/365.))
!  P2=4200.+ 650.*(1+cos(DoY*4*pi/365.))
!  bb2=64.+3.*sin(DoY*3*pi/365.)
!  kk2=1.6+0.5*sin(DoY*3*pi/365.)
! 
!  x1=16.0
!  y1=50
!  z1=7.0
!  aa=10.0
!  bb=20.0
!  cc=2.5
!
!  do iLat=1,nLats
!     do iLon=1,nLons
!        
!        !if (MLT(iLon,iLat,nAlts) < 0.0) then
!        !   MLT(iLon,iLat,nAlts) = MLT(iLon,iLat,nAlts) + 24.0
!        !endif
!        
!        !if (abs(MLT(iLon,iLat,nAlts)-x1)<12) then
!        !   temp=1-(MLT(iLon,iLat,nAlts)-x1)**2/aa**2- &
!        !        (abs(MLatitude(iLon,iLat,nAlts,iBlock))-y1)**2/bb**2
!        !else
!        !   temp=1-(24-abs(MLT(iLon,iLat,nAlts)-x1))**2/aa**2- &
!        !        (abs(MLatitude(iLon,iLat,nAlts,iBlock))-y1)**2/bb**2
!        !endif
!
!        !if (temp>0.0) then
!        !   eflux(iLon,iLat)=z1+cc*sqrt(temp)
!        !else
!        !   eflux(iLon,iLat)=z1
!        !endif
!
!        !eflux(iLon,iLat)=exp(log(10.)*eflux(iLon,iLat))
!        
!        !eflux(iLon,iLat)=eflux(iLon,iLat) * 1.602e-19 * 1e4 * 0.2
!!        eflux(iLon,iLat)=0.0  ! Assume not heat flux for electrons at Titan for now
!
!     enddo
!  enddo
!
! ! eflux = eflux*10.
! ! Lobe-like conditions Richard et al. [2011]
! ! electron fluxes = 
! ! approximate from Rymer et al. [2009].
! ! peak fluxes ~1.0e+05 e/cm^2/s
! ! 500 eV mean energy
!! eflux = (1.0e+05)*(1.0e+04)*(500.0)*Element_Charge ! J/m^2/s
! ! mlats = MLatitude(1:nLons,1:nLats,nAlts,iBlock)
! ! where(mlats .GE. 50. .AND. mlats .LE. -50.) eflux=1.e-5
!
! eflux(1:nLons,1:nLats) = (1.0e+05)*(1.0e+04)*(50.0)*Element_Charge ! J/m^2/s
!
!
!!sinI2 =  (B0(0:nLons+1,0:nLats+1,0:nAlts+1,iUp_,iBlock)/B0(0:nLons+1,0:nLats+1,0:nAlts+1,iMag_,iBlock))**2
!
!sinI2 =  1.0  ! Assume radial fields
!!where(sinI2 .LE. 0.01) sinI2=0.01
!sinI = sqrt(sinI2)
!
!te = eTemperature(1:nLons,1:nLats,0:nAlts+1,iBlock)
!ti = iTemperature(1:nLons,1:nLats,0:nAlts+1,iBlock) 
!
! !!!! Calculate electron temperature
!!lam_e = lame * sinI2
!lam_e = lame 
!
!JParaAlt = 0.0
!PartialJPara = 0.0
!
!!call calc_thermoelectric_current
!!UserData3D(:,:,:,2,iBlock) = 0.0
!!UserData3D(1:nLons,1:nLats,nAlts,2,iBlock) = JParaAlt(1:nLons,1:nLats)
!
!! Use tri-diagnal solver to solve the equation
!do iLon = 1, nLons
!   do iLat = 1, nLats 
!      
!      do iAlt = 1, nAlts 
!         zu = alts(iLon,iLat,iAlt+1) - alts(iLon,iLat,iAlt)
!         zl = alts(iLon,iLat,iAlt  ) - alts(iLon,iLat,iAlt-1)
!         
!         lamu = lam_e(iLon,iLat,iAlt+1) - lam_e(iLon,iLat,iAlt)
!         laml = lam_e(iLon,iLat,iAlt) - lam_e(iLon,iLat,iAlt-1)
!         nne  = ne(iLon,iLat,iAlt) 
!         tte  = te(iLon,iLat,iAlt)
!
!         neu = ne(iLon,iLat,iAlt+1) - ne(iLon,iLat,iAlt)
!         nel = ne(iLon,iLat,iAlt) - ne(iLon,iLat,iAlt-1)
!
!         uiu = uiup(iLon,iLat,iAlt+1) - uiup(iLon,iLat,iAlt)
!         uil = uiup(iLon,iLat,iAlt) - uiup(iLon,iLat,iAlt-1)
!
!         ilam = lam_e(iLon,iLat,iAlt)
!
!         hcoef = zu*zl*(zu+zl)  
!         xcoef = 1.5*Boltzmanns_Constant*nne/Dt   
!         fcoef = (zl/hcoef)**2*lamu + (zu/hcoef)**2*laml
!
!         b(iAlt) = -xcoef - 2*ilam/hcoef*(zu+zl) + fcoef*(zu**2-zl**2) - eHeatingp(iLon,iLat,iAlt)
!         a(iAlt) = 2*ilam/hcoef*zu - fcoef*zu**2
!         c(iAlt) = 2*ilam/hcoef*zl + fcoef*zl**2
!         d(iAlt) = -xcoef*tte - eHeatingm(iLon,iLat,iAlt)
!
!         eConduction(iLon,iLat,iAlt) = c(iAlt)*te(iLon,iLat,iAlt+1) + &
!              a(iAlt)*te(iLon,iLat,iAlt+1) + &
!              (-2*ilam/hcoef*(zu+zl) + fcoef*(zu**2-zl**2))*tte
!
!
!
!!         if (abs(MLatitude(iLon,iLat,nAlts,iBlock)) .GE. 45.) then
!!
!!            eThermo(iLon,iLat,iAlt) = &
!!                 5./2.*Boltzmanns_Constant &
!!                 /Element_Charge &
!!                 *JParaAlt(iLon,iLat) &
!!                 * ((te(iLon,iLat,iAlt+1)-te(iLon,iLat,iAlt))  *zl**2 &
!!                 +  (te(iLon,iLat,iAlt)  -te(iLon,iLat,iAlt-1))*zu**2)/hcoef
!!
!!            eHeatadv(iLon,iLat,iAlt) = eThermo(iLon,iLat,iAlt)/5.*3.
!!
!!
!!            eAdiab(iLon,iLat,iAlt) = &
!!                 - Boltzmanns_Constant * JParaAlt(iLon,iLat) &
!!                 /nne/Element_Charge * tte *&
!!                 (zl**2*neu+zu**2*nel)/hcoef
!!
!!
!!!!! use thermoelectric heating                                                
!!
!!            a(iAlt) = a(iAlt)-4*Boltzmanns_Constant &
!!                 /Element_Charge &
!!                 *JParaAlt(iLon,iLat)*zu**2/hcoef   
!!            !     *sinI(iLon,iLat,iAlt)
!!
!!            !     +1.5*nne*Boltzmanns_Constant* &
!!            !     uiup(iLon,iLat,iAlt) &
!!            !     *zu**2/hcoef
!!
!!            b(iAlt) = b(iAlt) + 4*Boltzmanns_Constant &
!!                 /Element_Charge * &
!!                 JParaAlt(iLon,iLat)*(zu**2-zl**2)/hcoef & 
!!            !     *sinI(iLon,iLat,iAlt) &
!!
!!                 - Boltzmanns_Constant * JParaAlt(iLon,iLat) &
!!                 /nne/Element_Charge * &
!!                 (zl**2*neu+zu**2*nel)/hcoef 
!!                 
!!            !     - Boltzmanns_Constant * nne * &
!!            !     (zl**2*uiu+zu**2*uil)/hcoef &
!!
!!            !     - 1.5*nne*Boltzmanns_Constant* &
!!            !     uiup(iLon,iLat,iAlt) * &
!!            !     (zu**2-zl**2)/hcoef
!!
!!
!!            c(iAlt) = c(iAlt)+4*Boltzmanns_Constant &
!!                 /Element_Charge &
!!                 *JParaAlt(iLon,iLat)*zl**2/hcoef 
!!            !     *sinI(iLon,iLat,iAlt)
!!                 
!!            !     -1.5*nne*Boltzmanns_Constant* &
!!            !     uiup(iLon,iLat,iAlt) &
!!            !     *zl**2/hcoef
!!                                   
!!         endif
!!
!      enddo
!
!     ! Bottom Bounday Te = Tn
!      a(1)=0
!      b(1)=1
!      c(1)=0
!      d(1)=Temperature(iLon,iLat,1,iBlock)*TempUnit(iLon,iLat,1)
!  
!      a(nAlts) =  -1
!      b(nAlts) =  1
!      c(nAlts) =  0
!      d(nAlts) =  dAlt_GB(iLon,iLat,nAlts,iBlock)*eflux(iLon,iLat) &
!           /lam_e(iLon,iLat,nAlts)
!
!      call tridag(a, b, c, d, u)
!
!      etemp(iLon,iLat,1:nAlts) = u(1:nAlts)
!      etemp(iLon,iLat,nAlts+1) = etemp(iLon,iLat,nAlts)
!      etemp(iLon,iLat,0) = etemp(iLon,iLat,1)
!       
!   enddo
!enddo
!
!lam_i = lami * sinI2
!
!!  do iAlt = 1, nAlts
!!    write(*,*) 'iAlt, Altitude, Lambda_i, Lambda_e =', &
!!               iAlt, Altitude_GB(1,1,iAlt,iBlock), lam_i(1,1,iAlt), lam_e(1,1,iAlt)
!!  enddo 
!!
!!stop
!! Use tri-diagnal solver to solve the equation
!
!do iLon = 1, nLons
!   do iLat = 1, nLats 
!      
!      do iAlt = 1, nAlts 
!         zu = alts(iLon,iLat,iAlt+1) - alts(iLon,iLat,iAlt)
!         zl = alts(iLon,iLat,iAlt  ) - alts(iLon,iLat,iAlt-1)
!         
!         lamu = lam_i(iLon,iLat,iAlt+1) - lam_i(iLon,iLat,iAlt)
!         laml = lam_i(iLon,iLat,iAlt) - lam_i(iLon,iLat,iAlt-1)
!         nni  = ni(iLon,iLat,iAlt) 
!         tti  = ti(iLon,iLat,iAlt)
!         ilam = lam_i(iLon,iLat,iAlt)
!         
!         hcoef = zu*zl*(zu+zl)  
!         xcoef = 1.5*Boltzmanns_Constant*nni/Dt    
!         fcoef = (zl/hcoef)**2*lamu + (zu/hcoef)**2*laml      
!
!         a(iAlt) = 2*ilam/hcoef*zu - fcoef*zu**2
!         b(iAlt) = -xcoef - 2*ilam/hcoef*(zu+zl) + fcoef*(zu**2-zl**2) - iHeatingp(iLon,iLat,iAlt)
!         c(iAlt) = 2*ilam/hcoef*zl + fcoef*zl**2
!         d(iAlt) = -xcoef*tti - iHeatingm(iLon,iLat,iAlt)
!
!         if (iAlt .EQ. nAlts) then
!            a(iAlt) = 0.
!            b(iAlt) = -xcoef - iHeatingp(iLon,iLat,iAlt)
!            c(iAlt) = 0.
!            d(iAlt) = -xcoef*tti - iHeatingm(iLon,iLat,iAlt)
!         endif
!         
!
!         iConduction(iLon,iLat,iAlt) = c(iAlt)*ti(iLon,iLat,iAlt+1) + &
!              a(iAlt)*ti(iLon,iLat,iAlt+1) + &
!              (-2*ilam/hcoef*(zu+zl) + fcoef*(zu**2-zl**2))*tti
!
!      enddo
!      
!      ! Bottom Bounday Ti = Tn
!      a(1)=0
!      b(1)=1
!      c(1)=0
!      d(1)=Temperature(iLon,iLat,1,iBlock)*TempUnit(iLon,iLat,1)
!      
!      call tridag(a, b, c, d, u)
!
!      itemp(iLon,iLat,1:nAlts) = u(1:nAlts) 
!      itemp(iLon,iLat,nAlts+1) = itemp(iLon,iLat,nAlts)
!      itemp(iLon,iLat,0) = itemp(iLon,iLat,1)
!   enddo
!enddo
!
!!!!! Check if reasonable temperatures
!where(etemp .LT. tn) etemp=tn
!where(etemp .GT. 6000.) etemp = 6000.

!eTemperature(1:nLons,1:nLats,0:nAlts+1,iBlock) = etemp(1:nLons,1:nLats,0:nAlts+1) 
!eTemperature(1:nLons,1:nLats,-1,iBlock) = etemp(1:nLons,1:nLats,0)
!eTemperature(1:nLons,1:nLats,nAlts+2,iBlock) = etemp(1:nLons,1:nLats,nAlts+1)
!
!where(itemp .GT. tipct * etemp) itemp = tipct * etemp
!where(itemp .LT. tn) itemp = tn
!where(itemp .GT. 6000.) itemp = 6000.
!
!iTemperature(1:nLons,1:nLats,0:nAlts+1,iBlock) = itemp(1:nLons,1:nLats,0:nAlts+1) 
!iTemperature(1:nLons,1:nLats,-1,iBlock) = itemp(1:nLons,1:nLats,0)
!iTemperature(1:nLons,1:nLats,nAlts+2,iBlock) = itemp(1:nLons,1:nLats,nAlts+1)

contains
  Subroutine tridag(a,b,c,d,u)

    use ModSizeGitm
    implicit none
    real(kind=8), intent(in)::a(nAlts),b(nAlts),c(nAlts),d(nAlts)
    real(kind=8), intent(out)::u(nAlts)
    real(kind=8):: cpp(nAlts), dpp(nAlts)
    real(kind=8) :: m
    
    ! initialize c-prime and d-prime
    cpp(1) = c(1)/b(1)
    dpp(1) = d(1)/b(1)
    ! solve for vectors c-prime and d-prime
    do iAlt = 2,nAlts
       m = b(iAlt)-cpp(iAlt-1)*a(iAlt)
       cpp(iAlt) = c(iAlt)/m
       dpp(iAlt) = (d(iAlt)-dpp(iAlt-1)*a(iAlt))/m
    enddo

    ! initialize u
    u(nAlts) = dpp(nAlts)
    ! solve for u from the vectors c-prime and d-prime
    do iAlt = nAlts-1, 1, -1
       u(iAlt) = dpp(iAlt)-cpp(iAlt)*u(iAlt+1)
    end do
    
  end subroutine tridag
   

!--------- calculate thermoelectric FAC currents-------                         

  Subroutine calc_thermoelectric_current

    real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2,3) :: iVelo, eVelo
    real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2,3) :: JuTotal, OverB0
    real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2) :: DivJPerp, JuTotalDotB
    real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: &
       Gradient_GC

    integer :: iDir

    iVelo = IVelocity(:,:,:,:,iBlock)
    eVelo = ExB(:,:,:,:)
    
    JuTotal = 0.0

    do iDir = 1, 3
       JuTotal(:,:,:,iDir) = IDensityS(:,:,:,ie_,iBlock)*Element_Charge* &
            (iVelo(:,:,:,iDir)-eVelo(:,:,:,iDir))

       OverB0(:,:,:,iDir) = B0(:,:,:,iDir,iBlock)/B0(:,:,:,iMag_,iBlock)**2
    enddo

     do iAlt = -1, nAlts+2
        do iLat = -1, nLats+2
           do iLon = -1, nLons+2
              JuTotalDotB(iLon, iLat, iAlt) = sum( &
                   JuTotal(iLon,iLat,iAlt,1:3)* &
                   B0(iLon,iLat,iAlt,1:3,iBlock))
           enddo
        enddo
     enddo


    DivJPerp = 0.0
    do iDir = 1, 3
       call UAM_Gradient_GC(JuTotal(:,:,:,iDir), Gradient_GC, iBlock)
       DivJPerp(:,:,:) = DivJPerp(:,:,:) + Gradient_GC(:,:,:,iDir)

       call UAM_Gradient_GC(JuTotalDotB(:,:,:), Gradient_GC, iBlock)
       DivJPerp(:,:,:) = DivJPerp(:,:,:) - Gradient_GC(:,:,:,iDir) &
            *B0(:,:,:,iDir,iBlock)/B0(:,:,:,iMag_,iBlock)**2

       call UAM_Gradient_GC(OverB0(:,:,:,iDir), Gradient_GC, iBlock)
       DivJPerp(:,:,:) = DivJPerp(:,:,:) - Gradient_GC(:,:,:,iDir) &
            *JuTotalDotB(:,:,:)
    enddo

    PartialJPara = - DivJPerp

    JParaAlt = 0.0
    do iAlt=1,nAlts
       JParaAlt(:,:) = JParaAlt(:,:) - DivJPerp(:,:,iAlt)*dAlt_GB(:,:,iAlt,iBlock)
    enddo

  end Subroutine calc_thermoelectric_current

end subroutine calc_electron_temperature

subroutine calc_electron_ion_sources(iBlock,eHeatingp,iHeatingp,eHeatingm,iHeatingm, iHeating, lame, lami)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModConstants
  use ModUserGITM
  use ModTime
  use ModInputs

  implicit none
  integer, intent(in) :: iBlock

  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1), intent(out) :: eHeatingp, iHeatingp, eHeatingm, &
       iHeatingm, iHeating, lame, lami
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: &
       nn, ni, ne, nn2, nch4, nh2, &
       tn, te, ti, temp, te_6000, te_exc, tn_exc, &
       nch5p, nc2h5p, nhcnhp, &
! Neutral-electron collision freqs
       nu_n2e, nu_ch4e, nu_h2e, &  ! Nu_en
! Ion-electron collision freqs
       nu_ch5pe, nu_c2h5pe, nu_hcnhpe, &  ! Nu_ie
! electron-electron collision freqs
       nu_ee, &  ! nu_ee
! electron-total collision freqs
       nu_etotal, &  ! nu_ee
!
! ion-neutral collision freqs
       nu_n2c2h5p, nu_n2ch5p, nu_n2hcnhp, & ! nu_in
       nu_h2c2h5p, nu_h2ch5p, nu_h2hcnhp, & ! nu_in
       nu_ch4c2h5p, nu_ch4ch5p, nu_ch4hcnhp, & !nu_in
!
! ion-ion self_collision freqs
       nu_ii_c2h5p, nu_ii_ch5p, nu_ii_hcnhp, & ! nu_ii
!
! ion-ion collision freqs
       nu_ch5pc2h5p, nu_ch5phcnhp, & ! nu_i-j, ion-ion collisions
       nu_c2h5phcnhp, &  ! nu_i-j, ion-ion collisions
       tr, dv2, dv2_en, dv2_ei 
  real(kind=8), dimension(-1:nLons+2,-1:nLats+2,0:nAlts+1) :: te_con, ti_con
  real(kind=8), dimension(-1:nLons+2,-1:nLats+2) :: dtedphe, dtidphe, dtedtheta, dtidtheta
  real(kind=8), dimension(-1:nLons+2,-1:nLats+2) :: dledphe, dlidphe, dledtheta, dlidtheta
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: Qphe, Qenc, Qeic, Qiec, &
       Qinc_t, Qinc_v, Qnic_t, Qnic_v, iAdvection, Qaurora, QprecipIon, &
       Qencp, Qeicp, Qiecp, Qinc_tp, Qrotp, Qfp, Qexcp, Qvib_o2p, Qvib_n2p, &
       Qencm, Qeicm, Qiecm, Qinc_tm, Qrotm, Qfm, Qexcm, Qvib_o2m, Qvib_n2m, &
       Qrot, Qf, Qeic_v, Qenc_v, Qexc, Qvib_o2, Qvib_n2, &
       Qeconhm, Qeconhp,Qiconhm, Qiconhp  !! Conduction perpendicular to field lines  
  !! Magnetic dip/declination angles, Conductivities
  real(kind=8), dimension(nLons,nLats,0:nAlts+1) :: cos2dip, magh2, sin2dec, cos2dec, sindec, cosdec 
  real(kind=8), dimension(nLons,nLats) :: sin2theta, sintheta, costheta !! polar angle
  real(kind=8) :: dLon, dLat
  real, dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: lam_op, lam_o2p, lam_n2p, lam_nop, lam_np
  real, dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: lam_c2h5p, lam_ch5p, lam_hcnhp
  integer :: Ao=16, Ao2=32, Ano=30, An=14
  real :: Ach4=16.0, Ah2=2.0, An2=28.0
  real :: Ach5p=17.0, Ac2h5p=29.0, Ahcnhp=28.0
  real:: Din_n2c2h5p, Din_n2ch5p, Din_n2hcnhp, &
          Din_h2c2h5p, Din_h2ch5p, Din_h2hcnhp, &
          Din_ch4c2h5p, Din_ch4ch5p, Din_ch4hcnhp

  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: x, epsilon, logx
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: dz_u, dz_l
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: alts 

! for O2 vibration
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1,10) :: Q0v
  real(kind=8), dimension(7) :: Av, Bv, Cv, Dv, Fv, Gv
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: logQ


! for V2 vibration
  real(kind=8), dimension(10) :: Av0, Bv0, Cv0, Dv0, Fv0
  real(kind=8), dimension(2:9) :: Av1, Bv1, Cv1, Dv1, Fv1
  real(kind=8) :: Av0L, Bv0L, Cv0L, Dv0L, Fv0L
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1,10)  :: logQv0
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1,2:9) :: logQv1 
  real(kind=8) :: tte, ttn, tte_6000
  integer :: iLevel


! for O fine structure
  real(kind=8), parameter :: s21 = 1.863e-11
  real(kind=8), parameter :: s20 = 1.191e-11  
  real, dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: s10, Dfine  

! for O excitation
  real(kind=8), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: dexc

! Aurora heaeting efficiency coefffient
  real(kind=8) :: auroheat=1.0

  integer :: iDir, iAlt, iLon, iLat
  real :: tipct = 1.1
  real :: mn, mi, e_esu, gamma_n, muin

  te  = eTemperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)
  ti  = iTemperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock) 
  te_con  = eTemperature(-1:nLons+2,-1:nLats+2,0:nAlts+1,iBlock)
  ti_con  = iTemperature(-1:nLons+2,-1:nLats+2,0:nAlts+1,iBlock) 
  tn  = Temperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)* &
        TempUnit(0:nLons+1,0:nLats+1,0:nAlts+1)
  nn  = NDensity(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)
  ne  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,ie_,iBlock)
  ni  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,ie_,iBlock)

  nh2 = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iH2_,iBlock)
  nn2 = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iN2_,iBlock)
 nch4 = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iCH4_,iBlock)
! Ions
  !nch5p  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iCH5P_,iBlock) 
  nch5p  = 1.0e+03
  nc2h5p  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iC2H5P_,iBlock) 
  nhcnhp  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iHCNHP_,iBlock) 


  alts = Altitude_GB(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)

 where(ti .LE. tn) ti = tn*1.0001
 where(te .LE. tn) te = tn*1.0001
 where(ti .GE. te*tipct) te=ti/tipct
 te_6000 = te
 te_exc = te
 where(te_6000 .GT. 6000.) te_6000 = 6000.0
 where(te_exc .GT. 18000.) te_exc = 18000.

! STEP1:  CALCULATE RELEVANT COLLISION FREQUENCIES
! Schunk and Nagy [2009], Matta et al. [2014]
! Electron Collision Frequencies
! Nu_ee, Nu_ei, Nu_en
      nu_ee = 54.5*(ne*1.0e-06)/sqrt(2.0)/te**1.5
  nu_ch5pe  = 54.5*(nch5p*1.0e-06)/Te**1.5
  nu_c2h5pe = 54.5*(nc2h5p*1.0e-06)/Te**1.5
  nu_hcnhpe = 54.5*(nhcnhp*1.0e-06)/Te**1.5
! Nu_en
! Schunk and Nagy [2009] -Table 4.6
  nu_n2e = (2.33e-11)*nn2*1.0e-06*&
           te*(1.0 - 1.21e-04*te)
 ! Rosenqvist et al. [2009], Kelley [1989]
 nu_ch4e = (5.4e-10)*nch4*1.0e-06*sqrt(te)
  nu_h2e = (5.4e-10)*nh2*1.0e-06*sqrt(te)
 
 nu_etotal = nu_ee + &
     (13.0/8.0)*(nu_ch5pe + nu_c2h5pe + nu_hcnhpe) + &
     ( 5.0/4.0)*(nu_n2e   + nu_ch4e   + nu_h2e   ) 

 
 ! lame is in units of [J/K/kg]*[J/m^3]/(1/s) -> s*(kg*m^2/s^2/K)*(J/m^3)/kg
 ! lame is then [m^2/s/K]*[J/m^3] -> J/m/K/s
 lame = (25.0/8.0)*(Boltzmanns_Constant/Mass_Electron)*&
       (te*ne*Boltzmanns_Constant)/nu_etotal 

! Next, CAlculate the Electron Conductivity
! lambda = (25/8)*kb*Pe/me/nu_etotal

  !lame = 7.7e5*te**2.5/(1+3.22e4*te**2/ne*nn*1.e-16)  
  ! lame Unit: eV cm-1 from Schunk and Nagy Page 147 eq 5.146
  !lame = lame *1.602e-19*100                      !Unit: J m-1

  Collision_ee(1:nLons,1:nLats,1:nAlts,iBlock) = &
         nu_ee(1:nLons,1:nLats,1:nAlts)

  Collision_ch5pe(1:nLons,1:nLats,1:nAlts,iBlock) = &
         nu_ch5pe(1:nLons,1:nLats,1:nAlts)

  Collision_c2h5pe(1:nLons,1:nLats,1:nAlts,iBlock) = &
         nu_c2h5pe(1:nLons,1:nLats,1:nAlts)

  Collision_hcnhpe(1:nLons,1:nLats,1:nAlts,iBlock) = &
         nu_hcnhpe(1:nLons,1:nLats,1:nAlts)

  Collision_n2e(1:nLons,1:nLats,1:nAlts,iBlock) = &
         nu_n2e(1:nLons,1:nLats,1:nAlts)

  Collision_ch4e(1:nLons,1:nLats,1:nAlts,iBlock) = &
         nu_ch4e(1:nLons,1:nLats,1:nAlts)

  Collision_h2e(1:nLons,1:nLats,1:nAlts,iBlock) = &
         nu_h2e(1:nLons,1:nLats,1:nAlts)

  Collision_Totale(1:nLons,1:nLats,1:nAlts,iBlock) = &
         nu_etotal(1:nLons,1:nLats,1:nAlts)

  eLambda_Temp(1:nLons,1:nLats,1:nAlts,iBlock) = &
         lame(1:nLons,1:nLats,1:nAlts)

! 1. Photoelectron heating from swartz and nisbet 1972
  Qphe=0.0
  !x = ne/(nn2+ no2 + no)
  x = ne/(nn2+ nh2 + nch4)
!!! fitting range
  where(x .GT. 10.) x=10.
  !  where(x .LT. 2.e-8) x=2.e-8

  ! This is from Smithro and Solomon [2008]
   logx=alog(x)
   epsilon = &
     exp(5.342 + 1.056*logx -4.392e-2*logx**2 -5.9e-2*logx**3 - &
        9.346e-3*logx**4 - 5.755e-4*logx**5 - 1.249e-5*logx**6) * &
       Element_Charge

   ! Qphe:  Parameterized Heating Rates (J/m^3/s)
   Qphe(1:nLons,1:nLats,1:nAlts) = epsilon(1:nLons,1:nLats,1:nAlts) * ( &
       N2TotalIonRateS(1:nLons,1:nLats,1:nAlts,iBlock) * &
                   nn2(1:nLons,1:nLats,1:nAlts)        + &
      CH4TotalIonRateS(1:nLons,1:nLats,1:nAlts,iBlock) * &
                  nch4(1:nLons,1:nLats,1:nAlts))

   PhotoEHeat(1:nLons,1:nLats,1:nAlts,iBlock) = Qphe(1:nLons,1:nLats,1:nAlts)

! Auroral heating
! No Auroral heating at Titan
   Qaurora = 0.0
! Ion Precipitation heating
! Ignore ion precipitation for now   
   QprecipIon = 0.0
! 2. ion-electron collisions From Schunk and Nagy 2009, and Bei-Chen Zhang and Y. Kamide 2003
!  Nu_ei = 54.5*Ni*Z^2.0/Te^1.5
! Qeic = (Nu_ei)*m_e*3.0*kb*(Ti - Tn)*(SUM(mu_ie))
  ! This is the ((3*kb*(Ti - te)) term
  !Qeic = ne*Mass_Electron*3.*Boltzmanns_Constant*(ti-te) * 5.45*1.e-5/(te**1.5) * &
  !      (nch5p/(Mass_Electron+MassI(iCH5P_))  &
  !     + nc2h5p/(Mass_Electron+MassI(iC2H5P_))   &
  !     + nhcnhp/(Mass_Electron+MassI(iHCNHP_)) )  
  ! Elastic electron-ion energy transfer
  ! J. Zhu:  Q_ei = 3.0*kb*ne*me*(sum(nu_ei)*(Ti-Te)/(me + mi)
  Qeic = 3.0*Boltzmanns_Constant*Mass_Electron*&
         ne*(ti-te) * (&
         nu_ch5pe/(Mass_Electron+MassI(iCH5P_)) + &
         nu_c2h5pe/(Mass_Electron+MassI(iC2H5P_)) + &
         nu_hcnhpe/(Mass_Electron+MassI(iHCNHP_)))

   Qei(1:nLons,1:nLats,1:nAlts,iBlock) = &
  Qeic(1:nLons,1:nLats,1:nAlts) 
  ! Units (J/K)*kg*(m^-3)*(K)*(Hz/kg) -> J/m^3/s


  ! Use this term for the (ue - ui)^2.0 joule heating term
!  Qeicp = ne*Mass_Electron*3.*Boltzmanns_Constant* 5.45*1.e-5/(te**1.5) * &
!        (nch5p/(Mass_Electron+MassI(iCH5P_))  &
!       + nc2h5p/(Mass_Electron+MassI(iC2H5P_))   &
!       + nhcnhp/(Mass_Electron+MassI(iHCNHP_)) )  
  Qeicp = 3.0*Boltzmanns_Constant*Mass_Electron*&
         ne*(&
         nu_ch5pe/(Mass_Electron+MassI(iCH5P_)) + &
         nu_c2h5pe/(Mass_Electron+MassI(iC2H5P_)) + &
         nu_hcnhpe/(Mass_Electron+MassI(iHCNHP_)))
  ! Units J/K/m^3/s
  Qeicm = Qeicp * ti
  Qiecp = Qeicp
  Qiecm = Qiecp * te

! 3. Elastic electron-neutral interactions From Shunk and Nagy 2009; actual value
! Approximate e-CH4 collisions with e-N2 for now
!  Qenc = ne*Mass_Electron*3.*Boltzmanns_Constant*(tn-te) *  &
!      ( 2.3276e-11 * nn2 * 1.e-6 * (1 - 1.21e-4 * te)* te / (Mass_Electron + Mass(iN2_)) + &
!        2.3276e-11 * nch4 * 1.e-6 * (1 - 1.21e-4 * te)* te / (Mass_Electron + Mass(iCH4_)) + &
!        5.2102e-10 * nh2 * 1.e-6 * sqrt(te)*(1.0 + 3.23e-2 * sqrt(te) - &
!          1.16e-06* (te**1.5) - 7.18e-9*(te**2.0)) / (Mass_Electron + Mass(iH2_))) 

  Qenc = 3.0*Boltzmanns_Constant*Mass_Electron*&
         ne*(tn-te) * (&
         nu_ch4e/(Mass_Electron+Mass(iCH4_)) + &
         nu_n2e/(Mass_Electron+Mass(iN2_)) + &
         nu_h2e/(Mass_Electron+Mass(iH2_)))


!  Qencp = ne*Mass_Electron*3.*Boltzmanns_Constant *  &
!      ( 2.3276e-11 * nn2 * 1.e-6 * (1 - 1.21e-4 * te)* te / (Mass_Electron + Mass(iN2_)) + &
!        2.3276e-11 * nch4 * 1.e-6 * (1 - 1.21e-4 * te)* te / (Mass_Electron + Mass(iCH4_)) + &
!        5.2102e-10 * nh2 * 1.e-6 * sqrt(te)*(1.0 + 3.23e-2 * sqrt(te) - &
!          1.16e-06* (te**1.5) - 7.18e-9*(te**2.0)) / (Mass_Electron + Mass(iH2_))) 

  Qencp = 3.0*Boltzmanns_Constant*Mass_Electron*&
         ne*(&
         nu_ch4e/(Mass_Electron+Mass(iCH4_)) + &
         nu_n2e/(Mass_Electron+Mass(iN2_)) + &
         nu_h2e/(Mass_Electron+Mass(iH2_)))

  Qencm = Qencp * tn

  Qen_elas(1:nLons,1:nLats,1:nAlts,iBlock) = &
     Qenc(1:nLons,1:nLats,1:nAlts)

!   real , Dimension(1:nLons,1:nLats,1:nAlts,nBlocksMax) :: Qen_inelas

   dv2_en = 0.0
!   do iDir = 1, 3
!      dv2_en = dv2_en + (Velocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock) &
!              -ExB(0:nLons+1,0:nLats+1,0:nAlts+1,iDir))**2
!   enddo

! Ignore the ExB part
!   do iDir = 1, 3
!      dv2_en = dv2_en + Velocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock) 
!   enddo
!
   ! Momentum coupling between electrons and Neutrals (U_e - U_n)
!   Qenc_v =ne*Mass_Electron*dv2_en *  &
!      ( 2.3276e-11 * nn2 * 1.e-6 * (1 - 1.21e-4 * te)* te / (Mass_Electron + Mass(iN2_)) + &
!        2.3276e-11 * nch4 * 1.e-6 * (1 - 1.21e-4 * te)* te / (Mass_Electron + Mass(iCH4_)) + &
!        5.2102e-10 * nh2 * 1.e-6 * sqrt(te)*(1.0 + 3.23e-2 * sqrt(te) - &
!          1.16e-06* (te**1.5) - 7.18e-9*(te**2.0)) / (Mass_Electron + Mass(iH2_))) 
   Qenc_v =0.0

   !dv2_ei = 0.0
   !do iDir = 1, 3
   !   dv2_ei = dv2_ei + (IVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock) &
   !        -ExB(0:nLons+1,0:nLats+1,0:nAlts+1,iDir))**2
   !enddo

   ! For Titan ignore the ExB stuff

!   dv2_ei = 0.0
!   do iDir = 1, 3
!      dv2_ei = dv2_ei + (IVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock)) 
!   enddo

  ! Old way
  !Qeic_v = ne*Mass_Electron*dv2_ei* 5.45*1.e-5/(te**1.5) * (nop+no2p+nn2p+nnop+nnp)

!  Qeic_v = ne*Mass_Electron*dv2_ei* 5.45*1.e-5/(te**1.5) * &
!        (nch5p + nc2h5p + nhcnhp)
  Qeic_v = 0.0


! ------
! Use only N2 rotation cooling
! 3. Neutral (N2, O2) rotation From Shunk & Nagy Page 277
  Qrot = 0.0
  !Qrot = 3.5e-14   * ne*1.e-6 * nn2*1.e-6 * (tn-te) / (te**0.5) 
  ! if Tn > Te, then you should heat, so Qrot > 0.0
  ! if Tn < Te, then Qrot < 0.0 --> Cooling
  Qrot = 3.5e-14   * ne*1.e-6 * nn2*1.e-6 * (tn-te)*sqrt(te)
  Qrot = Qrot * 1.6e-13     ! eV cm-3 s -> W m-3   

  !Qrotp =  3.5e-14   * ne*1.e-6 * nn2*1.e-6 / (te**0.5) 
  Qrotp = 3.5e-14   * ne*1.e-6 * nn2*1.e-6 *sqrt(te)
  Qrotp = Qrotp * 1.6e-13
  Qrotm = Qrotp * tn

  Qen_N2rot(1:nLons,1:nLats,1:nAlts,iBlock) = &
     Qrot(1:nLons,1:nLats,1:nAlts)


! 5. fine strcture heating rate  by Shunk and Nagy Page 282
! No fine structure contribution at Titan
  Qf    = 0.0
  Qf    = -Qf * 1.6e-13        ! eV cm-3 s-1 -> J m-3 s-1 
  Qfp = 0.
  Qfm = Qf

! 6. O(1D) excitation
  Qexc = 0.0
  Qexc = Qexc * 1.6e-13        ! eV cm-3 s-1 -> J m-3 s-1   

  Qexcp = 0.
  Qexcm = Qexc

! 7. O2 vibration
  Qvib_o2 = 0.0

! 8. N2 vibration from Pavlov 1998a
  Qvib_n2 =0.0
  !!! 1500K <= Te <=6000K
  Av0 = (/ 2.025, -7.066, -8.211, -9.713, -10.353, -10.819, -10.183, -12.698, -14.710, -17.538 /)
  Bv0 = (/ 8.782e-4, 1.001e-2, 1.092e-2, 1.204e-2, 1.243e-2, 1.244e-2, 1.185e-2, 1.309e-2, 1.409e-2, 1.6e-2 /)   
  Cv0 = (/ 2.954e-7, -3.066e-6, -3.369e-6, -3.732e-6, -3.850e-6, -3.771e-6, -3.570e-6, -3.952e-6, -4.249e-6, -4.916e-6 /)
  Dv0 = (/ -9.562e-11, 4.436e-10, 4.891e-10, 5.431e-10, 5.6e-10, 5.385e-10, 5.086e-10, 5.636e-10, 6.058e-10, 7.128e-10 /)
  Fv0 = (/ 7.252e-15, -2.449e-14, -2.706e-14, -3.008e-14, -3.1e-14, -2.936e-14, -2.769e-14, -3.071e-14, -3.3e-14, -3.941e-14 /)
 

  !!! 300K <= Te <= 1500K
  Av0L = -6.462
  Bv0L = 3.151e-2
  Cv0L = -4.075e-5
  Dv0L = 2.439e-8
  Fv0L = -5.479e-12

  !!! 1500K <= Te <=6000K
  Av1 = (/ -3.413, -4.16, -5.193, -5.939, -8.261, -8.185, -10.823, -11.273 /)
  Bv1 = (/ 7.326e-3, 7.803e-3, 8.36e-3, 8.807e-3, 1.01e-2, 1.01e-2, 1.199e-2, 1.283e-2 /)
  Cv1 = (/ -2.2e-6, -2.352e-6, -2.526e-6, -2.669e-6, -3.039e-6, -3.039e-6, -3.62e-6, -3.879e-6 /)
  Dv1 = (/ 3.128e-10, 3.352e-10, 3.606e-10, 3.806e-10, 4.318e-10, 4.318e-10, 5.159e-10, 5.534e-10 /)
  Fv1 = (/ -1.702e-14, -1.828e-14, -1.968e-14, -2.073e-14, -2.347e-14, -2.347e-14, -2.81e-14, -3.016e-14 /)
 
  !set small values for no n2 vibration
  logQv0 = -20.
  logQv1 = -20.

  do iLon=1, nLons
    do iLat=1, nLats
      do iAlt=1, nAlts 
      
        tte = te(iLon,iLat,iAlt)
        ttn = tn(iLon,iLat,iAlt)
        tte_6000 = te_6000(iLon,iLat,iAlt)
      
        if (tte .LE. 1500 .AND. tte .GT. 300) then
 !       if (tte .LE. 1500.) then
        
            logQv0(iLon,iLat,iAlt,1) = Av0L + Bv0L*tte + Cv0L*tte**2 + Dv0L*tte**3 + Fv0L*tte**4 -16.
            Qvib_n2(iLon,iLat,iAlt) = (1. - exp(-3353./ttn)) * &
                 10**logQv0(iLon,iLat,iAlt,1) * (1. - exp(1*3353.*(1./tte - 1./ttn)))
          
!         else if (tte .LE. 6000 .AND. tte .GT. 1500) then
         else if (tte .GT. 1500) then
           
           logQv0(iLon,iLat,iAlt,:) = Av0 + Bv0*tte_6000 + Cv0*tte_6000**2 + Dv0*tte_6000**3 + Fv0*tte_6000**4 - 16.
           logQv1(iLon,iLat,iAlt,:) = Av1 + Bv1*tte_6000 + Cv1*tte_6000**2 + Dv1*tte_6000**3 + Fv1*tte_6000**4 - 16.


           do iLevel = 1, 10
             Qvib_n2(iLon,iLat,iAlt) = Qvib_n2(iLon,iLat,iAlt) + &
               (1. - exp(-3353./ttn)) * 10**logQv0(iLon,iLat,iAlt,iLevel) * (1. - exp(iLevel*3353.*(1./tte - 1./ttn)))
           end do
 
           do iLevel = 2, 9
             Qvib_n2(iLon,iLat,iAlt) = Qvib_n2(iLon,iLat,iAlt) + &
               (1. - exp(-3353./ttn)) * exp(-3353./ttn) * 10**logQv1(iLon,iLat,iAlt,iLevel) * &
               (1. - exp((iLevel-1)*3353.*(1./tte-1./ttn)))
           end do         
   
        end if
       end do
    end do
  end do
  ! No N2-vibration below 300 K
 
  Qvib_n2  =  - ne * nn2 * 1.e-12 * Qvib_n2 * 1.6e-13

  Qvib_n2p = 0.
  Qvib_n2m = Qvib_n2

  Qen_N2vib(1:nLons,1:nLats,1:nAlts,iBlock) = &
    Qvib_n2(1:nLons,1:nLats,1:nAlts)


  Qen_inelas = Qen_N2vib + Qen_N2rot
!! END ELECTRON

!! \=========
!! BEGIN ION-----------------------------------------------------
!! /=========

! Need self-collision frequencies
  nu_ii_c2h5p = 1.27*nc2h5p*1.0e-06*sqrt(Ac2h5p/2.0)/Ac2h5p/ti**1.5
  nu_ii_ch5p  = 1.27*nch5p *1.0e-06*sqrt(Ach5p /2.0)/Ach5p /ti**1.5
  nu_ii_hcnhp = 1.27*nhcnhp*1.0e-06*sqrt(Ahcnhp/2.0)/Ahcnhp/ti**1.5

! Ion-Ion Collsion Frequencies

!
! Need Ion-Neutral momentum collision frequencies
! Ion Neutral Collision Frequencies
  tr = (tn+ti)/2.

  ! Nu-in = 2.21*pi*(mn/mn+mi)*Nn*sqrt(e^2*Gamma_n/muin)
  ! Nu_in in Schunk and Nagy, 4.8032045e-10 is "e" in esu units
  ! Gamma_N2 = 1.76*10^-24 (cm^3) polarizability

  !N2 - Ion
  e_esu = 4.8032045e-10  ! statvolts
  mn = Mass(iN2_)/AMU
  gamma_n = 1.76         ! polarizability
  mi = MassI(iC2H5P_)/AMU
  muin = mi*mn/(mn+mi)
!  nu_n2c2h5p = 2.21*PI*e_esu*(mn/(mi+mn))*nn2*1.0e-06*sqrt(gamma_n/muin)
  nu_n2c2h5p = 2.6e-09*nn2*1.0e-06*sqrt(gamma_n/muin)

  mi = MassI(iCH5P_)/AMU
  muin = mi*mn/(mn+mi)
!  nu_n2ch5p = 2.21*PI*e_esu*(mn/(mi+mn))*nn2*1.0e-06*sqrt(gamma_n/muin)
   nu_n2ch5p = 2.6e-09*nn2*1.0e-06*sqrt(gamma_n/muin)

  mi = MassI(iHCNHP_)/AMU
  muin = mi*mn/(mn+mi)
!  nu_n2hcnhp = 2.21*PI*e_esu*(mn/(mi+mn))*nn2*1.0e-06*sqrt(gamma_n/muin)
   nu_n2hcnhp = 2.6e-09*nn2*1.0e-06*sqrt(gamma_n/muin)

  ! H2-Ion
  mn = Mass(iH2_)/AMU
  gamma_n = 0.82         ! polarizability
  mi = MassI(iC2H5P_)/AMU
  muin = mi*mn/(mn+mi)
!  nu_h2c2h5p = 2.21*PI*e_esu*(mn/(mi+mn))*nn2*1.0e-06*sqrt(gamma_n/muin)
  nu_h2c2h5p = 2.6e-09*nh2*1.0e-06*sqrt(gamma_n/muin)

  mi = MassI(iCH5P_)/AMU
  muin = mi*mn/(mn+mi)
  nu_h2ch5p = 2.6e-09*nh2*1.0e-06*sqrt(gamma_n/muin)
!  nu_h2ch5p = 2.21*PI*e_esu*(mn/(mi+mn))*nn2*1.0e-06*sqrt(gamma_n/muin)

  mi = MassI(iHCNHP_)/AMU
  muin = mi*mn/(mn+mi)
  nu_h2hcnhp = 2.6e-09*nh2*1.0e-06*sqrt(gamma_n/muin)

  ! CH4-Ion
  mn = Mass(iCH4_)/AMU
  gamma_n = 2.59         ! polarizability
  mi = MassI(iC2H5P_)/AMU
  muin = mi*mn/(mn+mi)
  !nu_ch4c2h5p = 2.21*PI*e_esu*(mn/(mi+mn))*nn2*1.0e-06*sqrt(gamma_n/muin)
  nu_ch4c2h5p = 2.6e-09*nch4*1.0e-06*sqrt(gamma_n/muin)

  mi = MassI(iCH5P_)/AMU
  muin = mi*mn/(mn+mi)
  !nu_ch4ch5p = 2.21*PI*e_esu*(mn/(mi+mn))*nn2*1.0e-06*sqrt(gamma_n/muin)
  nu_ch4ch5p = 2.6e-09*nch4*1.0e-06*sqrt(gamma_n/muin)

  mi = MassI(iHCNHP_)/AMU
  muin = mi*mn/(mn+mi)
  !nu_ch4hcnhp = 2.21*PI*e_esu*(mn/(mi+mn))*nn2*1.0e-06*sqrt(gamma_n/muin)
  nu_ch4hcnhp = 2.6e-09*nch4*1.0e-06*sqrt(gamma_n/muin)

  dv2 = 0.0
  do iDir = 1, 3
     dv2 = dv2 + (IVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock) &
                  -Velocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock))**2
  enddo

   ! Next sum over the ion-neutral pairs
   ! Heating due to the Tn - Ti
   Qinc_t = nch5p*MassI(iCH5P_)*nu_n2ch5p*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iCH5P_) + Mass(iN2_)) &
          + nch5p*MassI(iCH5P_)*nu_h2ch5p*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iCH5P_) + Mass(iH2_)) &
          + nch5p*MassI(iCH5P_)*nu_ch4ch5p*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iCH5P_) + Mass(iCH4_)) &
!
          + nc2h5p*MassI(iC2H5P_)*nu_n2c2h5p*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iC2H5P_) + Mass(iN2_)) &
          + nc2h5p*MassI(iC2H5P_)*nu_h2c2h5p*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iC2H5P_) + Mass(iH2_)) &
          + nc2h5p*MassI(iC2H5P_)*nu_ch4c2h5p*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iC2H5P_) + Mass(iCH4_)) &
!
          + nhcnhp*MassI(iHCNHP_)*nu_n2hcnhp*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iHCNHP_) + Mass(iN2_)) &
          + nhcnhp*MassI(iHCNHP_)*nu_h2hcnhp*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iHCNHP_) + Mass(iH2_)) &
          + nhcnhp*MassI(iHCNHP_)*nu_ch4hcnhp*(3.*Boltzmanns_Constant*(tn-ti) )/&
               (MassI(iHCNHP_) + Mass(iCH4_)) 

   Qinc_tp = nch5p*MassI(iCH5P_)*nu_n2ch5p*(3.*Boltzmanns_Constant )/&
               (MassI(iCH5P_) + Mass(iN2_)) &
          + nch5p*MassI(iCH5P_)*nu_h2ch5p*(3.*Boltzmanns_Constant )/&
               (MassI(iCH5P_) + Mass(iH2_)) &
          + nch5p*MassI(iCH5P_)*nu_ch4ch5p*(3.*Boltzmanns_Constant )/&
               (MassI(iCH5P_) + Mass(iCH4_)) &
          + nc2h5p*MassI(iC2H5P_)*nu_n2c2h5p*(3.*Boltzmanns_Constant )/&
               (MassI(iC2H5P_) + Mass(iN2_)) &
          + nc2h5p*MassI(iC2H5P_)*nu_h2c2h5p*(3.*Boltzmanns_Constant )/&
               (MassI(iC2H5P_) + Mass(iH2_)) &
          + nc2h5p*MassI(iC2H5P_)*nu_ch4c2h5p*(3.*Boltzmanns_Constant )/&
               (MassI(iC2H5P_) + Mass(iCH4_)) &
          + nhcnhp*MassI(iHCNHP_)*nu_n2hcnhp*(3.*Boltzmanns_Constant )/&
               (MassI(iHCNHP_) + Mass(iN2_)) &
          + nhcnhp*MassI(iHCNHP_)*nu_h2hcnhp*(3.*Boltzmanns_Constant )/&
               (MassI(iHCNHP_) + Mass(iH2_)) &
          + nhcnhp*MassI(iHCNHP_)*nu_ch4hcnhp*(3.*Boltzmanns_Constant )/&
               (MassI(iHCNHP_) + Mass(iCH4_)) 

   Qinc_tm = Qinc_tp * tn

  ! For the Joule Heating Part
   Qinc_v = nch5p*MassI(iCH5P_)*nu_n2ch5p*dv2/(MassI(iCH5P_) + Mass(iN2_)) &
          + nch5p*MassI(iCH5P_)*nu_h2ch5p*dv2 /(MassI(iCH5P_) + Mass(iH2_)) &
          + nch5p*MassI(iCH5P_)*nu_ch4ch5p*dv2 /(MassI(iCH5P_) + Mass(iCH4_)) &
          + nc2h5p*MassI(iC2H5P_)*nu_n2c2h5p*dv2 /(MassI(iC2H5P_) + Mass(iN2_)) &
          + nc2h5p*MassI(iC2H5P_)*nu_h2c2h5p*dv2/(MassI(iC2H5P_) + Mass(iH2_)) &
          + nc2h5p*MassI(iC2H5P_)*nu_ch4c2h5p*dv2/(MassI(iC2H5P_) + Mass(iCH4_)) &
          + nhcnhp*MassI(iHCNHP_)*nu_n2hcnhp*dv2/(MassI(iHCNHP_) + Mass(iN2_)) &
          + nhcnhp*MassI(iHCNHP_)*nu_h2hcnhp*dv2/(MassI(iHCNHP_) + Mass(iH2_)) &
          + nhcnhp*MassI(iHCNHP_)*nu_ch4hcnhp*dv2/(MassI(iHCNHP_) + Mass(iCH4_)) 

!!!! Ion-neutral heating rate from Banks 1967
   Qnic_v = nch5p*MassI(iCH5P_)**2.*nu_n2ch5p*dv2/(MassI(iCH5P_) + Mass(iN2_)) &
          + nch5p*MassI(iCH5P_)**2.*nu_h2ch5p*dv2 /(MassI(iCH5P_) + Mass(iH2_)) &
          + nch5p*MassI(iCH5P_)**2.*nu_ch4ch5p*dv2 /(MassI(iCH5P_) + Mass(iCH4_)) &
!
          + nc2h5p*MassI(iC2H5P_)**2.*nu_n2c2h5p*dv2 /(MassI(iC2H5P_) + Mass(iN2_)) &
          + nc2h5p*MassI(iC2H5P_)**2.*nu_h2c2h5p*dv2/(MassI(iC2H5P_) + Mass(iH2_)) &
          + nc2h5p*MassI(iC2H5P_)**2.*nu_ch4c2h5p*dv2/(MassI(iC2H5P_) + Mass(iCH4_)) &
!
          + nhcnhp*MassI(iHCNHP_)**2.*nu_n2hcnhp*dv2/(MassI(iHCNHP_) + Mass(iN2_)) &
          + nhcnhp*MassI(iHCNHP_)**2.*nu_h2hcnhp*dv2/(MassI(iHCNHP_) + Mass(iH2_)) &
          + nhcnhp*MassI(iHCNHP_)**2.*nu_ch4hcnhp*dv2/(MassI(iHCNHP_) + Mass(iCH4_)) 

  Qnic_t = - Qinc_t

!!! Assume a radial field
!!! Thermal Conduction Perpendicular to Magnetic Field Lines

  !cos2dip =  1.- (B0(1:nLons,1:nLats,0:nAlts+1,iUp_,iBlock)/B0(1:nLons,1:nLats,0:nAlts+1,iMag_,iBlock))**2

!  magh2 =  B0(1:nLons,1:nLats,0:nAlts+1,iEast_,iBlock)**2 + B0(1:nLons,1:nLats,0:nAlts+1,iNorth_,iBlock)**2
!  cos2dec = B0(1:nLons,1:nLats,0:nAlts+1,iNorth_,iBlock)**2/magh2
!  sin2dec = B0(1:nLons,1:nLats,0:nAlts+1,iEast_,iBlock)**2/magh2  
!  sindec = B0(1:nLons,1:nLats,0:nAlts+1,iEast_,iBlock)/magh2**0.5
!  cosdec = B0(1:nLons,1:nLats,0:nAlts+1,iNorth_,iBlock)/magh2**0.5  

! Assume radial (so dip = Pi/2.0) sin dip = 1, cos2dip = 0.0
! Assume dec = 0.0, so cosdec = 1.0
  cos2dip = 0.0
  cos2dec = 1.0
  sin2dec = 0.0
  sindec  = 0.0
  cosdec  = 1.0

  dLat = Latitude(1,iBlock) - Latitude(0,iBlock)
  dLon = Longitude(1,iBlock) - Longitude(0,iBlock)

  ! Keep this conductivity for Titan as well
  ! Electron Conductivity 
  !lame = 7.7e5*te**2.5/(1+3.22e4*te**2/ne*nn*1.e-16)  
  ! lame Unit: eV cm-1 from Schunk and Nagy Page 147 eq 5.146
  !lame = lame *1.602e-19*100                      !Unit: J m-1
  
  Din_n2c2h5p  = ( 3.0*Ac2h5p**2.0 + An2**2.0  + (8.0/5.0)*Ac2h5p*An2)/(Ac2h5p  + An2)**2.0
  Din_n2ch5p   = ( 3.0* Ach5p**2.0 + An2**2.0  + (8.0/5.0)* Ach5p*An2)/(Ach5p   + An2)**2.0
  Din_n2hcnhp  = ( 3.0*Ahcnhp**2.0 + An2**2.0  + (8.0/5.0)*Ahcnhp*An2)/(Ahcnhp  + An2)**2.0
!
  Din_ch4c2h5p = ( 3.0*Ac2h5p**2.0 + Ach4**2.0 + (8.0/5.0)*Ac2h5p*Ach4)/(Ac2h5p + Ach4)**2.0
  Din_ch4ch5p  = ( 3.0* Ach5p**2.0 + Ach4**2.0 + (8.0/5.0)* Ach5p*Ach4)/(Ach5p  + Ach4)**2.0
  Din_ch4hcnhp = ( 3.0*Ahcnhp**2.0 + Ach4**2.0 + (8.0/5.0)*Ahcnhp*Ach4)/(Ahcnhp + Ach4)**2.0
!
  Din_h2c2h5p  = ( 3.0*Ac2h5p**2.0 + Ah2**2.0  + (8.0/5.0)*Ac2h5p*Ah2)/(Ac2h5p  + Ah2)**2.0
  Din_h2ch5p   = ( 3.0* Ach5p**2.0 + Ah2**2.0  + (8.0/5.0)* Ach5p*Ah2)/(Ach5p   + Ah2)**2.0
  Din_h2hcnhp  = ( 3.0*Ahcnhp**2.0 + Ah2**2.0  + (8.0/5.0)*Ahcnhp*Ah2)/(Ahcnhp  + Ah2)**2.0

  ! Ion Conductivity 
!  lam_ch5p = 3.1e4*ti**2.5/Ach5p**0.5 /(1+ 1.25* &
!       (nu_n2ch5p/nu_ii_ch5p)*(Din_n2ch5p + 1.5*(An2/(Ach5p + An2)))  +&
!       (nu_ch4ch5p/nu_ii_ch5p)*(Din_ch4ch5p + 1.5*(Ach4/(Ach5p + Ach4)))  + &
!       (nu_h2ch5p/nu_ii_ch5p)*(Din_h2ch5p + 1.5*(Ah2/(Ach5p + Ah2))) )


  ! Use ion -neutral collisions for now
  lam_ch5p = 3.1e4*ti**2.5/Ach5p**0.5 /(1+ 1.25* &
       ( (nu_n2ch5p/nu_ii_ch5p)*(Din_n2ch5p + 1.5*(An2/(Ach5p + An2)))  +&
        (nu_ch4ch5p/nu_ii_ch5p)*(Din_ch4ch5p + 1.5*(Ach4/(Ach5p + Ach4)))  + &
         (nu_h2ch5p/nu_ii_ch5p)*(Din_h2ch5p + 1.5*(Ah2/(Ach5p + Ah2)))) )

  lam_c2h5p = 3.1e4*ti**2.5/Ac2h5p**0.5 /(1+ 1.25*( &
       (nu_n2c2h5p/nu_ii_c2h5p)*(Din_n2c2h5p + 1.5*(An2/(Ac2h5p + An2)))  +&
       (nu_ch4c2h5p/nu_ii_c2h5p)*(Din_ch4c2h5p + 1.5*(Ach4/(Ac2h5p + Ach4)))  + &
       (nu_h2c2h5p/nu_ii_c2h5p)*(Din_h2c2h5p + 1.5*(Ah2/(Ac2h5p + Ah2))) ))

  lam_hcnhp = 3.1e4*ti**2.5/Ahcnhp**0.5 /(1+ 1.25* (&
       (nu_n2hcnhp/nu_ii_hcnhp)*(Din_n2hcnhp + 1.5*(An2/(Ahcnhp + An2)))  +&
       (nu_ch4hcnhp/nu_ii_hcnhp)*(Din_ch4hcnhp + 1.5*(Ach4/(Ahcnhp + Ach4)))  + &
       (nu_h2hcnhp/nu_ii_hcnhp)*(Din_h2hcnhp + 1.5*(Ah2/(Ahcnhp + Ah2))) ))


  lami = (nc2h5p*lam_c2h5p+nch5p*lam_ch5p+nhcnhp*lam_hcnhp) &
       /(nc2h5p+nch5p+nhcnhp)*1.602e-19*100                      !Unit: J m-1


  do iLat = 1, nLats
        sintheta(:, iLat) = sin(PI/2. - Latitude(iLat,iBlock))
        costheta(:, iLat) = cos(PI/2. - Latitude(iLat,iBlock))
  end do

  where (sintheta .LT. 0.1 .AND. sintheta .GE. 0.0) sintheta = 0.1
  where (sintheta .GT. -0.1 .AND. sintheta .LE. 0.0) sintheta = -0.1

  sin2theta = sintheta**2

  do iAlt = 0, nAlts+1
     dtedphe = 0.0
     dtidphe = 0.0
     dtedphe(0:nLons+1, 0:nLats+1) = &
            (te_con(1:nLons+2,0:nLats+1,iAlt) - te_con(-1:nLons,0:nLats+1,iAlt)) /2./dLon
     dtidphe(0:nLons+1, 0:nLats+1) = &
            (ti_con(1:nLons+2,0:nLats+1,iAlt) - ti_con(-1:nLons,0:nLats+1,iAlt)) /2./dLon

     dtedtheta = 0.0
     dtidtheta = 0.0
     dtedtheta(0:nLons+1, 0:nLats+1) = &
            (te_con(0:nLons+1,1:nLats+2,iAlt) - te_con(0:nLons+1,-1:nLats,iAlt)) /2./dLat
     dtidtheta(0:nLons+1, 0:nLats+1) = &
            (ti_con(0:nLons+1,1:nLats+2,iAlt) - ti_con(0:nLons+1,-1:nLats,iAlt)) /2./dLat

     dledphe = 0.0
     dlidphe = 0.0
     dledphe(1:nLons, 1:nLats) = &
             (lame(2:nLons+1,1:nLats,iAlt) - lame(0:nLons-1,1:nLats,iAlt)) /2./dLon
     dlidphe(1:nLons, 1:nLats) = &
             (lami(2:nLons+1,1:nLats,iAlt) - lami(0:nLons-1,1:nLats,iAlt)) /2./dLon

     dledtheta = 0.0
     dlidtheta = 0.0
     dledtheta(1:nLons, 1:nLats) = &
              (lame(1:nLons,2:nLats+1,iAlt) - lame(1:nLons,0:nLats-1,iAlt)) /2./dLat
     dlidtheta(1:nLons, 1:nLats) = &
              (lami(1:nLons,2:nLats+1,iAlt) - lami(1:nLons,0:nLats-1,iAlt)) /2./dLat

     Qeconhm(1:nLons,1:nLats,iAlt) = lame(1:nLons,1:nLats,iAlt) * cos2dip(1:nLons,1:nLats,iAlt) &
          / (6.37e6 + alts(1:nLons,1:nLats,iAlt))**2 * (cos2dec(1:nLons,1:nlats,iAlt) &
          * (te_con(1:nLons,2:nLats+1,iAlt) + te_con(1:nLons,0:nLats-1,iAlt)) / dLat**2 &
          + sin2dec(1:nLons,1:nLats,iAlt) / sin2theta  &
          * (te_con(2:nLons+1,1:nLats,iAlt) + te_con(0:nLons-1,1:nLats,iAlt)) / dLon**2 &
          + 2*sindec(1:nLons,1:nLats,iAlt)*cosdec(1:nLons,1:nLats,iAlt)/sintheta &
          * (dtedphe(1:nLons,2:nLats+1) - dtedphe(1:nLons,0:nLats-1))/2./dLat &
          - sindec(1:nLons,1:nLats,iAlt)*cosdec(1:nLons,1:nLats,iAlt)/sintheta**2*costheta * dtedphe(1:nLons,1:nLats) )&
          !!!!! Gradient in lambda !!!!!!
          + cos2dip(1:nLons,1:nLats,iAlt) /(6.37e6 + alts(1:nLons,1:nLats,iAlt))**2 * ( &
          cos2dec(1:nLons,1:nLats,iAlt) *dtedtheta(1:nLons,1:nLats)*dledtheta(1:nLons,1:nLats) + &
          + sin2dec(1:nLons,1:nLats,iAlt) / sin2theta * &
          dtedphe(1:nLons,1:nLats)*dledphe(1:nLons,1:nLats) + &
          + sindec(1:nLons,1:nLats,iAlt) * cosdec(1:nLons,1:nLats,iAlt) /sintheta * (&
          dtedtheta(1:nLons,1:nLats)*dledphe(1:nLons,1:nLats) + dtedtheta(1:nLons,1:nLats)*dledphe(1:nLons,1:nLats) ) )        
          

     Qeconhp(1:nLons,1:nLats,iAlt) = 2.*lame(1:nLons,1:nLats,iAlt) * cos2dip(1:nLons,1:nLats,iAlt) &
          / (6.37e6 + alts(1:nLons,1:nLats,iAlt))**2 * (cos2dec(1:nLons,1:nLats,iAlt)/ dLat**2 &
           + sin2dec(1:nLons,1:nLats,iAlt) / sin2theta / dLon**2 )


     Qiconhm(1:nLons,1:nLats,iAlt) = lami(1:nLons,1:nLats,iAlt) * cos2dip(1:nLons,1:nLats,iAlt) &
          / (6.37e6 + alts(1:nLons,1:nLats,iAlt))**2 * (cos2dec(1:nLons,1:nLats,iAlt) &
          * (ti_con(1:nLons,2:nLats+1,iAlt) + ti_con(1:nLons,0:nLats-1,iAlt)) / dLat**2 &
          + sin2dec(1:nLons,1:nLats,iAlt) / sin2theta  &
          * (ti_con(2:nLons+1,1:nLats,iAlt) + ti_con(0:nLons-1,1:nLats,iAlt)) / dLon**2 &
          + 2*sindec(1:nLons,1:nLats,iAlt)*cosdec(1:nLons,1:nLats,iAlt)/sintheta &
          * (dtidphe(1:nLons,2:nLats+1) - dtidphe(1:nLons,0:nLats-1))/2./dLat &
          - sindec(1:nLons,1:nLats,iAlt)*cosdec(1:nLons,1:nLats,iAlt)&
          /sintheta**2*costheta * dtidphe(1:nLons,1:nLats)) &
          !!!!! Gradient in lambda !!!!!!
          + cos2dip(1:nLons,1:nLats,iAlt) /(6.37e6 + alts(1:nLons,1:nLats,iAlt))**2 * ( &
          cos2dec(1:nLons,1:nLats,iAlt) *dtidtheta(1:nLons,1:nLats)*dlidtheta(1:nLons,1:nLats) + &
          + sin2dec(1:nLons,1:nLats,iAlt) / sin2theta * &
          dtidphe(1:nLons,1:nLats)*dlidphe(1:nLons,1:nLats) + &
          + sindec(1:nLons,1:nLats,iAlt) * cosdec(1:nLons,1:nLats,iAlt) /sintheta * (&
          dtidtheta(1:nLons,1:nLats)*dlidphe(1:nLons,1:nLats) + dtidtheta(1:nLons,1:nLats)*dlidphe(1:nLons,1:nLats) ) )        
          

     Qiconhp(1:nLons,1:nLats,iAlt) = 2.*lami(1:nLons,1:nLats,iAlt) * cos2dip(1:nLons,1:nLats,iAlt) &
          / (6.37e6 + alts(1:nLons,1:nLats,iAlt))**2 * (cos2dec(1:nLons,1:nLats,iAlt)/ dLat**2 &
           + sin2dec(1:nLons,1:nLats,iAlt) / sin2theta / dLon**2 )


     ! Qiconh(1:nLons,1:nLats,iAlt) = lami(1:nLons,1:nLats,iAlt) * cos2dip(1:nLons,1:nLats,iAlt) &
     ! / (6.37e6 + alts(1:nLons,1:nLats,iAlt))**2 * (cos2dec(1:nLons,1:nLats,iAlt) &
     !      * (ti_con(1:nLons,2:nLats+1,iAlt) + ti_con(1:nLons,0:nLats-1,iAlt) - 2*ti_con(1:nLons,1:nLats,iAlt)) / dLat**2 &
     !      + sin2dec(1:nLons,1:nLats,iAlt) / sin2theta  &
     !      * (ti_con(2:nLons+1,1:nLats,iAlt) + ti_con(0:nLons-1,1:nLats,iAlt) - 2*ti_con(1:nLons,1:nLats,iAlt)) / dLon**2 )
  end do

  
  if (UseJouleHeating .and. UseIonDrag) then

!     JouleHeating = (Qnic_t(1:nLons,1:nLats,1:nAlts) + Qnic_v(1:nLons,1:nLats,1:nAlts))/ &
!          TempUnit(1:nLons,1:nLats,1:nAlts) / &
!          cp(1:nLons,1:nLats,1:nAlts,iBlock) / Rho(1:nLons,1:nLats,1:nAlts,iBlock)
      
     JouleHeating = (Qnic_v(1:nLons,1:nLats,1:nAlts) * 2.)/ &
          TempUnit(1:nLons,1:nLats,1:nAlts) / &
          cp(1:nLons,1:nLats,1:nAlts,iBlock) / Rho(1:nLons,1:nLats,1:nAlts,iBlock)
    
 else

     JouleHeating = 0.0

 endif

  ElectronHeating = -Qenc(1:nLons,1:nLats,1:nAlts) / &
          TempUnit(1:nLons,1:nLats,1:nAlts) / &
          cp(1:nLons,1:nLats,1:nAlts,iBlock) / Rho(1:nLons,1:nLats,1:nAlts,iBlock)

!!!!!!  Qaurora = 0.0
! REmove the horizontal conduction
  Qeconhp = 0.
  Qiconhp = 0.
  Qeconhm = 0.
  Qiconhm = 0.


!  eHeatingm  = Qphe   ! This is in units of J/m^3/s

  !eHeatingm  = Qphe + Qencm + Qeicm + Qrotm + Qvib_n2 + Qenc_v + Qeic_v 
!  eHeatingm  = Qphe*1.0e-06 + Qencm + Qeicm + Qrotm + Qvib_n2 + Qenc_v + Qeic_v 
!  eHeatingp  = Qencp + Qeicp + Qrotp 

 ! eHeatingm  = 0.0
  eHeatingm  =  Qphe + Qencm + Qrotm + Qvib_n2  ! explicit coefficients
  eHeatingp  = Qencp + Qeicp+ Qrotp  ! coefficient for the implicit term

  iHeatingm  = Qiecm + Qinc_tm 
  iHeatingm(1:nLons,1:nLats,1:nAlts) = iHeatingm(1:nLons,1:nLats,1:nAlts) + &
                                       ChemicalHeatingRateIon
  iHeatingp  = Qiecp + Qinc_tp

!       + Qeconhm

!  eHeatingm  = Qphe + Qencm + Qeicm + Qrotm + Qf + Qexc + Qvib_o2 + Qvib_n2 + Qaurora + QprecipIon + Qenc_v + Qeic_v &
!       + Qeconhm
!  eHeatingp  = Qencp + Qeicp + Qrotp + Qeconhp
!  iHeatingm = Qiecm + Qinc_tm + Qinc_v + Qiconhm 
!  iHeatingm(1:nLons,1:nLats,1:nAlts) = iHeatingm(1:nLons,1:nLats,1:nAlts) + ChemicalHeatingRateIon
!  iHeatingp = Qiecp + Qinc_tp + Qiconhp
!  iHeating = Qiec+ Qinc_t + Qinc_v

end subroutine calc_electron_ion_sources
