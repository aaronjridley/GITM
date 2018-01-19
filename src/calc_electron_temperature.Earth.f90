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
  real(kind=8), intent(in), dimension(0:nLons+1,0:nLats+1,0:nAlts+1) :: &
       eHeatingp, iHeatingp, eHeatingm, &
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

!!!!   Upper Boundary Conditions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! [Reference: 1. J.T. Hastings and R.G. Roble, Planet. Space Sci., Vol.25, pp.209, 1977.
!!!!!! Equation(31)
!!!!!!   2. M.W. Liemohn et al, J. Geohys. Res., Vol.105, pp.27,767, 2000   Plate 1.

  Kp=2.0
  Dst=-Kp*20.0
  
  call time_real_to_int(CurrentTime, iTime)
  DoY = jday(iTime(1), &
       iTime(2), &
       iTime(3)) 
  
  bb1=2365-508*exp(-((DoY-jday(1998,6,12))/180.)**2)
  kk1=7.8+6.2*sin(DoY*pi/365.)
  P1=1200.+ 500.*(1+cos(DoY*4*pi/365.))
  P2=4200.+ 650.*(1+cos(DoY*4*pi/365.))
  bb2=64.+3.*sin(DoY*3*pi/365.)
  kk2=1.6+0.5*sin(DoY*3*pi/365.)
 
  x1=16.0
  y1=50
  z1=7.0
  aa=10.0
  bb=20.0
  cc=2.5

  do iLat=1,nLats
     do iLon=1,nLons
        
        if (MLT(iLon,iLat,nAlts) < 0.0) then
           MLT(iLon,iLat,nAlts) = MLT(iLon,iLat,nAlts) + 24.0
        endif
        
        if (abs(MLT(iLon,iLat,nAlts)-x1)<12) then
           temp=1-(MLT(iLon,iLat,nAlts)-x1)**2/aa**2- &
                (abs(MLatitude(iLon,iLat,nAlts,iBlock))-y1)**2/bb**2
        else
           temp=1-(24-abs(MLT(iLon,iLat,nAlts)-x1))**2/aa**2- &
                (abs(MLatitude(iLon,iLat,nAlts,iBlock))-y1)**2/bb**2
        endif

        if (temp>0.0) then
           eflux(iLon,iLat)=z1+cc*sqrt(temp)
        else
           eflux(iLon,iLat)=z1
        endif

        eflux(iLon,iLat)=exp(log(10.)*eflux(iLon,iLat))
        
        eflux(iLon,iLat)=eflux(iLon,iLat) * 1.602e-19 * 1e4 * 0.2

     enddo
  enddo

 ! eflux = eflux*10.
 ! eflux = 1.e-9
 ! mlats = MLatitude(1:nLons,1:nLats,nAlts,iBlock)
 ! where(mlats .GE. 50. .AND. mlats .LE. -50.) eflux=1.e-5


  sinI2 =  (B0(0:nLons+1,0:nLats+1,0:nAlts+1,iUp_,iBlock) / &
       B0(0:nLons+1,0:nLats+1,0:nAlts+1,iMag_,iBlock))**2
  where(sinI2 .LE. 0.01) sinI2=0.01
  sinI = sqrt(sinI2)

  te = eTemperature(1:nLons,1:nLats,0:nAlts+1,iBlock)
  ti = iTemperature(1:nLons,1:nLats,0:nAlts+1,iBlock) 

  ! Calculate electron temperature
    
  lam_e = lame * sinI2

  JParaAlt = 0.0
  PartialJPara = 0.0

  call calc_thermoelectric_current
  UserData3D(:,:,:,2,iBlock) = 0.0
  UserData3D(1:nLons,1:nLats,nAlts,2,iBlock) = JParaAlt(1:nLons,1:nLats)

  ! Use tri-diagnal solver to solve the equation
  do iLon = 1, nLons
     do iLat = 1, nLats 
      
        do iAlt = 1, nAlts 
           zu = alts(iLon,iLat,iAlt+1) - alts(iLon,iLat,iAlt)
           zl = alts(iLon,iLat,iAlt  ) - alts(iLon,iLat,iAlt-1)
         
           lamu = lam_e(iLon,iLat,iAlt+1) - lam_e(iLon,iLat,iAlt)
           laml = lam_e(iLon,iLat,iAlt) - lam_e(iLon,iLat,iAlt-1)
           nne  = ne(iLon,iLat,iAlt) 
           tte  = te(iLon,iLat,iAlt)

           neu = ne(iLon,iLat,iAlt+1) - ne(iLon,iLat,iAlt)
           nel = ne(iLon,iLat,iAlt) - ne(iLon,iLat,iAlt-1)

           uiu = uiup(iLon,iLat,iAlt+1) - uiup(iLon,iLat,iAlt)
           uil = uiup(iLon,iLat,iAlt) - uiup(iLon,iLat,iAlt-1)

           ilam = lam_e(iLon,iLat,iAlt)

           hcoef = zu*zl*(zu+zl)  
           xcoef = 1.5*Boltzmanns_Constant*nne/Dt   
           fcoef = (zl/hcoef)**2*lamu + (zu/hcoef)**2*laml

           c(iAlt) = 2*ilam/hcoef*zl + fcoef*zl**2
           b(iAlt) = -xcoef - 2*ilam/hcoef*(zu+zl) + fcoef*(zu**2-zl**2) - eHeatingp(iLon,iLat,iAlt)
           a(iAlt) = 2*ilam/hcoef*zu - fcoef*zu**2
           d(iAlt) = -xcoef*tte - eHeatingm(iLon,iLat,iAlt)

           eConduction(iLon,iLat,iAlt) = c(iAlt)*te(iLon,iLat,iAlt+1) + &
                a(iAlt)*te(iLon,iLat,iAlt+1) + &
                (-2*ilam/hcoef*(zu+zl) + fcoef*(zu**2-zl**2))*tte

           if (abs(MLatitude(iLon,iLat,nAlts,iBlock)) .GE. 45.) then

              eThermo(iLon,iLat,iAlt) = &
                   5./2.*Boltzmanns_Constant &
                   /Element_Charge &
                   *JParaAlt(iLon,iLat) &
                   * ((te(iLon,iLat,iAlt+1)-te(iLon,iLat,iAlt))  *zl**2 &
                   +  (te(iLon,iLat,iAlt)  -te(iLon,iLat,iAlt-1))*zu**2)/hcoef

              eHeatadv(iLon,iLat,iAlt) = eThermo(iLon,iLat,iAlt)/5.*3.

              eAdiab(iLon,iLat,iAlt) = &
                   - Boltzmanns_Constant * JParaAlt(iLon,iLat) &
                   /nne/Element_Charge * tte *&
                   (zl**2*neu+zu**2*nel)/hcoef
              
              ! use thermoelectric heating                                                

              a(iAlt) = a(iAlt)-4*Boltzmanns_Constant &
                   /Element_Charge &
                   *JParaAlt(iLon,iLat)*zu**2/hcoef   
              !     *sinI(iLon,iLat,iAlt)

              !     +1.5*nne*Boltzmanns_Constant* &
              !     uiup(iLon,iLat,iAlt) &
              !     *zu**2/hcoef

              b(iAlt) = b(iAlt) + 4*Boltzmanns_Constant &
                   /Element_Charge * &
                   JParaAlt(iLon,iLat)*(zu**2-zl**2)/hcoef & 
                   !     *sinI(iLon,iLat,iAlt) &

                   - Boltzmanns_Constant * JParaAlt(iLon,iLat) &
                   /nne/Element_Charge * &
                   (zl**2*neu+zu**2*nel)/hcoef 
                 
              !     - Boltzmanns_Constant * nne * &
              !     (zl**2*uiu+zu**2*uil)/hcoef &

              !     - 1.5*nne*Boltzmanns_Constant* &
              !     uiup(iLon,iLat,iAlt) * &
              !     (zu**2-zl**2)/hcoef

              c(iAlt) = c(iAlt)+4*Boltzmanns_Constant &
                   /Element_Charge &
                   *JParaAlt(iLon,iLat)*zl**2/hcoef 
              !     *sinI(iLon,iLat,iAlt)
                 
              !     -1.5*nne*Boltzmanns_Constant* &
              !     uiup(iLon,iLat,iAlt) &
              !     *zl**2/hcoef
                                   
           endif

        enddo

        ! Bottom Bounday Te = Tn
        a(1)=0
        b(1)=1
        c(1)=0
        d(1)=Temperature(iLon,iLat,1,iBlock)*TempUnit(iLon,iLat,1)
  
        a(nAlts) =  -1
        b(nAlts) =  1
        c(nAlts) =  0
        d(nAlts) =  dAlt_GB(iLon,iLat,nAlts,iBlock)*eflux(iLon,iLat) &
             /lam_e(iLon,iLat,nAlts)

        call tridag(a, b, c, d, u)

        etemp(iLon,iLat,1:nAlts) = u(1:nAlts)
        etemp(iLon,iLat,nAlts+1) = etemp(iLon,iLat,nAlts)
        etemp(iLon,iLat,0) = etemp(iLon,iLat,1)
       
     enddo
  enddo

  lam_i = lami * sinI2
  ! Use tri-diagnal solver to solve the equation

  do iLon = 1, nLons
     do iLat = 1, nLats 
      
        do iAlt = 1, nAlts 
           zu = alts(iLon,iLat,iAlt+1) - alts(iLon,iLat,iAlt)
           zl = alts(iLon,iLat,iAlt  ) - alts(iLon,iLat,iAlt-1)
           
           lamu = lam_i(iLon,iLat,iAlt+1) - lam_i(iLon,iLat,iAlt)
           laml = lam_i(iLon,iLat,iAlt) - lam_i(iLon,iLat,iAlt-1)
           nni  = ni(iLon,iLat,iAlt) 
           tti  = ti(iLon,iLat,iAlt)
           ilam = lam_i(iLon,iLat,iAlt)
         
           hcoef = zu*zl*(zu+zl)  
           xcoef = 1.5*Boltzmanns_Constant*nni/Dt    
           fcoef = (zl/hcoef)**2*lamu + (zu/hcoef)**2*laml      

           a(iAlt) = 2*ilam/hcoef*zu - fcoef*zu**2
           b(iAlt) = -xcoef - 2*ilam/hcoef*(zu+zl) + fcoef*(zu**2-zl**2) - iHeatingp(iLon,iLat,iAlt)
           c(iAlt) = 2*ilam/hcoef*zl + fcoef*zl**2
           d(iAlt) = -xcoef*tti - iHeatingm(iLon,iLat,iAlt)

           if (iAlt .EQ. nAlts) then
              a(iAlt) = 0.
              b(iAlt) = -xcoef - iHeatingp(iLon,iLat,iAlt)
              c(iAlt) = 0.
              d(iAlt) = -xcoef*tti - iHeatingm(iLon,iLat,iAlt)
           endif
         
           iConduction(iLon,iLat,iAlt) = c(iAlt)*ti(iLon,iLat,iAlt+1) + &
                a(iAlt)*ti(iLon,iLat,iAlt+1) + &
                (-2*ilam/hcoef*(zu+zl) + fcoef*(zu**2-zl**2))*tti
        enddo
      
        ! Bottom Bounday Ti = Tn
        a(1)=0
        b(1)=1
        c(1)=0
        d(1)=Temperature(iLon,iLat,1,iBlock)*TempUnit(iLon,iLat,1)
      
        call tridag(a, b, c, d, u)

        itemp(iLon,iLat,1:nAlts) = u(1:nAlts) 
        itemp(iLon,iLat,nAlts+1) = itemp(iLon,iLat,nAlts)
        itemp(iLon,iLat,0) = itemp(iLon,iLat,1)
     enddo
  enddo

  ! Check if reasonable temperatures
  where(etemp .LT. tn) etemp=tn
  where(etemp .GT. 6000.) etemp = 6000.

  eTemperature(1:nLons,1:nLats,0:nAlts+1,iBlock) = etemp(1:nLons,1:nLats,0:nAlts+1) 
  eTemperature(1:nLons,1:nLats,-1,iBlock) = etemp(1:nLons,1:nLats,0)
  eTemperature(1:nLons,1:nLats,nAlts+2,iBlock) = etemp(1:nLons,1:nLats,nAlts+1)

  where(itemp .GT. tipct * etemp) itemp = tipct * etemp
  where(itemp .LT. tn) itemp = tn
  where(itemp .GT. 6000.) itemp = 6000.

  iTemperature(1:nLons,1:nLats,0:nAlts+1,iBlock) = itemp(1:nLons,1:nLats,0:nAlts+1) 
  iTemperature(1:nLons,1:nLats,-1,iBlock) = itemp(1:nLons,1:nLats,0)
  iTemperature(1:nLons,1:nLats,nAlts+2,iBlock) = itemp(1:nLons,1:nLats,nAlts+1)

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
       nn, ni, ne, nh, nhe, no, nn2, no2, nnr, nno, &
       tn, te, ti, temp, te_6000, te_exc, tn_exc, &
       nop, no2p, nn2p, nnop, nhp, nhep, nnp, &
       nu_oop, nu_nnp, nu_o2o2p, nu_o2op, nu_n2op, nu_n2o2p, &
       nu_oo2p, nu_n2n2p, tr, dv2, dv2_en, dv2_ei 
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
  integer :: Ao=16, Ao2=32, An2=28, Ano=30, An=14

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

  te  = eTemperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)
  ti  = iTemperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock) 
  te_con  = eTemperature(-1:nLons+2,-1:nLats+2,0:nAlts+1,iBlock)
  ti_con  = iTemperature(-1:nLons+2,-1:nLats+2,0:nAlts+1,iBlock) 
  tn  = Temperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)* &
        TempUnit(0:nLons+1,0:nLats+1,0:nAlts+1)
  nn  = NDensity(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)
  ne  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,ie_,iBlock)
  ni  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,ie_,iBlock)
  no  = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iO_3P_,iBlock) &
       + NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iO_1D_,iBlock)

  !Neutral Nitrgen
  nnr = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iN_4S_,iBlock) &
       +NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iN_2D_,iBlock) &
       +NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iN_2P_,iBlock)

  nh  = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iH_,iBlock)
  nhe = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iHe_,iBlock)
  nn2 = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iN2_,iBlock)
  no2 = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iO2_,iBlock)
  nno = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iNO_,iBlock)
! Ions

  nop  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iO_4SP_,iBlock) + &
         IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iO_2DP_,iBlock) + &
         IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iO_2PP_,iBlock)
  no2p = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iO2P_,iBlock) 
  nn2p = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iN2P_,iBlock)
  nnop = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iNOP_,iBlock)
  nhp  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iHP_,iBlock)
  nhep = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iHeP_,iBlock)
  nnp  = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iNP_,iBlock)
 
  alts = Altitude_GB(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)

 where(ti .LE. tn) ti = tn*1.0001
 where(te .LE. tn) te = tn*1.0001
 where(ti .GE. te*tipct) te=ti/tipct
 te_6000 = te
 te_exc = te
 where(te_6000 .GT. 6000.) te_6000 = 6000.0
 where(te_exc .GT. 18000.) te_exc = 18000.

! 1. Photoelectron heating from swartz and nisbet 1972

  Qphe=0.0
  x = ne/(nn2+ no2 + no)
!!! fitting range
  where(x .GT. 10.) x=10.
  !  where(x .LT. 2.e-8) x=2.e-8

   logx=log(x)
   epsilon = exp(5.342 + 1.056*logx -4.392e-2*logx**2 -5.9e-2*logx**3 -9.346e-3*logx**4 &
        -5.755e-4*logx**5 - 1.249e-5*logx**6) * 1.6e-19

 !!! change epsilon  
!   where(alts .LE. 20000.) epsilon = 0.0001

   Qphe(1:nLons,1:nLats,1:nAlts) = epsilon(1:nLons,1:nLats,1:nAlts) * ( &
        EuvIonRateS(1:nLons,1:nLats,1:nAlts,iO_4SP_,iBlock)*no(1:nLons,1:nLats,1:nAlts)  &
        + EuvIonRateS(1:nLons,1:nLats,1:nAlts,iO2P_,iBlock)*no2(1:nLons,1:nLats,1:nAlts)  &
        + EuvIonRateS(1:nLons,1:nLats,1:nAlts,iN2P_,iBlock)*nn2(1:nLons,1:nLats,1:nAlts)  &
        + EuvIonRateS(1:nLons,1:nLats,1:nAlts,iNP_,iBlock)*nnr(1:nLons,1:nLats,1:nAlts)  &
        + EuvIonRateS(1:nLons,1:nLats,1:nAlts,iNOP_,iBlock)*nno(1:nLons,1:nLats,1:nAlts)  &
        + EuvIonRateS(1:nLons,1:nLats,1:nAlts,iO_2DP_,iBlock)*no(1:nLons,1:nLats,1:nAlts)  &
        + EuvIonRateS(1:nLons,1:nLats,1:nAlts,iO_2PP_,iBlock)*no(1:nLons,1:nLats,1:nAlts)  &
        ) 

! Auroral heating
   Qaurora = 0.0
   Qaurora(1:nLons,1:nLats,1:nAlts) = auroheat &
        * epsilon(1:nLons,1:nLats,1:nAlts) * ( &
        AuroralIonRateS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) &
        + AuroralIonRateS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock) &
        + AuroralIonRateS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock) &
   )
   
! Ion Precipitation heating
   
   QprecipIon = 0.0
   QprecipIon(1:nLons,1:nLats,1:nAlts) = epsilon(1:nLons,1:nLats,1:nAlts) * ( &
        IonPrecipIonRateS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) &
        + IonPrecipIonRateS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock) &
        + IonPrecipIonRateS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock) &
        )
        
! 2. ion-electron collisions From Schunk and Nagy 2009, and Bei-Chen Zhang and Y. Kamide 2003

  Qeic = ne*Mass_Electron*3.*Boltzmanns_Constant*(ti-te) * 5.45*1.e-5/(te**1.5) * &
        (nop/(Mass_Electron+MassI(iO_4SP_))  &
       + no2p/(Mass_Electron+MassI(iO2P_))   &
       + nn2p/(Mass_Electron+MassI(iN2P_))   &
       + nnop/(Mass_Electron+MassI(iNOP_))   &
  !    + nhp/(Mass_Electron+MassI(iHP_))   &
  !    + nhep/(Mass_Electron+MassI(iHeP_))   &
       + nnp/(Mass_Electron+MassI(iNP_)) &
       )


  Qeicp = ne*Mass_Electron*3.*Boltzmanns_Constant * 5.45*1.e-5/(te**1.5) * &
        (nop/(Mass_Electron+MassI(iO_4SP_))  &
       + no2p/(Mass_Electron+MassI(iO2P_))   &
       + nn2p/(Mass_Electron+MassI(iN2P_))   &
       + nnop/(Mass_Electron+MassI(iNOP_))   &
  !    + nhp/(Mass_Electron+MassI(iHP_))   &
  !    + nhep/(Mass_Electron+MassI(iHeP_))   &
       + nnp/(Mass_Electron+MassI(iNP_)) &
       )

  Qeicm = Qeicp * ti
  Qiecp = Qeicp
  Qiecm = Qiecp * te

! 3. Elastic electron-neutral interactions From Shunk and Nagy 2009; actual value
  Qenc = ne*Mass_Electron*3.*Boltzmanns_Constant*(tn-te) *  &
         ( 2.33e-11 * nn2 * 1.e-6 * (1 - 1.21e-4 * te)* te / (Mass_Electron + Mass(iN2_))         &
        + 1.82e-10 * no2 * 1.e-6 * (1 + 3.60e-2 * te**0.5) * te**0.5  / (Mass_Electron + Mass(iO2_))         &
        + 8.90e-11 * no  * 1.e-6 * (1 + 5.70e-4 * te) * te**0.5  / (Mass_Electron + Mass(iO_3P_))     &
 !      +  4.6e-10  * nhe * 1.e-6 * te**0.5  / (Mass_Electron + Mass(iH_)) + &
 !      +  4.5e-9   * nh  * 1.e-6 * (1 - 1.35e-4 * te) * te**0.5  / (Mass_Electron + Mass(iHe_))  &
         )


  Qencp = ne*Mass_Electron*3.*Boltzmanns_Constant *  &
         ( 2.33e-11 * nn2 * 1.e-6 * (1 - 1.21e-4 * te)* te / (Mass_Electron + Mass(iN2_))         &
        + 1.82e-10 * no2 * 1.e-6 * (1 + 3.60e-2 * te**0.5) * te**0.5  / (Mass_Electron + Mass(iO2_))         &
        + 8.90e-11 * no  * 1.e-6 * (1 + 5.70e-4 * te) * te**0.5  / (Mass_Electron + Mass(iO_3P_))     &
 !      +  4.6e-10  * nhe * 1.e-6 * te**0.5  / (Mass_Electron + Mass(iH_)) + &
 !      +  4.5e-9   * nh  * 1.e-6 * (1 - 1.35e-4 * te) * te**0.5  / (Mass_Electron + Mass(iHe_))  &
         )

  Qencm = Qencp * tn


   dv2_en = 0.0
   do iDir = 1, 3
      dv2_en = dv2_en + (Velocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock) &
              -ExB(0:nLons+1,0:nLats+1,0:nAlts+1,iDir))**2
   enddo


   Qenc_v = ne*Mass_Electron*dv2_en*  &
          ( 2.33e-11 * nn2 * 1.e-6 * (1 - 1.21e-4 * te)* te * Mass(iN2_)/ (Mass_Electron + Mass(iN2_))         &
         + 1.82e-10 * no2 * 1.e-6 * (1 + 3.60e-2 * te**0.5) * te**0.5 * Mass(iO2_) / (Mass_Electron + Mass(iO2_))         &
         + 8.90e-11 * no  * 1.e-6 * (1 + 5.70e-4 * te) * te**0.5  * Mass(iO2_) / (Mass_Electron + Mass(iO_3P_))     &
 ! !      +  4.6e-10  * nhe * 1.e-6 * te**0.5  / (Mass_Electron + Mass(iH_)) + &
 ! !      +  4.5e-9   * nh  * 1.e-6 * (1 - 1.35e-4 * te) * te**0.5  / (Mass_Electron + Mass(iHe_))  &
          )


   dv2_ei = 0.0
   do iDir = 1, 3
      dv2_ei = dv2_ei + (IVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock) &
           -ExB(0:nLons+1,0:nLats+1,0:nAlts+1,iDir))**2
   enddo

  Qeic_v = ne*Mass_Electron*dv2_ei* 5.45*1.e-5/(te**1.5) * (nop+no2p+nn2p+nnop+nnp)
!        (nop/(Mass_Electron+MassI(iO_4SP_))  &
!       + no2p/(Mass_Electron+MassI(iO2P_))   &
!       + nn2p/(Mass_Electron+MassI(iN2P_))   &
!       + nnop/(Mass_Electron+MassI(iNOP_))   &
!  !    + nhp/(Mass_Electron+MassI(iHP_))   &                                                                                                                          
!  !    + nhep/(Mass_Electron+MassI(iHeP_))   &                                                                                                                      !  
!       + nnp/(Mass_Electron+MassI(iNP_)) &
!       )




! 3. Neutral (N2, O2) rotation From Shunk & Nagy Page 277
  Qrot = 0.0
  Qrot = 3.5e-14   * ne*1.e-6 * nn2*1.e-6 * (tn-te) / (te**0.5) + & ! N2 
         5.2e-15   * ne*1.e-6 * no2*1.e-6 * (tn-te) / (te**0.5)     ! O2
  Qrot = Qrot * 1.6e-13     ! eV cm-3 s -> J m-3   


  Qrotp =  3.5e-14   * ne*1.e-6 * nn2*1.e-6 / (te**0.5) + & ! N2 
         5.2e-15   * ne*1.e-6 * no2*1.e-6 / (te**0.5)     ! O2
  Qrotp = Qrotp * 1.6e-13
  Qrotm = Qrotp * tn


! 5. fine strcture heating rate  by Shunk and Nagy Page 282
  Qf    = 0.0
  Dfine = 5. + exp(-326.6/tn) + 3.*exp(-227.7/tn)
  s10   = 8.249e-16 * te**0.6 * exp(-227.7/tn)  
  Qf    = ne * no * 1.e-12 / Dfine * (s10*(1.-exp( 98.9*(1/te-1/tn))) + &
                                      s20*(1.-exp(326.6*(1/te-1/tn))) + &
                                      s21*(1.-exp(227.7*(1/te-1/tn))))
  Qf    = -Qf * 1.6e-13        ! eV cm-3 s-1 -> J m-3 s-1 

  Qfp = 0.
  Qfm = Qf

! 6. O(1D) excitation

  Qexc = 0.0
  where(tn .GE. te_exc) tn = te_exc
  dexc = 2.4e4 + 0.3*(te_exc-1500.) - 1.947e-5*(te_exc-1500.)*(te_exc-4000.)
  Qexc = 1.57e-12*ne*no*1.e-12* exp(dexc* (te_exc-3000.)/3000./te_exc) * (exp(-22713*(te_exc-tn)/te_exc/tn)-1)
  Qexc = Qexc * 1.6e-13        ! eV cm-3 s-1 -> J m-3 s-1   

  Qexcp = 0.
  Qexcm = Qexc



! 7. O2 vibration

  Qvib_o2 = 0.0

  ! set a small value for no O2 vibration
  logQ=-20.

  where(te_6000 .GE. 300) logQ = 5.0148e-31*te_6000**9 - 1.5346e-26*te_6000**8 &
       + 2.0127e-22*te_6000**7 - 1.4791e-18*te_6000**6 + 6.6865e-15*te_6000**5 &
       - 1.9228e-11*te_6000**4 + 3.5187e-8*te_6000**3 - 3.996e-5*te_6000**2 &
          + 0.0267*te_6000 - 19.9171
    
  Qvib_o2 = ne * no2 * 1.e-12 * 10**(logQ) * (1 - exp(2239.*(1./te - 1./tn)))
  Qvib_o2 = - Qvib_o2 * 1.6e-13

  Qvib_o2p = 0.
  Qvib_o2m = Qvib_o2

 

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
 
  Qvib_n2  =  - ne * nn2 * 1.e-12 * Qvib_n2 * 1.6e-13

  Qvib_n2p = 0.
  Qvib_n2m = Qvib_n2

!!!!! Ion-neutral heating rate
  tr = (tn+ti)/2.
  nu_nnp   = 0.0
  nu_oop   = 0.0
  nu_o2o2p = 0.0
  nu_n2n2p = 0.0
  nu_o2op  = 0.0
  nu_n2op  = 0.0
  nu_n2o2p = 0.0
  nu_oo2p  = 0.0
  where(tr .GT. 275.) nu_nnp = 3.83e-11 * nnr*1.e-6 * tr**0.5 * (1 - 0.063*log10(tr))**2
  where(tr <= 275) nu_nnp = 1.0e-9 * nnr*1.0e-6 ! T < 550

  where(tr .GT. 235.) nu_oop = 3.67e-11 * no*1.e-6 * tr**0.5 * (1 - 0.064*log10(tr))**2
  where(tr <= 235) nu_oop = 8.6e-10 * no*1.0e-6 ! T < 470

  where(tr .GT. 800.) nu_o2o2p = 2.59e-11 * no2*1.e-6 * tr**0.5 * (1 - 0.073*log10(tr))**2
  where(tr <= 800) nu_o2o2p = 8.2e-10 * no2*1.0e-6  ! T < 1600

  !where(tr .GT. 170.) nu_n2n2p = 5.14e-11 * nn2*1.e-6 * tr**0.5 * (1 - 0.069*log10(tr))**2
  nu_n2n2p = 5.14e-11 * nn2*1.e-6 * tr**0.5 * (1 - 0.069*log10(tr))**2
  ! Banks:

  nu_o2op = 6.64 * 1.e-10 * no2*1.e-6
  nu_n2op = 6.82*1.e-10 * nn2*1.e-6
  nu_n2o2p = 4.13*1.e-10 * nn2*1.e-6
  nu_oo2p = 2.31*1.e-10 * no*1.e-6
  
  dv2 = 0.0
  do iDir = 1, 3
     dv2 = dv2 + (IVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock) &
             -Velocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock))**2
  enddo
  
  Qinc_t = nop*MassI(iO_4SP_)*nu_oop * (3.*Boltzmanns_Constant*(tn-ti) )/(MassI(iO_4SP_) + Mass(iO_3P_)) &
       + nnp*MassI(iNP_)*nu_nnp * (3.*Boltzmanns_Constant*(tn-ti) )/(MassI(iNP_) + Mass(iN_2P_)) &
       + nop*MassI(iO_4SP_)*nu_o2op * (3.*Boltzmanns_Constant*(tn-ti) )/(MassI(iO_4SP_) + Mass(iO2_)) &
       + nop*MassI(iO_4SP_)*nu_n2op * (3.*Boltzmanns_Constant*(tn-ti) )/(MassI(iO_4SP_) + Mass(iN2_)) &
       + no2p*MassI(iO2P_)*nu_o2o2p * (3.*Boltzmanns_Constant*(tn-ti))/(MassI(iO2P_) + Mass(iO2_)) &
       + no2p*MassI(iO2P_)*nu_n2o2p * (3.*Boltzmanns_Constant*(tn-ti))/(MassI(iO2P_) + Mass(iN2_)) &
       + nn2p*MassI(iN2P_) * nu_n2n2p *(3.*Boltzmanns_Constant*(tn-ti))/(MassI(iN2P_) + Mass(iN2_)) &
       + no2p*MassI(iO2P_)*nu_oo2p * (3.*Boltzmanns_Constant*(tn-ti))/(MassI(iO2P_) + Mass(iO_3P_)) 

  Qinc_tp = nop*MassI(iO_4SP_)*nu_oop * 3.*Boltzmanns_Constant/(MassI(iO_4SP_) + Mass(iO_3P_)) &
       + nnp*MassI(iNP_)*nu_nnp * 3.*Boltzmanns_Constant/(MassI(iNP_) + Mass(iN_2P_)) &
       + nop*MassI(iO_4SP_)*nu_o2op * 3.*Boltzmanns_Constant/(MassI(iO_4SP_) + Mass(iO2_)) &
       + nop*MassI(iO_4SP_)*nu_n2op * 3.*Boltzmanns_Constant/(MassI(iO_4SP_) + Mass(iN2_)) &
       + no2p*MassI(iO2P_)*nu_o2o2p * 3.*Boltzmanns_Constant/(MassI(iO2P_) + Mass(iO2_)) &
       + no2p*MassI(iO2P_)*nu_n2o2p * 3.*Boltzmanns_Constant/(MassI(iO2P_) + Mass(iN2_)) &
       + nn2p*MassI(iN2P_) * nu_n2n2p *3.*Boltzmanns_Constant/(MassI(iN2P_) + Mass(iN2_)) &
       + no2p*MassI(iO2P_)*nu_oo2p * 3.*Boltzmanns_Constant/(MassI(iO2P_) + Mass(iO_3P_)) 

  Qinc_tm = Qinc_tp * tn

  Qinc_v = nop*MassI(iO_4SP_)*nu_oop * Mass(iO_3P_)*dv2/(MassI(iO_4SP_) + Mass(iO_3P_)) &
       + nnp*MassI(iNP_)*nu_nnp * Mass(iN_2P_)*dv2/(MassI(iNP_) + Mass(iN_2P_)) &
       + nop*MassI(iO_4SP_)*nu_o2op * Mass(iO2_)*dv2/(MassI(iO_4SP_) + Mass(iO2_)) &
       + nop*MassI(iO_4SP_)*nu_n2op * Mass(iN2_)*dv2/(MassI(iO_4SP_) + Mass(iN2_)) &
       + no2p*MassI(iO2P_)*nu_o2o2p * Mass(iO2_)*dv2/(MassI(iO2P_) + Mass(iO2_)) &
       + no2p*MassI(iO2P_)*nu_n2o2p * Mass(iN2_)*dv2/(MassI(iO2P_) + Mass(iN2_)) &
       + nn2p*MassI(iN2P_) * nu_n2n2p * Mass(iN2_)*dv2/(MassI(iN2P_) + Mass(iN2_))&
       + no2p*MassI(iO2P_)*nu_oo2p * Mass(iO_3P_)*dv2/(MassI(iO2P_) + Mass(iO_3P_)) 

  ! Ion-neutral heating rate from Banks 1967
  ! because ns * ms * nu_st = nt * mt * nu_ts, the only thing that changes is the neutral mass to ion mass:

  Qnic_v = nop * MassI(iO_4SP_)**2 *nu_oop *dv2/(MassI(iO_4SP_) + Mass(iO_3P_)) &
       + nnp * MassI(iNP_)**2 *nu_nnp *dv2/(MassI(iNP_) + Mass(iN_2P_)) &
       + nop * MassI(iO_4SP_)**2 *nu_o2op *dv2/(MassI(iO_4SP_) + Mass(iO2_)) &
       + nop *MassI(iO_4SP_)**2 *nu_n2op *dv2/(MassI(iO_4SP_) + Mass(iN2_)) &
       + no2p *MassI(iO2P_)**2 *nu_o2o2p *dv2/(MassI(iO2P_) + Mass(iO2_)) &
       + no2p *MassI(iO2P_)**2 *nu_n2o2p *dv2/(MassI(iO2P_) + Mass(iN2_)) &
       + nn2p *MassI(iN2P_)**2 * nu_n2n2p *dv2/(MassI(iN2P_) + Mass(iN2_)) &
       + no2p *MassI(iO2P_)**2 *nu_oo2p *dv2/(MassI(iO2P_) + Mass(iO_3P_)) 

  ! because ns * ms * nu_st = nt * mt * nu_ts, this is symmetric:
  Qnic_t = - Qinc_t
 
!!! Thermal Conduction Perpendicular to Magnetic Field Lines
   cos2dip =  1.- (B0(1:nLons,1:nLats,0:nAlts+1,iUp_,iBlock)/B0(1:nLons,1:nLats,0:nAlts+1,iMag_,iBlock))**2
   magh2 =  B0(1:nLons,1:nLats,0:nAlts+1,iEast_,iBlock)**2 + B0(1:nLons,1:nLats,0:nAlts+1,iNorth_,iBlock)**2
   cos2dec = B0(1:nLons,1:nLats,0:nAlts+1,iNorth_,iBlock)**2/magh2
   sin2dec = B0(1:nLons,1:nLats,0:nAlts+1,iEast_,iBlock)**2/magh2  
   sindec = B0(1:nLons,1:nLats,0:nAlts+1,iEast_,iBlock)/magh2**0.5
   cosdec = B0(1:nLons,1:nLats,0:nAlts+1,iNorth_,iBlock)/magh2**0.5  
   dLat = Latitude(1,iBlock) - Latitude(0,iBlock)
   dLon = Longitude(1,iBlock) - Longitude(0,iBlock)

  ! Electron Conductivity 
  lame = 7.7e5*te**2.5/(1+3.22e4*te**2/ne*nn*1.e-16)  !Unit: eV cm-1 from Schunk and Nagy Page 147 eq 5.146
  lame = lame *1.602e-19*100                      !Unit: J m-1
  
  ! Ion Conductivity 
  lam_op = 3.1e4*ti**2.5/Ao**0.5 /(1+ 1.75* &
       (no2p/nop*(Ao2/(Ao2+Ao))**0.5*(3*Ao**2+1.6*Ao2*Ao+1.3*Ao2**2)/(Ao+Ao2)**2 &
       +nn2p/nop*(An2/(An2+Ao))**0.5*(3*Ao**2+1.6*An2*Ao+1.3*An2**2)/(Ao+An2)**2 &
       +nnop/nop*(Ano/(Ano+Ao))**0.5*(3*Ao**2+1.6*Ano*Ao+1.3*Ano**2)/(Ao+Ano)**2 &
       +nnp /nop*(An/(An+Ao))**0.5*(3*Ao**2+1.6*An*Ao+1.3*An**2)/(Ao+An)**2 &
       ) )      !Unit: eV cm-1 s-1 K-1 from Shunk P153 Assume O atom is dominant, and ignore ...
  
  lam_o2p = 3.1e4*ti**2.5/Ao2**0.5 /(1+ 1.75* &
       (nop/no2p*(Ao/(Ao2+Ao))**0.5*(3*Ao2**2+1.6*Ao2*Ao+1.3*Ao**2)/(Ao+Ao2)**2 &
       +nn2p/no2p*(An2/(An2+Ao2))**0.5*(3*Ao2**2+1.6*An2*Ao2+1.3*An2**2)/(Ao2+An2)**2 &
       +nnop/no2p*(Ano/(Ano+Ao2))**0.5*(3*Ao2**2+1.6*Ano*Ao2+1.3*Ano**2)/(Ao2+Ano)**2 &
       +nnp/no2p*(An/(An+Ao2))**0.5*(3*Ao2**2+1.6*An*Ao2+1.3*An**2)/(Ao2+An)**2 &
       ))
  
  lam_n2p = 3.1e4*ti**2.5/An2**0.5 /(1+ 1.75* &
       (nop/nn2p*(Ao/(An2+Ao))**0.5*(3*An2**2+1.6*An2*Ao+1.3*Ao**2)/(Ao+An2)**2 &
       +no2p/nn2p*(Ao2/(An2+Ao2))**0.5*(3*An2**2+1.6*An2*Ao2+1.3*Ao2**2)/(Ao2+An2)**2 &
       +nnop/nn2p*(Ano/(An2+Ano))**0.5*(3*An2**2+1.6*An2*Ano+1.3*Ano**2)/(Ano+An2)**2 &
       +nnp/nn2p*(An/(An2+An))**0.5*(3*An2**2+1.6*An2*An+1.3*An**2)/(An+An2)**2 &
       ) )

  lam_nop = 3.1e4*ti**2.5/Ano**0.5 /(1+ 1.75* &
       (nop/nnop*(Ao/(Ano+Ao))**0.5*(3*Ano**2+1.6*Ano*Ao+1.3*Ao**2)/(Ao+Ano)**2 &
       +no2p/nnop*(Ao2/(Ano+Ao2))**0.5*(3*Ano**2+1.6*Ano*Ao2+1.3*Ao2**2)/(Ao2+Ano)**2 &
       +nn2p/nnop*(An2/(Ano+An2))**0.5*(3*Ano**2+1.6*Ano*An2+1.3*An2**2)/(An2+Ano)**2 &
       +nnp/nnop*(An/(Ano+An))**0.5*(3*Ano**2+1.6*Ano*An+1.3*An**2)/(An+Ano)**2 &
       ) )
  
  lam_np = 3.1e4*ti**2.5/An**0.5 /(1+ 1.75* &
       (nop/nnp*(Ao/(An+Ao))**0.5*(3*An**2+1.6*An*Ao+1.3*Ao**2)/(Ao+An)**2 &
       +no2p/nnp*(Ao2/(An+Ao2))**0.5*(3*An**2+1.6*An*Ao2+1.3*Ao2**2)/(Ao2+An)**2 &
       +nn2p/nnp*(An2/(An+An2))**0.5*(3*An**2+1.6*An*An2+1.3*An2**2)/(An2+An)**2 &
       +nnop/nnp*(Ano/(An+Ano))**0.5*(3*An**2+1.6*An*Ano+1.3*Ano**2)/(Ano+An)**2 &
       ) )
  lami = (nop*lam_op+no2p*lam_o2p+nnop*lam_nop+nn2p*lam_n2p+nnp*lam_np) &
       /(nop+no2p+nnop+nn2p+nnp)*1.602e-19*100                      !Unit: J m-1

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
     dtedphe(0:nLons+1, 0:nLats+1) = (te_con(1:nLons+2,0:nLats+1,iAlt) - te_con(-1:nLons,0:nLats+1,iAlt)) /2./dLon
     dtidphe(0:nLons+1, 0:nLats+1) = (ti_con(1:nLons+2,0:nLats+1,iAlt) - ti_con(-1:nLons,0:nLats+1,iAlt)) /2./dLon

     dtedtheta = 0.0
     dtidtheta = 0.0
     dtedtheta(0:nLons+1, 0:nLats+1) = (te_con(0:nLons+1,1:nLats+2,iAlt) - te_con(0:nLons+1,-1:nLats,iAlt)) /2./dLat
     dtidtheta(0:nLons+1, 0:nLats+1) = (ti_con(0:nLons+1,1:nLats+2,iAlt) - ti_con(0:nLons+1,-1:nLats,iAlt)) /2./dLat

     dledphe = 0.0
     dlidphe = 0.0
     dledphe(1:nLons, 1:nLats) = (lame(2:nLons+1,1:nLats,iAlt) - lame(0:nLons-1,1:nLats,iAlt)) /2./dLon
     dlidphe(1:nLons, 1:nLats) = (lami(2:nLons+1,1:nLats,iAlt) - lami(0:nLons-1,1:nLats,iAlt)) /2./dLon

     dledtheta = 0.0
     dlidtheta = 0.0
     dledtheta(1:nLons, 1:nLats) = (lame(1:nLons,2:nLats+1,iAlt) - lame(1:nLons,0:nLats-1,iAlt)) /2./dLat
     dlidtheta(1:nLons, 1:nLats) = (lami(1:nLons,2:nLats+1,iAlt) - lami(1:nLons,0:nLats-1,iAlt)) /2./dLat

     ! Qeconh(:,:,iAlt) = lame(:,:,iAlt) * cos2dip(:,:,iAlt) / (6.37e6 + alts(:,:,iAlt))**2 * (cos2dec(:,:,iAlt) &
     !      * (te_con(1:nLons,2:nLats+1,iAlt) + te_con(1:nLons,0:nLats-1,iAlt) - 2*te_con(1:nLons,1:nLats,iAlt)) / dLat**2 &
     !      + sin2dec(:,:,iAlt) / sin2theta  &
     !      * (te_con(2:nLons+1,1:nLats,iAlt) + te_con(0:nLons-1,1:nLats,iAlt) - 2*te_con(1:nLons,1:nLats,iAlt)) / dLon**2 &
     !      + 2*sindec(:,:,iAlt)*cosdec(:,:,iAlt)/sintheta * (dtedphe(1:nLons,2:nLats+1) - dtedphe(1:nLons,0:nLats-1))/2./dLat &
     !      - sindec(:,:,iAlt)*cosdec(:,:,iAlt)/sintheta**2*costheta * dtedphe(1:nLons,1:nLats) &
     !      )

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

  
  JouleHeating2d = 0.0
  HeatTransfer2d = 0.0

  if (UseJouleHeating .and. UseIonDrag) then

     JouleHeating = (Qnic_t(1:nLons,1:nLats,1:nAlts) + Qnic_v(1:nLons,1:nLats,1:nAlts))/ &
          TempUnit(1:nLons,1:nLats,1:nAlts) / &
          cp(1:nLons,1:nLats,1:nAlts,iBlock) / Rho(1:nLons,1:nLats,1:nAlts,iBlock)

     ! This is the old approximation: 
     !     JouleHeating = (Qnic_v(1:nLons,1:nLats,1:nAlts) * 2.)/ &
     !          TempUnit(1:nLons,1:nLats,1:nAlts) / &
     !          cp(1:nLons,1:nLats,1:nAlts,iBlock) / Rho(1:nLons,1:nLats,1:nAlts,iBlock)

     do iAlt=1,nAlts
        JouleHeating2d(1:nLons, 1:nLats) = &
             JouleHeating2d(1:nLons, 1:nLats) + &
             Qnic_v(1:nLons,1:nLats,iAlt) * dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
        HeatTransfer2d(1:nLons, 1:nLats) = &
             HeatTransfer2d(1:nLons, 1:nLats) + &
             Qnic_t(1:nLons,1:nLats,iAlt) * dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
     enddo
    
 else

     JouleHeating = 0.0

  endif

  ElectronHeating = -Qenc(1:nLons,1:nLats,1:nAlts) / &
          TempUnit(1:nLons,1:nLats,1:nAlts) / &
          cp(1:nLons,1:nLats,1:nAlts,iBlock) / Rho(1:nLons,1:nLats,1:nAlts,iBlock)

!!!!!!  Qaurora = 0.0
!  Qeconhp = 0.
!  Qiconhp = 0.
!  Qeconhm = 0.
!  Qiconhm = 0.

  eHeatingm  = Qphe + Qencm + Qeicm + Qrotm + Qf + Qexc + Qvib_o2 + Qvib_n2 + Qaurora + QprecipIon + Qenc_v + Qeic_v &
       + Qeconhm
  eHeatingp  = Qencp + Qeicp + Qrotp + Qeconhp
  iHeatingm = Qiecm + Qinc_tm + Qinc_v + Qiconhm 
  iHeatingm(1:nLons,1:nLats,1:nAlts) = iHeatingm(1:nLons,1:nLats,1:nAlts) + ChemicalHeatingRateIon
  iHeatingp = Qiecp + Qinc_tp + Qiconhp
  iHeating = Qiec+ Qinc_t + Qinc_v

end subroutine calc_electron_ion_sources
