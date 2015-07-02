!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_etemp_sources(Heating,Cooling,iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources, only: eEuvHeating
  use ModConstants

  implicit none

  integer, external :: jday

  integer, intent(in) :: iBlock
  real, intent(out)::Heating(nLons, nLats, nAlts)
  real, intent(out)::Cooling(nLons, nLats, nAlts)


  integer :: i,j,k,N, CLAWIter

  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens

  real :: R,NA,ddalt,c,h,kbc,omega
  real :: c1,c2,c3,E1,E2,E3,AA1,AA2,AA3,B1,B2,B3 
  

  integer, dimension(7) :: iTime

  integer :: nSteps, isteps
  real :: DtSub

  integer :: iError,DoY

  real, dimension(nLons,                 &
                  nLats,                  &
                  1:nAlts) :: TOld,        &
                   SolarHeating,eJouleHeating,expansion,        &
                  f,g,hh,T0,T1,  &
                  ZZ,Dx1,Dx2,Dx3,Ex1,Ex2,Ex3,ABT,d,Ki2,Ke2,   &
                  lnA,Lrot_E_N2,Lrot_E_O2,Lvib_e_N2,         &
                  Lvib_e_O2,Lf_e_O, L_e_O1D, L_e_i,Temp_T
                                   

real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2, &
                  nSpeciesTotal) :: Temp_NDensity

real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2, &
                  nIons) :: Temp_IDensity

!!!!!!!!!!!!!!!!!!!!!!!!!!

!  call report("Heating Terms",1)

  c=2.998e8
  h=6.626e-34
  kbc=1.381e-23
  NA=6.0221e23
  R= 8.3145 

Temp_NDensity(:,:,:,:)=NDensityS(:,:,:,:,iBlock)*1e-6
Temp_IDensity(:,:,:,:)=IDensityS(:,:,:,:,iBlock)*1e-6

  
!!$ Ke(i,j,k,iBlock) = 7.7e5 * eTemperature(i,j,k,iBlock)**2.5 

!!!!!! change the units eV cm-1 s-1 deg-1 to Joule m-1 s-1 deg-1 !!!!!!!!!1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i=1,nLons
   do j=1,nLats
      do k = 1, nAlts
!     write(*,*)'It is Ok----------1' 


  Lrot_E_N2(i,j,k)=2.9e-14 * Temp_IDensity(i,j,k,ie_) *       &
                             Temp_NDensity(i,j,k,iN2_) *               &
                            (eTemperature(i,j,k,iBlock) -        &
                             (Temperature(i,j,k,iBlock)*TempUnit(i,j,k)))/         &
                            sqrt(eTemperature(i,j,k,iBlock))

  Lrot_E_O2(i,j,k) =6.9e-14 * Temp_IDensity(i,j,k,ie_) *       &
                              Temp_NDensity(i,j,k,iO2_) *               &
                             (eTemperature(i,j,k,iBlock) -        &
                              (Temperature(i,j,k,iBlock)*TempUnit(i,j,k)))/         &
                             sqrt(eTemperature(i,j,k,iBlock))


  f(i,j,k) = 1.06e4 + 7.51e3 * tanh(1.10e-3 *                     &
      (eTemperature(i,j,k,iBlock)-1800.))

  g(i,j,k) = 3300. + 1.233 * (eTemperature(i,j,k,iBlock) - 1000.) -   &
           2.056e-4 * (eTemperature(i,j,k,iBlock) - 1000.) *     &
                      (eTemperature(i,j,k,iBlock) - 4000.)
        
  hh(i,j,k) = 3300. - 839. * sin(1.91e-4 *                       &
            (eTemperature(i,j,k,iBlock)-2700.))

     
  Lvib_e_N2(i,j,k) = -2.99e-12 *Temp_IDensity(i,j,k,ie_) *            &
                               Temp_NDensity(i,j,k,iN2_) *                    &
                    exp(f(i,j,k)*(eTemperature(i,j,k,iBlock) - 2000)/   & 
                               (2000.*eTemperature(i,j,k,iBlock)))*     &
                   (exp(-g(i,j,k)*(eTemperature(i,j,k,iBlock)-          &
                                   Temperature(i,j,k,iBlock)*TempUnit(i,j,k))/ &
                                  (eTemperature(i,j,k,iBlock)*          &
                                   Temperature(i,j,k,iBlock)*TempUnit(i,j,k)))-1)
 

!     write(*,*)'It is Ok----------1bb'

  Lvib_e_O2(i,j,k) = -5.196e-13 *Temp_IDensity(i,j,k,ie_) *           &
                               Temp_NDensity(i,j,k,iO2_) *                    &
                    exp(hh(i,j,k)*(eTemperature(i,j,k,iBlock) - 700)/    & 
                               (700.*eTemperature(i,j,k,iBlock)))*      &
                          (exp(-2770*(eTemperature(i,j,k,iBlock)-       &
                                   (Temperature(i,j,k,iBlock)*TempUnit(i,j,k)))/ &
                                  (eTemperature(i,j,k,iBlock)*          &
                                   (Temperature(i,j,k,iBlock)*TempUnit(i,j,k))))-1)                              
!     write(*,*)'It is Ok----------1b'
  T0(i,j,k) =  Temperature(i,j,k,iBlock)*TempUnit(i,j,k)
  T1(i,j,k) =  Temperature(i,j,k,iBlock)*TempUnit(i,j,k)
  ZZ(i,j,k) = 5. + 3.* exp(-228./T1(i,j,k)) + exp(-326./T0(i,j,k))
  c1 = 0.02
  c2 = 0.028
  c3 = 0.008
  Dx1(i,j,k) = exp(-228./T1(i,j,k))
  Dx2(i,j,k) = exp(-326./T0(i,j,k))
  Dx3(i,j,k) = exp(-326./T0(i,j,k))
  Ex1(i,j,k) = exp(-228./eTemperature(i,j,k,iBlock))
  Ex2(i,j,k) = exp(-326./eTemperature(i,j,k,iBlock))
  Ex3(i,j,k) = exp(-98./eTemperature(i,j,k,iBlock)-228./T1(i,j,k))
  E1 = 228.
  E2 = 326.
  E3 = 98.
  AA1 = 8.58e-6
  AA2 = 7.201e-6
  AA3 = 2.463e-7
  B1 = 1.008
  B2 = 0.9617
  B3 = 1.1448


  ABT(i,j,k) =AA1*B1*eTemperature(i,j,k,iBlock)**(B1-0.5)*          &
             ( c1*(Dx1(i,j,k)-Ex1(i,j,k))+5.91e-9*(TempUnit(i,j,k)* &
             Temperature(i,j,k,iBlock)-     &
             eTemperature(i,j,k,iBlock))*((1+B1)*Dx1(i,j,k)+               &
             (E1/eTemperature(i,j,k,iBlock)+1+B1)*Ex1(i,j,k)))

  ABT(i,j,k) =ABT(i,j,k)+AA2*B2*eTemperature(i,j,k,iBlock)**(B2-0.5)*&
             ( c2*(Dx2(i,j,k)-Ex2(i,j,k))+5.91e-9*(TempUnit(i,j,k)* &
             Temperature(i,j,k,iBlock)-     &
             eTemperature(i,j,k,iBlock))*((1+B2)*Dx2(i,j,k)+               &
             (E2/eTemperature(i,j,k,iBlock)+1+B2)*Ex2(i,j,k)))

  ABT(i,j,k) =ABT(i,j,k)+AA3*B3*eTemperature(i,j,k,iBlock)**(B3-0.5)*&
             ( c3*(Dx3(i,j,k)-Ex3(i,j,k))+5.91e-9*(TempUnit(i,j,k)*Temperature(i,j,k,iBlock)-     &
             eTemperature(i,j,k,iBlock))*((1+B3)*Dx3(i,j,k)+               &
             (E3/eTemperature(i,j,k,iBlock)+1+B3)*Ex3(i,j,k)))


  Lf_e_O(i,j,k) = -8.629e-6 * Temp_IDensity(i,j,k,ie_) *           &
                             Temp_NDensity(i,j,k,iO_3P_) *            &
                             ABT(i,j,k)/ZZ(i,j,k)

  Lf_e_O(i,j,k) = 3.4e-12*(1-7e-5*eTemperature(i,j,k,iBlock))*    &
                  Temp_IDensity(i,j,k,ie_) *           &
                  Temp_NDensity(i,j,k,iO_3P_) *            &
                  (150./eTemperature(i,j,k,iBlock)+0.4)* &
                  (eTemperature(i,j,k,iBlock)- &
                   Temperature(i,j,k,iBlock)*TempUnit(i,j,k))/ &
                   (Temperature(i,j,k,iBlock)*TempUnit(i,j,k))


!     write(*,*)'It is Ok----------1c'


  d(i,j,k) = 2.4e4 + 0.3 * (eTemperature(i,j,k,iBlock)-1500.) -      &
             1.947e-5 * (eTemperature(i,j,k,iBlock)-1500.) *         &
                        (eTemperature(i,j,k,iBlock)-4000.)

  L_e_O1D(i,j,k) = -1.57e-12 *Temp_IDensity(i,j,k,ie_) *              &
                               Temp_NDensity(i,j,k,iO_3P_) *                     &
                    exp(d(i,j,k)*(eTemperature(i,j,k,iBlock) -3000)/    & 
                               (3000.*eTemperature(i,j,k,iBlock)))*     &
                          (exp(-22713*(eTemperature(i,j,k,iBlock)-      &
                                   TempUnit(i,j,k)*Temperature(i,j,k,iBlock))/     &
                                  (eTemperature(i,j,k,iBlock)*          &
                                   TempUnit(i,j,k)*Temperature(i,j,k,iBlock)))-1)
!     write(*,*)'It is Ok----------1d'

  Ki2(i,j,k) = 4 * 3.14159 * Temp_IDensity(i,j,k,ie_) *              &
               (1.602e-19)**2/(kbc*ITemperature(i,j,k,iBlock))
  Ke2(i,j,k) = 4 * 3.14159 * Temp_IDensity(i,j,k,ie_) *              &
               (1.602e-19)**2/(kbc*eTemperature(i,j,k,iBlock))

lnA(i,j,k) = LOG(4. * kbc * eTemperature(i,j,k,iBlock) /                &
                 ((1.602e-19)**2 * sqrt(Ke2(i,j,k)))) - 2 * 0.577 -     &
                 (Ki2(i,j,k) + Ke2(i,j,k)) / Ki2(i,j,k) *               &
                 LOG(sqrt((Ki2(i,j,k) + Ke2(i,j,k)) / Ke2(i,j,k)))

lnA(i,j,k) = 15

  L_e_i(i,j,k) = 3.2e-8/3. * Temp_IDensity(i,j,k,ie_) *                 &
                         (eTemperature(i,j,k,iBlock) -                  &
                          iTemperature(i,j,k,iBlock)) /                      &
                          (eTemperature(i,j,k,iBlock))**1.5*            &
                          lnA(i,j,k) * (Temp_IDensity(i,j,k,iO_4SP_)+         &
!                          4.* Temp_IDensity(i,j,k,iHeP_)+                    &
!                          16.*Temp_IDensity(i,j,k,iHP_) +                    &
                          0.5 * Temp_IDensity(i,j,k,iO2P_) +                 &
                          0.53 * Temp_IDensity(i,j,k,iNOP_))
!     write(*,*)'It is Ok----------1e'

  TOld(i,j,k) = eTemperature(i,j,k,iBlock)

!     write(*,*)'It is Ok----------2'           

        
!!!!  Heating !!!!!!!!!!!!!!!!!

           SolarHeating(i,j,k)= eEuvHeating(i,j,k,iBlock)       
           eJouleHeating(i,j,k)=0.0
           heating(i,j,k) = SolarHeating(i,j,k) + &
                            eJouleHeating(i,j,k)
           heating(i,j,k) = heating(i,j,k) * 0.5
!     write(*,*)'It is Ok----------3'

!!!!  cooling !!!!!!!!!!!!!!!!


           cooling(i,j,k) = Lvib_E_O2(i,j,k)             &
                           +Lvib_E_N2(i,j,k)            &
                           + L_e_O1D(i,j,k)             &
!                           +  L_e_i(i,j,k)              &
                           +  Lrot_e_N2(i,j,k)          &
                           + Lrot_e_O2(i,j,k)!           &
!                           + Lf_e_O(i,j,k) 


           cooling(i,j,k) = cooling(i,j,k) * 1.602e-19 * 1e6 * 0.5
!!!!!! change the units eV cm-3 s-1  to Joule m-3 s-1 !!!!!!!!!

           if (altitude(k)/1000. <250.) then
              cooling(i,j,k) = 0.99*Heating(i,j,k)
           endif

       enddo
           

        do k = 1, nAlts           

           if (cooling(i,j,k)> 1.5*Heating(i,j,k)) &
           cooling(i,j,k)= 1.5*Heating(i,j,k)
           if (cooling(i,j,k)< 0.6*Heating(i,j,k)) &
           cooling(i,j,k)= 0.6*Heating(i,j,k)


        enddo

  enddo
enddo

end subroutine calc_etemp_sources







