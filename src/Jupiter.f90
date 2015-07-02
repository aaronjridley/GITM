!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine fill_photo(photoion, photoabs)

  use ModPlanet
  use ModEUV

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpecies)

  photoabs           = 0.0
  photoabs(:,iH_)    = PhotoAbs_H
  photoabs(:,iH2_)   = PhotoAbs_H2

  photoion           = 0.0
  photoion(:,iHP_)   = PhotoIon_H
  photoion(:,iH2P_)  = PhotoIon_H2

end subroutine fill_photo


!---------------------------------------------------------------------
! Initialize Heating Rates
!---------------------------------------------------------------------

subroutine init_heating_efficiency

  use ModEUV, only: HeatingEfficiency, eHeatingEfficiency

  implicit none

  !! CHANGE
  HeatingEfficiency  = 0.3
  eHeatingEfficiency = 0.0

  call init_radcool

end subroutine init_heating_efficiency

!---------------------------------------------------------+
subroutine calc_planet_sources(iBlock)

!  All new source code in here

  use ModInputs
  use ModSources
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

! New sources specificly for Mars include: 
! (1) calc_radcooling(iBlock):  added 1/31/07 (BOUGHER)
! (2) calc_radcode(iBlock)   :  to be added later (NELLI)

  call calc_radcooling(iBlock)

end subroutine calc_planet_sources

!---------------------------------------------------------+
subroutine calc_radcooling(iBlock)

!  ---------------------------------------------------------------------
!  Purpose: to calculate NLTE CO2 15-micron cooling in Mars upper atmosphere
!           (above 80 km, where NLTE physics becomes important)
!  Source: nltecool.F from mtgcm16 
!  Coders: M. Lopez-Valverde (2001)
!          S. W. Bougher (2004-MTGCM15, 2007-Mars GITM)
!  Inputs: from MarsGITM  (Pressure, Temperature, VMRCO2, VMRO, VMRN2CO)
!          all on MarsGITM grid (nLons,nLats,nAlts)
!  Inputs: from init_radcool (pnbr, ef1, ef2, co2vmr,o2pvmr,n2covmr)
!          all on np=68 1-D NLTE code grid
!  Output: COOLTOT (K/sec) and RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock)
!  =======================================================================
!
!     subroutine nltecool(nlayer,player,tlayer,dt)
!
! This code was designed as a delivery for the "Martian Environment Models" 
! project ( ESA contract 11369/95/nl/jg CCN2 )
! Computes non-LTE heating rates from CO2 emission at 15 um
! in the Martian upper atmosphere.
! Uses a simplified model consisting of two excited levels with two
! emission bands, one of them stronger than the other, which correspond
! to the behaviours of the 626 fundamental band and the isotopic fund.bands.
! It uses a cool-to-space approximation with tabulated escape functions.
! These escape functions have been precomputed for the strong and weak bands,
! and are given as a function of pressure in separate files.
! The output values are the heating rates (actually, cooling, since they
! are always negative) for the two bands, i.e., the total cooling is the
! sum of them.
! Miguel A. Lopez Valverde
! Instituto de Astrofisica de Andalucia (CSIC), Granada, Spain
!
! Version 1b.  See description above.  22-March-2000.
! Adapted as a subroutine for use in GCM -- PLR/SRL 6/2000
! Version 1c.  Inclusion of VMR in the tabulation of escape functions. 
!              Table contains now only 1 input file -- Miguel 11/Jul/2000
! Nov-2001   Adapt to the MZ1D model (1D version & extended parametp_table)
! Nov-2005   Adapt to the MTGCM16 model (3D version & extended parametp_table)
! Jan-2007   Adapt to the MarsGITM model (3D version & extended parametp_table)
!  -----------------------------------------------------------------

  use ModInputs, only:  iDebugLevel
  use ModSources, only: RadCoolingRate
  use ModPlanet
  use ModGITM
  use ModConstants, only:  Boltzmanns_Constant, Speed_Light, &
       Planck_Constant
  use ModIndicesInterfaces

  implicit none

! INPUTS:----------------------------------------------------------

  integer,intent(in):: iBlock 

!  -----------------------------------------------------------------
! Input variables
  
  integer, parameter :: nlayer = nAlts ! no. of atmospheric layers
  integer, parameter :: nlon = nLons ! no. of long gridpoints
  integer, parameter :: nlat = nLats ! no. of latitude gridpoints
  real :: player(nlayer)            ! input pressure grid [Pa]
  real :: tlayer(nlayer)            ! input temperatures [K]

!  -----------------------------------------------------------------

! Standard atmosphere variables (overlay with MarsGITM values)
 
  real :: nt , co2(nlayer) , o3p(nlayer) , n2co(nlayer)
 
  ! Vectors/indexes for the datavol2 tabulation of escape functions and VMR
 
  !     Number data points in Tabulation (in ModMars module)
  !     integer,parameter :: np=68 
  !     Dimensioned above from ModMars module
  !     real ::  pnbr(np), ef1(np), ef2(np), co2vmr(np), o3pvmr(np),n2covmr(np)
  !     Reference Pressure (Pascals) from 1-D NLTE code

  real ::  pnb(np) 

  !     Interpolated escape functions

  real ::  escf1(nlayer) , escf2(nlayer) 

! Local Constants and variables
 
  real  ::  n1, n2, co2t , l1, p1, p12 , l2, p2, p21
  real  ::  tt, c1, c2, ae1, ae2,  a1, a2, a12, a21
  real  ::  pl1, pl2, el1, el2, hr1, hr2, x 
  real  ::  hr(nlayer)  

  ! Indexes
 
  integer  :: i,j,ii
  integer :: iLat, iLon, iAlt
 
  ! Rate coefficients

  real  ::  k19xca, k19xcb, k19cap1, k19cap2, k19cbp1, k19cbp2
  real  ::  d19c, d19cp1, d19cp2, k20xc, k20cp1, k20cp2
  real  ::  k21xc, k21cp2

  !
  ! Data and settings
  !
  real,parameter :: nu1=667.38, nu2=662.3734, hplanck=6.6261e-27,  &
       gamma1=1.191e-5, vlight=3.e10, ee=1.438769
  real,parameter :: imr1=0.987, imr2=0.00408 + 0.0112,  &
       rfvt=0.1, rfvto3p=1.0, rfvv=0.1

  !  -----------------------------------------------------------------
  ! Internal fields recast from MarsGITM

  real,dimension(1:nLons,1:nLats,1:nAlts) ::    &
       T, P, vmrco2, vmro, vmrn2co, cooltot, mnd 
  
  !     Mars GITM real temperature (on its grid)  ----------------------

  T(1:nLons,1:nLats,1:nAlts) = &
       Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*  &
       TempUnit(1:nLons,1:nLats,1:nAlts)

  !     Mars GITM pressure in pascals  (on its grid) -------------------

  P(1:nLons,1:nLats,1:nAlts) = &
       Pressure(1:nLons,1:nLats,1:nAlts,iBlock)

  !    Volume Mixing Ratio Calculations

  mnd(1:nLons,1:nLats,1:nAlts) = &  
       NDensity(1:nLons,1:nLats,1:nAlts,iBlock)+1.0

  vmro(1:nLons,1:nLats,1:nAlts)  = & 
       NdensityS(1:nLons,1:nLats,1:nAlts,iO_,iBlock)/ &
       mnd(1:nLons,1:nLats,1:nAlts)
  vmrn2co(1:nLons,1:nLats,1:nAlts) = & 
       (NdensityS(1:nLons,1:nLats,1:nAlts,iCO_,iBlock)  &
       + NdensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock))/ &
       mnd(1:nLons,1:nLats,1:nAlts) 
  vmrco2(1:nLons,1:nLats,1:nAlts) = &  
       NdensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock) &
       /mnd(1:nLons,1:nLats,1:nAlts)

  !  -----------------------------------------------------------------
  !
  ! **** ARRAYS made available via init_radcool 
  ! **** Convert pnb array to pascal units here
  !
  do i=1,np
     pnb(i)=1.0e-4*exp(pnbr(i)) ! p into Pa
  enddo

  ! ****  Initialize to zero
  cooltot(1:nLons,1:nLats,1:nAlts) = 0.0 

  !-------------------------------------------------------------

  ! MAIN LOOP, for MarsGITM horizontal grid (Latitude,Longitude) 
 
  do iLat = 1, nlat
     do iLon = 1, nlon

        ! *****  populate 1-D arrays(nlayer) for calc_radcooling from MarsGITM 
        do  iAlt=1,nlayer
           player(iAlt)= P(iLon,iLat,iAlt)
           tlayer(iAlt)= T(iLon,iLat,iAlt)
           co2(iAlt)=vmrco2(iLon,iLat,iAlt)
           o3p(iAlt)=vmro(iLon,iLat,iAlt)
           n2co(iAlt)=vmrn2co(iLon,iLat,iAlt)
        enddo
 
        !        interpolate escape functions (only) to nlayer grid
        call interp1(escf2,player,nlayer,ef2,pnb,np)
        call interp1(escf1,player,nlayer,ef1,pnb,np)
 
        ! ALTITUDE LOOP, for each altitude (for 1-latitude/longitude profile):
        ! (from old nltecool.F routine)
 
        do i=1,nlayer  

           if (player(i).gt.3.5 .or. player(i).lt.4.0e-6) then 

              hr(i)=0.0
              cooltot(iLon,iLat,i)= 0.0

           else

              nt = player(i)/(1.381e-17*tlayer(i)) 
              co2(i)=co2(i)*nt                
              o3p(i)=o3p(i)*nt                
              n2co(i)=n2co(i)*nt
              n1 = co2(i) * imr1
              n2 = co2(i) * imr2
              co2t = n1 + n2

              tt = tlayer(i)*tlayer(i)
              k19xca = 4.2e-12 * exp( -2988.0/tlayer(i) + 303930.0/tt )
              k19xcb = 2.1e-12 * exp( -2659.0/tlayer(i) + 223052.0/tt )
              if (tlayer(i).le.175.0) then
                 k19xca = 3.3e-15
                 k19xcb = 7.6e-16
              endif
              k19xca = k19xca * rfvt
              k19xcb = k19xcb * rfvt
              k19cap1 = k19xca * 2.0 * exp( -ee*nu1/tlayer(i) )
              k19cap2 = k19xca * 2.0 * exp( -ee*nu2/tlayer(i) )
              k19cbp1 = k19xcb * 2.0 * exp( -ee*nu1/tlayer(i) )
              k19cbp2 = k19xcb * 2.0 * exp( -ee*nu2/tlayer(i) )
              d19c = k19xca*co2t + k19xcb*n2co(i)
              d19cp1 = k19cap1*co2t + k19cbp1*n2co(i)
              d19cp2 = k19cap2*co2t + k19cbp2*n2co(i)
                                !
              k20xc = 3.e-12 * rfvto3p
              k20cp1 = k20xc * 2.0 * exp( -ee/tlayer(i) * nu1 )
              k20cp2 = k20xc * 2.0 * exp( -ee/tlayer(i) * nu2 )
                                !
              k21xc = 2.49e-11 * 0.5 * rfvv
              k21cp2 = k21xc * exp( - ee/tlayer(i) * (nu2-nu1) )
                                !
              l1 = d19c + k20xc*o3p(i) + k21cp2*n2
              p1 = ( d19cp1 + k20cp1*o3p(i) ) * n1
              p12 = k21xc*n1
                                !
              l2 = d19c + k20xc*o3p(i) + k21xc*n1
              p2 = ( d19cp2 + k20cp2*o3p(i) ) * n2
              p21 = k21cp2*n2

              ae1 = 1.3546 * 1.66 / 4.0 * escf1(i)
              ae2 = ( 1.3452 + 1.1878 ) * 1.66 / 4.0 * escf2(i)
              l1 = l1 + ae1
              l2 = l2 + ae2

              c1 = gamma1*nu1**3. * 0.5
              c2 = gamma1*nu2**3. * 0.5
              a1 = c1 * p1 / (n1*l1)
              a2 = c2 * p2 / (n2*l2)
              a12 = (nu1/nu2)**3. * n2/n1 * p12/l1
              a21 = (nu2/nu1)**3. * n1/n2 * p21/l2
              el2 = (a2 + a21 * a1 ) / ( 1.0 - a21 * a12 )
              el1 = a1 + a12 * el2
              pl1 = el1 * n1 / c1
              pl2 = el2 * n2 / c2

              hr1 = - hplanck*vlight * nu1 * ae1 * pl1
              hr2 = - hplanck*vlight * nu2 * ae2 * pl2
              hr(i) = hr1 + hr2
              !        Cooling Rate (K/sec units)
              cooltot(iLon,iLat,i)= 0.1*hr(i)*tlayer(i)/(4.4*player(i))

           endif
        enddo  !-------- END OF MAIN ALTITUDE LOOP  

     enddo  !-------- END OF MAIN LONGITUDE LOOP  
  enddo  !-------- END OF MAIN LATITUDE LOOP  

  !-------------------------------------------------------------

  RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock) = &
       -cooltot(1:nLons,1:nLats,1:nAlts) 

  !-------------------------------------------------------------

end subroutine calc_radcooling

!---------------------------------------------------------+

subroutine interp1(escout,p,nlayer,escin,pin,nl)
  implicit none
 
! subroutine to perform linear interpolation in pressure from 1D profile 
! escin(nl) sampled on pressure grid pin(nl) to profile
! escout(nlayer) on pressure grid p(nlayer).
!
! Input args:
!     integer,intent(in):: iBlock 

  integer,intent(in) :: nlayer,nl
  real,dimension(nl),intent(in) :: escin,pin
  real,dimension(nlayer),intent(in) :: p
  real,dimension(nlayer),intent(out) :: escout

!
! Local :
  integer :: n,nm,np,n1
  real :: wm,wp
!
  do n1=1,nlayer
     if(p(n1) .gt. 800.0 .or. p(n1) .lt. 5.0e-9) then
        escout(n1) = 0.0
     else
        do n = 1,nl-1
           if (p(n1).le.pin(n).and.p(n1).ge.pin(n+1)) then
              nm=n
              np=n+1
              wm=abs(pin(np)-p(n1))/(pin(nm)-pin(np))
              wp=1.0 - wm
           endif
        enddo
        !write(*,*) 'nm =',nm
        escout(n1) = escin(nm)*wm + escin(np)*wp
     endif
  enddo

end subroutine interp1

!---------------------------------------------------------+

subroutine init_radcool
!
! ---------------------------------------------------------
! Call to init_radcool will return corresponding set of six 1-D arrays:
!    pnb(np),ef1(np),ef2(np),co2vmr(np),o3pvmr(np), n2covmr(np)
!    These are 1-D NLTE cooling code inputs from Miguel Lopez-Valverde
! ---------------------------------------------------------

  use ModPlanet
  use ModGITM

!\
! Below we prescribe the 6-arrays needed
!/

  pnbr = (/ 12.0000, 11.0000, 10.8000, 10.6000, 10.4000, 10.2000, &
       10.0000, 9.80000, 9.60000, 9.40000, 9.20000, 9.00000, 8.80000,   & 
       8.60000, 8.40000, 8.20000, 8.00000, 7.80000, 7.60000, 7.40000,   &
       7.20000, 7.00000, 6.80000, 6.60000, 6.40000, 6.20000, 6.00000,   &
       5.80000, 5.60000, 5.40000, 5.20000, 5.00000, 4.80000, 4.60000,   &
       4.40000, 4.20000, 4.00000, 3.80000, 3.60000, 3.40000, 3.20000,   &
       3.00000, 2.80000, 2.60000, 2.40000, 2.20000, 2.00000, 1.80000,   &
       1.60000, 1.40000, 1.20000, 1.00000, 0.80000, 0.59999, 0.40000,   &
       0.20000, 0.00000,-0.20000,-0.40000,-0.60000,-0.80000,-1.00000,   &
      -1.20000,-1.40000,-1.60000,-1.80000,-2.00000,-3.00000/)

  ef1 = (/  0.000458112, 0.000458112,                              &
       0.000458112, 0.000458112, 0.000458112, 0.000458112, 0.000458112,  &
       0.000458112, 0.000458112, 0.000458112, 0.000458112, 0.000458112,  &
       0.000458112, 0.000458112, 0.000458112, 0.000458112, 0.000458112,  &
       0.000458112, 0.000458112, 0.000458112, 0.000458112, 0.000458112,  &
       0.000461707, 0.000476886, 0.000495638, 0.000520935, 0.000555511,  &
       0.000601219, 0.000663734, 0.000750691, 0.000863474, 0.001009000,  &
       0.001196420, 0.001426900, 0.001713980, 0.002066630, 0.002489740,  &
       0.003015780, 0.003643500, 0.004403230, 0.005320660, 0.006404560,  &
       0.007720690, 0.009256840, 0.011090500, 0.013237400, 0.015764300,  &
       0.018738800, 0.022207200, 0.026309900, 0.031061400, 0.036694800,  &
       0.043237300, 0.051502200, 0.062145500, 0.077721200, 0.099202700,  &
       0.131155000, 0.179470000, 0.258913000, 0.380549000, 0.530450000,  &
       0.643180000, 0.741061000, 0.826336000, 0.922787000, 0.997203000,  &
       1.000000000/)

  ef2 = (/  0.00198992,                                            &
       0.001989920, 0.001989920, 0.001989920, 0.001989920, 0.001989920,  &
       0.001989920, 0.001989920, 0.001989920, 0.001989920, 0.001989920,  &
       0.002013760, 0.002094500, 0.002229930, 0.002420560, 0.002680180,  &
       0.003043980, 0.003438960, 0.003802820, 0.004206220, 0.004761210,  &
       0.008016980, 0.011994700, 0.016914900, 0.022449700, 0.028524400,  &
       0.035481300, 0.043926400, 0.054624800, 0.067536700, 0.082993100,  &
       0.101717000, 0.123422000, 0.148468000, 0.177096000, 0.208816000,  &
       0.244003000, 0.282013000, 0.322559000, 0.365542000, 0.410518000,  &
       0.457384000, 0.505358000, 0.553627000, 0.600472000, 0.644807000,  &
       0.687185000, 0.727429000, 0.764734000, 0.798562000, 0.828699000,  &
       0.854797000, 0.877717000, 0.897874000, 0.915258000, 0.929904000,  &
       0.942381000, 0.952906000, 0.962173000, 0.970191000, 0.976437000,  &
       0.981501000, 0.985406000, 0.988560000, 0.991111000, 0.993653000,  &
       0.995561000, 1.000000000/)

  co2vmr =  (/                                                 &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.949619, 0.947694, 0.945830, 0.944016, 0.940557, 0.937068,  &
       0.932366, 0.893661/)

  o3pvmr =  (/ 0.000000000, 0.000000000, 0.000000000, 0.000000000,   &
       0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000,   &
       0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000,   &
       0.000100148, 0.000118846, 0.000139681, 0.000164909, 0.000193617,   &
       0.000225161, 0.000260834, 0.000301501, 0.000344953, 0.000391011,   &
       0.000440377, 0.000490820, 0.000543200, 0.000595335, 0.000645420,   &
       0.000693166, 0.000743729, 0.000793710, 0.000844394, 0.000894318,   &
       0.000944732, 0.000994964, 0.001049010, 0.001100080, 0.001163020,   &
       0.001229890, 0.001300260, 0.001371310, 0.001455560, 0.001551860,   &
       0.001663280, 0.001778020, 0.001915460, 0.002075030, 0.002249030,   &
       0.002471170, 0.002717280, 0.002997390, 0.003335820, 0.003735070,   &
       0.004208190, 0.004768870, 0.005425580, 0.006208150, 0.007144730,   &
       0.008285450, 0.009517790, 0.010814000, 0.012235900, 0.013687000,   &
       0.015149500, 0.016719600, 0.018548500, 0.033625200/)

  n2covmr =  (/  0.027145, 0.027145,                                &
       0.0271450, 0.0271450, 0.0271450, 0.0271450, 0.0271450, 0.0271450, &
       0.0271450, 0.0271450, 0.0271450, 0.0271450, 0.0271450, 0.0271450, &
       0.0272661, 0.0272848, 0.0273054, 0.0273279, 0.0273514, 0.0273775, &
       0.0274048, 0.0274345, 0.0274672, 0.0275021, 0.0275404, 0.0275826, &
       0.0276340, 0.0277013, 0.0278220, 0.0279707, 0.0281759, 0.0284339, &
       0.0287587, 0.0291600, 0.0296561, 0.0302558, 0.0309922, 0.0318062, &
       0.0327010, 0.0335635, 0.0344388, 0.0353327, 0.0362143, 0.0370941, &
       0.0379315, 0.0387626, 0.0395524, 0.0403071, 0.0410071, 0.0416229, &
       0.0421231, 0.0425167, 0.0427964, 0.0429773, 0.0430488, 0.0429638, &
       0.0428049, 0.0426788, 0.0426822, 0.0429426, 0.0434634, 0.0442559, &
       0.0453038, 0.0465879, 0.0480262, 0.0496303, 0.0514885,            &
       0.0691651/)
!
end subroutine init_radcool
! ----------------------------------------------------------------------
