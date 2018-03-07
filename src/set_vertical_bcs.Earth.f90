!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! ------------------------------------------------------------
! set_boundary
! ------------------------------------------------------------
!/

subroutine set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp, LogINS, iVel, VertVel)

  ! Fill in ghost cells at the top and bottom

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect, Mass, nIons, IsEarth,&
                       iN2_,iNO_,iO_3P_, iN_4S_
  use ModGITM, only: iEast_, iNorth_, iUp_
  use ModInputs
  use ModConstants
  use ModTime, only: UTime, iJulianDay,currenttime
  use ModVertical, only: &
       Lat, Lon, Gravity_G, &
       Altitude_G, dAlt_F, &
       MeanMajorMass_1d, &
       iLon1D, iLat1D, iBlock1D, &
       Centrifugal, InvRadialDistance_C, Coriolis, &
       MLatVertical, SZAVertical
  use ModIndicesInterfaces, only: get_HPI
  use ModTides, only: TidesNorth, TidesEast, TidesTemp

  use EUA_ModMsis90, ONLY: meter6

  implicit none

  real, intent(inout) :: &
       LogRho(-1:nAlts+2), &
       LogNS(-1:nAlts+2,nSpecies), &
       LogINS(-1:nAlts+2,nIons), &
       Vel_GD(-1:nAlts+2,3), &
       IVel(-1:nAlts+2,3), &
       Temp(-1:nAlts+2), &
       VertVel(-1:nAlts+2,nSpecies)

  integer :: iSpecies, iAlt, iDir
  real    :: InvScaleHeightS, InvScaleHgt, Alt, Lst, Ap = 4.0, dn, dt
  logical :: IsFirstTime = .true., UseMsisBCs = .false.
  real    :: HP, v(2)
  integer :: ierror
  real    :: temptemp
  real    :: logNS_Species(nSpecies)
  real :: n0, n1, n2, n3, n4, n5
  real :: tec
  
  integer, dimension(25) :: sw

!! JMB Added 11-03-2014 (New BC Variables)
!! Useful for Upper Boundary Conditions--Hydrostatic
  real :: InvAtmScaleHeight
  real :: NS(-1:nAlts+2,1:nSpecies), NT(-1:nAlts+2)
  real :: SumRho
  real :: EffectiveGravity(-1:nAlts+2)
  real :: MeanGravity, MeanMass, MeanTemp
  logical :: IsHydrostatic(1:nSpecies), IsPhotoChemical(1:nSpecies)

  logical :: UsePlasmasphereBC
  
 ! Gradient Terms
  real :: dLogNS, dTemp, dVertVel
  real :: dLogINS


!!! Use for 4-th Order Forward Differences
!!! Need a 5-point Stencil
  real :: h1, h2, h3, h4
  real :: MeshH1, MeshH2, MeshH3, MeshH4
  real :: MeshCoef0, MeshCoef1, &
          MeshCoef2, MeshCoef3, &
          MeshCoef4

!!! Use for 4-th Order Backward Differences
!!! Need a 5-point Stencil
  real :: hm1, hm2, hm3, hm4
  real :: MeshHm1, MeshHm2, MeshHm3, MeshHm4
  real :: MeshCoefm0, MeshCoefm1, &
          MeshCoefm2, MeshCoefm3, &
          MeshCoefm4

  !-----------------------------------------------------------
  ! Bottom
  !-----------------------------------------------------------

  IsPhotoChemical(1:nSpecies) = .false.
  IsPhotoChemical(iN_4S_) = .true.
  IsPhotoChemical(iNO_) = .true.
!  IsPhotoChemical(iO_3P_) = .true.

  NS(-1:nAlts+2,1:nSpecies) = exp(LogNS(-1:nAlts+2,1:nSpecies))

!  do iAlt = -1, nAlts+2
!     EffectiveGravity(iAlt) = Gravity_G(iAlt) 
!  enddo 
!
  do iAlt = -1, nAlts + 2
     EffectiveGravity(iAlt) = &
        Gravity_G(iAlt) + &
        Centrifugal / InvRadialDistance_C(iAlt) 
  enddo 

  if (IsEarth) UseMsisBCs = UseMsis

  if (IsFirstTime .and. UseMsisBCs) then
     call meter6(.true.)
     sw = 1
     IsFirstTime = .true.
  endif

  if (UseMsisBCs) then
     call get_HPI(CurrentTime, HP, iError)  
     if (iError > 0) hp = 40.0
     Ap = min(200.,max(-40.72 + 1.3 * HP, 10.))
     do iAlt = -1, 0
        Alt = Altitude_G(iAlt)/1000.0
        Lst = mod(UTime/3600.0+Lon/15.0,24.0)

        call msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS_Species, temptemp, &
             LogRho(iAlt),v)

        LogNS(iAlt,:) = logNS_Species
        
        if (.not. DuringPerturb) temp(iAlt) = temptemp

        vel_gd(iAlt,iEast_) = v(iEast_)
        vel_gd(iAlt,iNorth_) = v(iNorth_)

     enddo
  else
     ! Don't Let the winds blow
     Vel_GD(-1:0,iEast_)  = 0.0
     Vel_GD(-1:0,iNorth_) = 0.0
  endif

  if (.not. DuringPerturb) then
     Vel_GD(-1:0,iUp_)    = 0.0
     VertVel(-1:0,:)      = 0.0
  endif

  if (UseGSWMTides) then
     Vel_GD(-1:0,iEast_)  = TidesEast(iLon1D,iLat1D,1:2,iBlock1D)
     Vel_GD(-1:0,iNorth_) = TidesNorth(iLon1D,iLat1D,1:2,iBlock1D)
     Temp(-1:0)           = TidesTemp(iLon1D,iLat1D,1:2,iBlock1D) + Temp(-1:0)
  endif
  if (UseWACCMTides) then
     Vel_GD(-1:0,iEast_)  = TidesEast(iLon1D,iLat1D,1:2,iBlock1D)
     Vel_GD(-1:0,iNorth_) = TidesNorth(iLon1D,iLat1D,1:2,iBlock1D)
     Temp(-1:0)           = TidesTemp(iLon1D,iLat1D,1:2,iBlock1D)
  endif
  ! Update the -1 Cell only for Temp, LogNS, and Vel_GD
  ! The 0 Cell is set by MSIS or User-supplied Settings
  ! Need to Calculate 0ne-sided first derivative
  ! at the Cell 0.

  ! Calculate the non-uniform mesh coefficients
  iAlt = -1
  h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
  h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
  h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
  h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

  ! Mesh Coefficients are summations over the individual mesh scales
  MeshH1 = h1                 
  MeshH2 = h1 + h2            
  MeshH3 = h1 + h2 + h3
  MeshH4 = h1 + h2 + h3 + h4

  !!! 3rd Order Mesh Coef
  MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                     MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                   (MeshH1*MeshH2*MeshH3*MeshH4) 
  MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/(h1*h2*(h2 + h3)*(h2 + h3 + h4))
  MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
  MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
  MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/(MeshH4*(h2+h3+h4)*(h3+h4)*h4)

  !  dTemp = MeshCoef0*Temp(iAlt+1) + &  
  !          MeshCoef1*Temp(iAlt+2) + &  
  !          MeshCoef2*Temp(iAlt+3) + &  
  !          MeshCoef3*Temp(iAlt+4) + &  
  !          MeshCoef4*Temp(iAlt+5)      
!
!    Temp(iAlt) = Temp(iAlt+1)&
!                 - dAlt_F(iAlt+1)*dTemp 

! For the MSIS BCS
iAlt = -1
    Temp(iAlt) = Temp(iAlt+1)
    do iDir = 1, 3
    Vel_GD(iAlt  ,iDir) = Vel_GD(iAlt+1,iDir) 
    enddo 


!    do iSpecies = 1, nSpecies
!       dLogNS = MeshCoef0*LogNS(iAlt+1,iSpecies) + &  
!                MeshCoef1*LogNS(iAlt+2,iSpecies) + &  
!                MeshCoef2*LogNS(iAlt+3,iSpecies) + &  
!                MeshCoef3*LogNS(iAlt+4,iSpecies) + &  
!                MeshCoef4*LogNS(iAlt+5,iSpecies)      
!
!       LogNS(iAlt,iSpecies) = LogNS(iAlt+1,iSpecies)&
!                 - dAlt_F(iAlt+1)*dLogNS 
!    enddo 

    ! Set the Neutral Bulk Winds

!    do iDir = 1, 3
!!       do iAlt = 0, -1, -1
!
!          h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
!          h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
!          h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
!          h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes
!
!          ! Mesh Coefficients are summations over the individual mesh scales
!          MeshH1 = h1                 
!          MeshH2 = h1 + h2            
!          MeshH3 = h1 + h2 + h3
!          MeshH4 = h1 + h2 + h3 + h4
!
!          !!! 3rd Order Mesh Coef
!          MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
!                             MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
!                           (MeshH1*MeshH2*MeshH3*MeshH4) 
!          MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/&
!                            (h1*h2*(h2 + h3)*(h2 + h3 + h4))
!          MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
!          MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
!          MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/&
!                            (MeshH4*(h2+h3+h4)*(h3+h4)*h4)
!
!          dVertVel = Vel_GD(iAlt+1,iDir)*MeshCoef0 + & 
!                     Vel_GD(iAlt+2,iDir)*MeshCoef1 + &
!                     Vel_GD(iAlt+3,iDir)*MeshCoef2 + &
!                     Vel_GD(iAlt+4,iDir)*MeshCoef3 + &
!                     Vel_GD(iAlt+5,iDir)*MeshCoef4
!
!          Vel_GD(iAlt  ,iDir) = Vel_GD(iAlt+1,iDir) - &
!                              dAlt_F(iAlt+1)*dVertVel 
!       enddo !iAlt = 0, -1, -1
!    enddo  ! iDir
!  endif !(.not. UseMSISBCs) then

  
  ! Do the following if we DO have MSIS BCS
  ! For MSIS, N0, T0, and V(bulk) are set in cell 0
  do iSpecies = 1, nSpecies

     if (.not. IsPhotoChemical(iSpecies)) then
       ! This is what we do when (1) We're using MSIS and 
       ! (2) the species is NOT Photochemical
       iAlt = -1
       MeanGravity = -0.5*(EffectiveGravity(iAlt  ) + &
                           EffectiveGravity(iAlt+1))
       MeanTemp =  0.5*( Temp(iAlt+1) + Temp(iAlt) )
       MeanMass = 0.5*(MeanMajorMass_1d(iAlt+1) + MeanMajorMass_1d(iAlt))
       InvScaleHeightS =  MeanGravity * MeanMass / &
                          (MeanTemp*Boltzmanns_Constant)

       NS(iAlt,iSpecies) = NS( iAlt+1,iSpecies)*&
                          (Temp(iAlt+1)/Temp(iAlt))*&
                     exp( +1.0*InvScaleHeightS*dAlt_F(iAlt)) 

      LogNS(iAlt  ,iSpecies) = alog(NS(iAlt,iSpecies))

      if (iSpecies .eq. iO_3P_) then
         LogNS(iAlt,iSpecies) = LogNS(0,iSpecies)  ! Assume 0 gradient below boundary
      endif 
!       iAlt = -1
!       h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
!       h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
!       h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
!       h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes
!
!       ! Mesh Coefficients are summations over the individual mesh scales
!       MeshH1 = h1                 
!       MeshH2 = h1 + h2            
!       MeshH3 = h1 + h2 + h3
!       MeshH4 = h1 + h2 + h3 + h4
!
!       !!! 3rd Order Mesh Coef
!       MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
!                          MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
!                        (MeshH1*MeshH2*MeshH3*MeshH4) 
!       MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/&
!                         (h1*h2*(h2 + h3)*(h2 + h3 + h4))
!       MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
!       MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
!       MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/&
!                         (MeshH4*(h2+h3+h4)*(h3+h4)*h4)
!!
!       dLogNS = LogNS(iAlt+1,iSpecies)*MeshCoef0 + & 
!                LogNS(iAlt+2,iSpecies)*MeshCoef1 + &
!                LogNS(iAlt+3,iSpecies)*MeshCoef2 + &
!                LogNS(iAlt+4,iSpecies)*MeshCoef3 + &
!                LogNS(iAlt+5,iSpecies)*MeshCoef4
!!
!       LogNS(iAlt  ,iSpecies) = LogNS(iAlt+1,iSpecies) - &
!                                  dAlt_F(iAlt+1)*dLogNS 
     else 
        do iAlt = 0,-1,-1

          h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
          h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
          h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
          h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

          ! Mesh Coefficients are summations over the individual mesh scales
          MeshH1 = h1                 
          MeshH2 = h1 + h2            
          MeshH3 = h1 + h2 + h3
          MeshH4 = h1 + h2 + h3 + h4

          !!! 3rd Order Mesh Coef
          MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                             MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                           (MeshH1*MeshH2*MeshH3*MeshH4) 
          MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/&
                            (h1*h2*(h2 + h3)*(h2 + h3 + h4))
          MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
          MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
          MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/&
                            (MeshH4*(h2+h3+h4)*(h3+h4)*h4)

          dLogNS = LogNS(iAlt+1,iSpecies)*MeshCoef0 + & 
                   LogNS(iAlt+2,iSpecies)*MeshCoef1 + &
                   LogNS(iAlt+3,iSpecies)*MeshCoef2 + &
                   LogNS(iAlt+4,iSpecies)*MeshCoef3 + &
                   LogNS(iAlt+5,iSpecies)*MeshCoef4

          LogNS(iAlt  ,iSpecies) = LogNS(iAlt+1,iSpecies) - &
                                     dAlt_F(iAlt+1)*dLogNS 

          !!! Limit the Gradients through the lower boundary
          if (dLogNS .le. 0.0) then
              LogNS(iAlt,iSpecies) = LogNS(iAlt+1,iSpecies)
          endif 
        enddo !iAlt = 0,-1,-1

     endif ! PhotoChemical Check
  enddo  ! iSpecies loop

  ! Ion Lower Boundaries
  do iAlt = 0,-1,-1

    h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
    h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
    h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
    h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

    ! Mesh Coefficients are summations over the individual mesh scales
    MeshH1 = h1                 
    MeshH2 = h1 + h2            
    MeshH3 = h1 + h2 + h3
    MeshH4 = h1 + h2 + h3 + h4

    !!! 3rd Order Mesh Coef
    MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                       MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                     (MeshH1*MeshH2*MeshH3*MeshH4) 
    MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/(h1*h2*(h2 + h3)*(h2 + h3 + h4))
    MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
    MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
    MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/(MeshH4*(h2+h3+h4)*(h3+h4)*h4)

    ! Ions Float at the LBC
    do iSpecies = 1, nIonsAdvect

       ! Simply assume no gradient, since this is a region in which chemitry
       ! is dominant.  It shouldn't matter, really:
       LogINS(iAlt,iSpecies) = LogINS(iAlt+1,iSpecies)

!        dLogINS = MeshCoef0*LogINS(iAlt+1,iSpecies) + &  ! LogNS(0)
!                  MeshCoef1*LogINS(iAlt+2,iSpecies) + &  ! LogNS(1)
!                  MeshCoef2*LogINS(iAlt+3,iSpecies) + &  ! LogNS(2)
!                  MeshCoef3*LogINS(iAlt+4,iSpecies) + &  ! LogNS(2)
!                  MeshCoef4*LogINS(iAlt+5,iSpecies)      ! LogNS(2)
!
!        ! Make sure that the ions decrease in the ghost cells!
!        if (dLogINS .gt. 0.0) then
!          LogINS(iAlt,iSpecies) = LogINS(iAlt+1,iSpecies)&
!                              - dAlt_F(iAlt+1)*dLogINS 
!        else
!          LogINS(iAlt,iSpecies) = LogINS(iAlt+1,iSpecies)
!        endif 
    enddo ! Ions Advect 

    do iSpecies = 1, nSpecies
         ! For Photochemical Species
         ! We enforce zero flux at the LBC

         dVertVel = VertVel(iAlt+1,iSpecies)*MeshCoef0 + & 
                    VertVel(iAlt+2,iSpecies)*MeshCoef1 + &
                    VertVel(iAlt+3,iSpecies)*MeshCoef2 + &
                    VertVel(iAlt+4,iSpecies)*MeshCoef3 + &
                    VertVel(iAlt+5,iSpecies)*MeshCoef4

         VertVel(iAlt  ,iSpecies) = VertVel(iAlt+1,iSpecies) - &
                                     dAlt_F(iAlt+1)*dVertVel 
    enddo ! nSpecies

    ! Set the Bulk Winds
    Vel_GD(iAlt,iUp_) = 0.0
    do iSpecies = 1, nSpecies
       Vel_GD(iAlt,iUp_) = Vel_GD(iAlt,iUp_) + &
        NS(iAlt,iSpecies)*Mass(iSpecies)*&
        VertVel(iAlt,iSpecies)/exp(LogRho(iAlt))
    enddo 

    ! Set the Ion Bulk Winds
    do iDir = 1, 3
       dVertVel = IVel(iAlt+1,iDir)*MeshCoef0 + & 
                  IVel(iAlt+2,iDir)*MeshCoef1 + &
                  IVel(iAlt+3,iDir)*MeshCoef2 + &
                  IVel(iAlt+4,iDir)*MeshCoef3 + &
                  IVel(iAlt+5,iDir)*MeshCoef4

       IVel(iAlt  ,iDir) = IVel(iAlt+1,iDir) - &
                           dAlt_F(iAlt+1)*dVertVel 
    enddo  ! iDir

  enddo ! End Outer IAlt Loop (0, -1, -1) 

  ! Special Case for PhotoChemical Species
  ! Allow O to flow upward
  do iSpecies = 1, nSpecies
     if (IsPhotoChemical(iSpecies) .and. (VertVel(1,iSpecies) .gt. 0.0) &
          .and. (iSpecies .ne. iO_3P_) ) then 
        ! Do NOT allow photochemical species to flow upward!
        ! The problem becomes unconstrained
        VertVel( 0,iSpecies) = -1.0*VertVel(1,iSpecies)
        VertVel(-1,iSpecies) = -1.0*VertVel(1,iSpecies)
     else 
        ! Do nothing (already taken care of)
     endif  

  enddo 

  ! Neutral Bulk Winds (East and West)
  if (UseMsisBcs) then

     iAlt = -1
     h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
     h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
     h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
     h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

    ! Mesh Coefficients are summations over the individual mesh scales
     MeshH1 = h1                 
     MeshH2 = h1 + h2            
     MeshH3 = h1 + h2 + h3
     MeshH4 = h1 + h2 + h3 + h4

    !!! 3rd Order Mesh Coef
     MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                        MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                      (MeshH1*MeshH2*MeshH3*MeshH4) 
     MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/(h1*h2*(h2 + h3)*(h2 + h3 + h4))
     MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
     MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
     MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/(MeshH4*(h2+h3+h4)*(h3+h4)*h4)

     do iDir = iEast_, iNorth_
        dVertVel = Vel_GD(iAlt+1,iDir)*MeshCoef0 + & 
                   Vel_GD(iAlt+2,iDir)*MeshCoef1 + &
                   Vel_GD(iAlt+3,iDir)*MeshCoef2 + &
                   Vel_GD(iAlt+4,iDir)*MeshCoef3 + &
                   Vel_GD(iAlt+5,iDir)*MeshCoef4

        Vel_GD(iAlt  ,iDir) = Vel_GD(iAlt+1,iDir) - &
                                 dAlt_F(iAlt+1)*dVertVel 
     enddo 
  endif 

  !-----------------------------------------------------------
  ! Top
  !-----------------------------------------------------------

  ! Slip flow at the top
  ! Assume zero gradients in the velocities & temps
  Vel_GD(nAlts+1:nAlts+2,iEast_)  = Vel_GD(nAlts,iEast_)
  Vel_GD(nAlts+1:nAlts+2,iNorth_) = Vel_GD(nAlts,iNorth_)

  IVel(nAlts+1:nAlts+2,iEast_)  = IVel(nAlts,iEast_)
  IVel(nAlts+1:nAlts+2,iNorth_) = IVel(nAlts,iNorth_)

  ! Things can go up or down in the ions
  !IVel(nAlts+1,iUp_)   =  IVel(nAlts  ,iUp_)
  !IVel(nAlts+2,iUp_)   =  IVel(nAlts-1,iUp_)

  do iSpecies = 1, nSpecies
     VertVel(nAlts+1,iSpecies) =  1.0*VertVel(nAlts,iSpecies)
     VertVel(nAlts+2,iSpecies) =  1.0*VertVel(nAlts,iSpecies)

!    ! Can't have photochemical species flowing downward
    if (IsPhotoChemical(iSpecies) .and. (VertVel(nAlts,iSpecies) .lt. 0.0) ) then
       VertVel(nAlts+1,iSpecies) = -1.0*VertVel(nAlts  ,iSpecies)
       VertVel(nAlts+2,iSpecies) = -1.0*VertVel(nAlts-1,iSpecies)
    endif

  enddo 
!
  if (IVel(nAlts,iUp_) .lt. 0.0) then
     IVel(nAlts+1,iUp_) = -1.0*IVel(nAlts,iUp_)
     IVel(nAlts+2,iUp_) = -1.0*IVel(nAlts,iUp_)
  else
     IVel(nAlts+1,iUp_) = IVel(nAlts,iUp_)
     IVel(nAlts+2,iUp_) = IVel(nAlts,iUp_)
  endif

  if(Vel_GD(nAlts,iUp_) .lt. 0.0) then
     Vel_GD(nAlts+1,iUp_) = -1.0*Vel_GD(nAlts,iUp_)
     Vel_GD(nAlts+2,iUp_) = -1.0*Vel_GD(nAlts,iUp_)
  else
     Vel_GD(nAlts+1,iUp_) =  1.0*Vel_GD(nAlts,iUp_)
     Vel_GD(nAlts+2,iUp_) =  1.0*Vel_GD(nAlts,iUp_)
  endif

  ! Constant temperature (zero gradient)

  Temp(nAlts+1) = Temp(nAlts)
  Temp(nAlts+2) = Temp(nAlts)

  dn = (LogRho(nAlts) - LogRho(nAlts-1))
  LogRho(nAlts+1) = LogRho(nAlts) + dn
  LogRho(nAlts+2) = LogRho(nAlts+1) + dn

  ! Hydrostatic pressure for the neutrals

  do iSpecies=1,nSpecies
     do iAlt = nAlts+1, nAlts+2
        MeanGravity = -0.5*(EffectiveGravity(iAlt  ) + &
                            EffectiveGravity(iAlt-1))
        MeanTemp =  0.5*( Temp(iAlt-1) + Temp(iAlt) )

        InvScaleHeightS =  MeanGravity * Mass(iSpecies) / &
                           (MeanTemp*Boltzmanns_Constant)

        NS(iAlt,iSpecies) = NS( iAlt-1,iSpecies)*&
                           (Temp(iAlt-1)/Temp(iAlt))*&
              exp( -1.0*InvScaleHeightS*dAlt_F(iAlt)) 

        LogNS(iAlt,iSpecies) = alog(NS(iAlt,iSpecies))

        if (LogNS(nAlts+1,iSpecies) > 75.0 .or. &
             LogNS(nAlts+2,iSpecies) > 75.0) then
           write(*,*) "======> bcs : ", iSpecies, 1.0e-3/InvScaleHeightS, &
                Gravity_G(nAlts), Mass(iSpecies), Temp(nAlts), &
                LogNS(nAlts,iSpecies), LogNS(nAlts+1,iSpecies), &
                dAlt_F(nAlts), LogNS(nAlts+2,iSpecies)
        endif
     enddo
  enddo

  UsePlasmasphereBC = .false.
  if (UseNighttimeIonBCs) then
     if ( SZAVertical > Pi/2 .and. &
          (( MLatVertical > 30.0 .and. MLatVertical <70.0) .or. &
          ( MLatVertical > -60.0 .and. MLatVertical < -30.0 ))) then
        UsePlasmasphereBC = .true.
     endif
  endif

  if (UseImprovedIonAdvection) then
     tec = sum(dAlt_f(1:nAlts) * LogINS(1:nAlts,1))/1e16
  else
     tec = sum(dAlt_f(1:nAlts) * exp(LogINS(1:nAlts,1)))/1e16
  endif

  do iAlt = nAlts + 1, nAlts + 2

     hm1 = dAlt_F(iAlt-1) ! 
     hm2 = dAlt_F(iAlt-2) ! 
     hm3 = dAlt_F(iAlt-3) ! 
     hm4 = dAlt_F(iAlt-4) ! 

     ! Mesh Coefficients are summations over the individual mesh scales
     MeshHm1 = hm1                 
     MeshHm2 = hm1 + hm2            
     MeshHm3 = hm1 + hm2 + hm3
     MeshHm4 = hm1 + hm2 + hm3 + hm4

     ! 3rd Order Mesh Coef
     MeshCoefm0 =  1.0*( MeshHm2*MeshHm3*MeshHm4 + MeshHm1*MeshHm3*MeshHm4 + &
          MeshHm1*MeshHm2*MeshHm4 + MeshHm1*MeshHm2*MeshHm3)/&
          (MeshHm1*MeshHm2*MeshHm3*MeshHm4) 
     MeshCoefm1 = -1.0*( MeshHm2*MeshHm3*MeshHm4)/&
          (hm1*hm2*(hm2 + hm3)*(hm2 + hm3 + hm4))
     MeshCoefm2 =  1.0*( MeshHm1*MeshHm3*MeshHm4)/(MeshHm2*hm2*hm3*(hm3+hm4))
     MeshCoefm3 = -1.0*( MeshHm1*MeshHm2*MeshHm4)/(MeshHm3*(hm3+hm2)*hm3*hm4)
     MeshCoefm4 =  1.0*( MeshHm1*MeshHm2*MeshHm3)/&
          (MeshHm4*(hm2+hm3+hm4)*(hm3+hm4)*hm4)
           
     do iSpecies=1,nIonsAdvect

        if (UseImprovedIonAdvection) then
           n0 = alog(LogINS(iAlt  ,iSpecies))
           n1 = alog(LogINS(iAlt-1,iSpecies))
           n2 = alog(LogINS(iAlt-2,iSpecies))
           n3 = alog(LogINS(iAlt-3,iSpecies))
           n4 = alog(LogINS(iAlt-4,iSpecies))
           n5 = alog(LogINS(iAlt-5,iSpecies))
           if (isnan(n0)) write(*,*) 'n0 :',iAlt, LogINS(iAlt  ,iSpecies)
           if (isnan(n1)) write(*,*) 'n1 :',iAlt-1, LogINS(iAlt-1  ,iSpecies)
           if (isnan(n2)) write(*,*) 'n2 :',iAlt-2, LogINS(iAlt-2  ,iSpecies)
           if (isnan(n3)) write(*,*) 'n3 :',iAlt-3, LogINS(iAlt-3  ,iSpecies)
           if (isnan(n4)) write(*,*) 'n4 :',iAlt-4, LogINS(iAlt-4  ,iSpecies)
           if (isnan(n5)) write(*,*) 'n5 :',iAlt-5, LogINS(iAlt-5  ,iSpecies)
        else
           n0 = LogINS(iAlt  ,iSpecies)
           n1 = LogINS(iAlt-1,iSpecies)
           n2 = LogINS(iAlt-2,iSpecies)
           n3 = LogINS(iAlt-3,iSpecies)
           n4 = LogINS(iAlt-4,iSpecies)
           n5 = LogINS(iAlt-5,iSpecies)
        endif

        dn = MeshCoefm0*n1 + &  
             MeshCoefm1*n2 + &  
             MeshCoefm2*n3 + &  
             MeshCoefm3*n4 + &  
             MeshCoefm4*n5      

        if (tec < MinTEC .and. UsePlasmasphereBC) then
           ! do nothing
        else

           ! Limit the slope of the gradient to be negative, since the ion density
           ! should be decreasing at the top of the model.
           
           if (dn > -0.001*n0/dAlt_F(nAlts)) dn = -0.001*n0/dAlt_F(nAlts)
           n0 = n1 + dn*dAlt_F(iAlt)

           if (UseImprovedIonAdvection) then
              LogINS(iAlt,iSpecies) = exp(n0)
           else
              LogINS(iAlt,iSpecies) = n0
           endif

        endif
              
        if (isnan(LogINS(iAlt,1))) write(*,*) 'svbc ',iAlt,LogINS(iAlt,1), n0, dn, &
             n1, n2, n3, n4, n5, tec, MinTEC, UsePlasmasphereBC, &
             dAlt_F(iAlt)

     enddo

  enddo

  do iAlt = nAlts+1, nAlts+2
     SumRho = 0.0     
     do iSpecies=1,nSpecies
        SumRho  = SumRho  + &
             Mass(iSpecies)*exp(LogNS(iAlt,iSpecies))     
     enddo
     LogRho(iAlt) = alog(SumRho)
  enddo

end subroutine set_vertical_bcs

