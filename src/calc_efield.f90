!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_efield(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock

  integer :: i, j, k, imax, jmax, kmax, ku, kd
  real :: maxi, edotb
  real :: altu, alt, altd, du, dd, r, r2, pu, p, pd, dpda
  real :: bdir(3)

  call report("Electric Field",2)
  call start_timing("calc_efield")

  EField = 0.0
  ExB    = 0.0

  maxi = 0.0

  !!! This is only first order accurate for stretched grids ???
  do k=0,nAlts+1
     do i=0,nLats+1
        do j=0,nLons+1

           EField(j,i,k,iEast_) = &
               -(Potential(j+1,i,k,iBlock)-Potential(j-1,i,k,iBlock))/ &
               ((Longitude(j+1,iBlock) - Longitude(j-1,iBlock)) * &
               0.5*(RadialDistance_GB(j+1, i, k, iBlock) &
               +    RadialDistance_GB(j-1, i, k, iBlock)) &
               *max(abs(CosLatitude(i,iBlock)),0.17)) ! 0.17 is 80 degrees.

           EField(j,i,k,iNorth_)= &
                -(Potential(j,i+1,k,iBlock)-Potential(j,i-1,k,iBlock)) &
                *0.5*InvdLatDist_GB(j,i,k,iBlock)

           ! vertical direction, take a larger stencil:  
           ku = k+7
           kd = k-7

           if (ku > nAlts+2) ku = nAlts+2
           if (kd < -1) kd = -1
              
           altu = Altitude_GB(j,i,ku,iBlock)
           alt  = Altitude_GB(j,i,k ,iBlock)
           altd = Altitude_GB(j,i,kd,iBlock)
           du   = altu - alt
           dd   = alt - altd

           r = du/dd
           r2 = r**2

           pu = potential(j,i,ku,iBlock)
           p  = potential(j,i,k ,iBlock)
           pd = potential(j,i,kd,iBlock)

           dpda = (pu - r2*pd - (1-r2)*p) / (du + r2*dd)

           EField(j,i,k,iUp_) = -dpda

!           ! the electric field in the vertical direction should not be very
!           ! large.  Let's limit it, since it seems to blow up:
!
!           if (EField(j,i,k,iUp_) >  MaxEField) EField(j,i,k,iUp_) =  MaxEField
!           if (EField(j,i,k,iUp_) < -MaxEField) EField(j,i,k,iUp_) = -MaxEField

           bdir = B0(j,i,k,1:3,iBlock) / B0(j,i,k,iMag_,iBlock)

           edotb = &
                efield(j,i,k,iEast_) * bdir(iEast_) + &
                efield(j,i,k,iNorth_) * bdir(iNorth_) + &
                efield(j,i,k,iUp_) * bdir(iUp_)
                
           efield(j,i,k,iEast_)  = efield(j,i,k,iEast_)  - edotb * bdir(iEast_)
           efield(j,i,k,iNorth_) = efield(j,i,k,iNorth_) - edotb * bdir(iNorth_)
           efield(j,i,k,iUp_)    = efield(j,i,k,iUp_)    - edotb * bdir(iUp_)

        enddo
     enddo
  enddo

  ExB(:,:,:,iEast_)  =    EField(:,:,:,iNorth_) * B0(:,:,:,iUp_,iBlock)    - &
                          EField(:,:,:,iUp_)    * B0(:,:,:,iNorth_,iBlock)

  ExB(:,:,:,iNorth_) = - (EField(:,:,:,iEast_)  * B0(:,:,:,iUp_,iBlock)    - &
                        EField(:,:,:,iUp_)    * B0(:,:,:,iEast_,iBlock))

  ExB(:,:,:,iUp_)    =    EField(:,:,:,iEast_)  * B0(:,:,:,iNorth_,iBlock) - &
                        EField(:,:,:,iNorth_) * B0(:,:,:,iEast_,iBlock)

  ExB(:,:,:,iEast_)  = ExB(:,:,:,iEast_)  / (B0(:,:,:,iMag_,iBlock)**2)
  ExB(:,:,:,iNorth_) = ExB(:,:,:,iNorth_) / (B0(:,:,:,iMag_,iBlock)**2)
  ExB(:,:,:,iUp_)    = ExB(:,:,:,iUp_)    / (B0(:,:,:,iMag_,iBlock)**2)

  if (iDebugLevel > 5) then

     kmax = 1
     jmax = 1
     imax = 1
     maxi = 0.0

     do k=0,nAlts+1
        do i=0,nLats+1
           do j=0,nLons+1

              if (abs(ExB(j,i,k,iNorth_)) > maxi) then
                 write(*,*) "======> efield : ", i,j,k, &
                      EField(j,i,k,iEast_) * 1000.0, &
                      Potential(j+1,i,k,iBlock)/1000.0, &
                      Potential(j-1,i,k,iBlock)/1000.0, &
                      dLonDist_GB(j, i, k, iBlock), &
                      ExB(j,i,k,iNorth_), B0(j,i,k,iMag_,iBlock)
                 maxi = abs(ExB(j,i,k,iNorth_))
              endif

           enddo
        enddo
     enddo

  endif

  if (iDebugLevel > 3) then

     write(*,*) "====> Min ExB : ", minval(ExB(:,:,:,iEast_)),  &
          minval(ExB(:,:,:,iNorth_)),  minval(ExB(:,:,:,iUp_))
     write(*,*) "====> Max ExB : ", maxval(ExB(:,:,:,iEast_)),  &
          maxval(ExB(:,:,:,iNorth_)),  maxval(ExB(:,:,:,iUp_))

  endif

  call end_timing("calc_efield")

end subroutine calc_efield
