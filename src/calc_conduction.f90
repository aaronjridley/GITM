! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

! -------------------------------------------------------------------------------------
! This file contains the following subroutines:
! - calc_conduction: calculates the thermal conduction across all lats/lons
! - calc_conduction_1d: a generalized solver for the Poisson Equation
! -------------------------------------------------------------------------------------

! -------------------------------------------------------------------------------------
! This calculates the thermal conduction in the vertical direction
! JMB Update 07/13/2017
! Updated Conduction Routines Require intermediate
! Steps outlined in Harwood et al. [2016]
! Oscillation-free method for semi-linear diffusion equations
! digitalcommons.georgefox.edu  
! -------------------------------------------------------------------------------------

subroutine calc_thermal_conduction(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  implicit none

  integer, intent(in) :: iBlock
  integer :: iErr
  
  real :: Rho110(nLons, nLats, 0:nAlts+1)
  real :: tmp2(nLons, nLats, 0:nAlts+1)

  ! NeuBCS used for Conduction and Bulk Winds
  logical ::NeuBCS                      ! Set to true if you want Neumann BCs

  real :: Prandtl(nLons,nLats,0:nalts+1)

  ! Sub-timestep used in the multi-step conduction
  real :: DtCSLocal

  ! Temperature variables
  real :: TemperatureStage1(1:nLons,1:nLats,-1:nAlts+2)
  real :: TemperatureH(1:nLons,1:nLats,-1:nAlts+2)
  real :: TemperatureF(1:nLons,1:nLats,-1:nAlts+2)

  real :: TmpTemp(nLons, nLats, -1:nAlts+2)
  real :: TmpDiff(nLons, nLats, 0:nAlts+1)
  real :: TmpMulFac(nLons, nLats, 0:nAlts+1)

  if (iDebugLevel > 4) write(*,*) "=====> conduction", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iErr)

  Rho110 = Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock)
  
  ! Note: Neumann = .true. is needed if you want 
  ! the gradient to be zero at the top.

  NeuBCS = .true.  ! Use Neumann Boundary Conditions

  tmp2 = Rho110 * cp(1:nLons, 1:nLats,0:nAlts+1, iBlock)

  if (UseTurbulentCond) then
     Prandtl = &
          KappaEddyDiffusion(1:nLons, 1:nLats,0:nAlts+1, iBlock) * &
          Rho110 * &
          Cp(1:nLons, 1:nLats,0:nAlts+1, iBlock)
  else 
     Prandtl = 0.0
  endif

  DtCSLocal = Dt/2.0
  TmpTemp = Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
       TempUnit(1:nLons, 1:nLats,-1:nAlts+2)
  TmpDiff = KappaTemp(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
       Prandtl(1:nLons, 1:nLats, 0:nAlts+1)

  call calc_conduction_1d(iBlock, DtCSLocal, NeuBCS, &
       TmpTemp, TmpDiff, tmp2, MoleConduction)

  TemperatureStage1(1:nLons,1:nLats,0:nAlts+1) = &
       Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock) + &
       MoleConduction(1:nLons,1:nLats,0:nAlts+1)/&
       TempUnit(1:nLons,1:nLats,0:nAlts+1)

  DtCSLocal = Dt/2.0
  TmpTemp = TemperatureStage1(1:nLons, 1:nLats,-1:nAlts+2) * &
       TempUnit(1:nLons, 1:nLats,-1:nAlts+2)
  TmpDiff = KappaTemp(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
       Prandtl(1:nLons, 1:nLats, 0:nAlts+1)

  call calc_conduction_1d(iBlock, DtCSLocal, NeuBCS, &
       TmpTemp, TmpDiff, tmp2, MoleConduction)

  TemperatureH(1:nLons,1:nLats, 0:nAlts+1) = &
       TemperatureStage1(1:nLons,1:nLats,0:nAlts+1) + &
       MoleConduction(1:nLons,1:nLats,0:nAlts+1)/&
       TempUnit(1:nLons,1:nLats,0:nAlts+1)

  ! Full Time Step Update
  DtCSLocal = Dt
  TmpTemp = Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
       TempUnit(1:nLons, 1:nLats,-1:nAlts+2)
  TmpDiff = KappaTemp(1:nLons, 1:nLats, 0:nAlts+1, iBlock) + &
       Prandtl(1:nLons, 1:nLats, 0:nAlts+1)

  call calc_conduction_1d(iBlock, DtCSLocal, NeuBCS, &
       TmpTemp, TmpDiff, tmp2, MoleConduction)

  TemperatureF(1:nLons,1:nLats, 0:nAlts+1) = &
       Temperature(1:nLons,1:nLats, 0:nAlts+1,iBlock) + &
       MoleConduction(1:nLons,1:nLats, 0:nAlts+1) / &
       TempUnit(1:nLons,1:nLats, 0:nAlts+1)

  ! Note that Conduction is the net temperature update
  Conduction(1:nLons,1:nLats,0:nAlts+1) = &
       (2.0*TemperatureH(1:nLons,1:nLats,0:nAlts+1) - &
       TemperatureF(1:nLons,1:nLats,0:nAlts+1)) - &
       Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock) 

end subroutine calc_thermal_conduction


! -------------------------------------------------------------------------------------
! General Poisson Equation solver in 1d (vertical), but looping across lons / lats
! -------------------------------------------------------------------------------------

subroutine calc_conduction_1d(iBlock, DtIn, NeuBCS, Quantity, Diff, MulFac, dTdt_cond)

  use ModSizeGitm
  use ModGITM, only: dAlt_GB, Latitude, Longitude, Altitude_GB, RadialDistance_GB
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock
  real, intent(in)    :: DtIn
  logical, intent(in) :: NeuBCS
  real, intent(in)    :: Quantity( nLons, nLats, -1:nAlts+2)
  real, intent(in)    :: Diff(     nLons, nLats,  0:nAlts+1)
  real, intent(in)    :: MulFac(   nLons, nLats,  0:nAlts+1)
  real, intent(out)   :: dTdt_cond(nLons, nLats,  0:nAlts+1)

  real, dimension(0:nAlts+1) :: &
       m, du, r, du12, du22, dl, lou, dlou, di, r2

  real, dimension(0:nAlts+1) :: tempold, temp
  real, dimension(0:nAlts+1) :: a, b, c, d, cp, dp
  integer :: iLon, iLat, iAlt
  integer :: i

  call start_timing("conduction_1d")
  call report("calc_conduction_1d",3)

  do iLon = 1, nLons
     do iLat = 1, nLats

        tempold = Quantity(iLon, iLat, 0:nAlts+1)

        r2 = RadialDistance_GB(iLon,iLat, 0:nAlts+1,iBlock)**2
        di = diff(iLon,iLat,:)*r2

        m = DtIn/(MulFac(iLon,iLat, 0:nAlts+1)*r2)
        du = Altitude_GB(iLon,iLat, 1:nAlts+2,iBlock) - &
             Altitude_GB(iLon,iLat, 0:nAlts+1,iBlock)
        dl = Altitude_GB(iLon,iLat, 0:nAlts+1,iBlock) - &
             Altitude_GB(iLon,iLat,-1:nAlts+0,iBlock)
        r = du/dl

        du12 = du*du * (1+r)*(1+r)
        du22 = 0.5 * (dl*du + du*du)

        lou  = di/du22
        dlou = di/du22

        dl(1:nAlts) =                            di(2:nAlts+1) - &
                         r(1:nAlts)*r(1:nAlts) * di(0:nAlts-1) - &
                      (1-r(1:nAlts)*r(1:nAlts))* di(1:nAlts  )

        dl(0) = dl(1)
        dl(nAlts+1) = dl(nAlts)

        ! Do a google search for a tri-diagnal solver and you will come up
        ! with this:

        a =  di/du22*r - dl/du12 * r*r

        c =  di/du22 + dl/du12

        b = -1/m - di/du22*(1+r) - dl/du12*(1-r*r) 

        d = -tempold/m

        ! Boundary Conditions:

        d(0) = -tempold(0)
        a(0) = 0.0
        b(0) = -1.0
        c(0) = 0.0

!        i = nAlts+1
!        a(i) = -1.0*( -r(i)*(1.0+r(i))*di(i)/du22(i))
!        b(i) = -1.0*( 1.0/m(i) + r(i)*(1+r(i))*di(i)/du22(i))
!        c(i) = 0.0
!        d(i) = -tempold(i)/m(i)

        if(NeuBCS) then
           i = nAlts+1 
           !a(nAlts+1) = 1.0
           !b(nAlts+1) = -1.0
           !c(nAlts+1) = 0.0
           !d(nAlts+1) = 0.0
           a(i) =  1.0*( r(i)*(1.0+r(i))*di(i)*m(i)/du22(i))
           b(i) = -1.0*( 1.0 + r(i)*(1+r(i))*di(i)*m(i)/du22(i))
           c(i) = 0.0
           d(i) = -tempold(i)
        else
           a(nAlts+1) = 0.0
           b(nAlts+1) = -1.0
           c(nAlts+1) = 0.0
           d(nAlts+1) = -tempold(nAlts+1)
        endif

        cp(0) = c(0)/b(0)
        do iAlt = 1, nAlts+1
           cp(iAlt) = c(iAlt)/(b(iAlt)-cp(iAlt-1)*a(iAlt))
        enddo
        dp(0) = d(0)/b(0)
        do iAlt = 1, nAlts+1
           dp(iAlt) = (d(iAlt)-dp(iAlt-1)*a(iAlt))/(b(iAlt)-cp(iAlt-1)*a(iAlt))
        enddo
        temp(nAlts+1) = dp(nAlts+1)
        do iAlt=nAlts,0,-1
           temp(iAlt) = dp(iAlt)-cp(iAlt)*temp(iAlt+1)
        enddo

        dTdt_cond(iLon,iLat,0:nAlts+1) = &
             temp(0:nAlts+1) - Quantity(iLon,iLat,0:nAlts+1)

     enddo
  enddo

  call end_timing("conduction_1d")

end subroutine calc_conduction_1d

