!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! MODULE ------------------------------------------------
!/

module ModVertical

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect, nSpeciesTotal, nIons

  implicit none

  real, dimension(-1:nAlts+2)   :: dAlt_F, InvDAlt_F, Altitude_G, Gravity_G
  real, dimension(-1:nAlts+2)   :: InvRadialDistance_C, dAlt_C,Cv_1D
  real, dimension(-1:nAlts+2)   :: LogRho, Temp,MeanMajorMass_1d,Gamma_1d
  real, dimension(-1:nAlts+2,3) :: Vel_GD
  real, dimension(-1:nAlts+2,3) :: IVel
  
  real, dimension(-1:nAlts+2)             :: NewLogRho, NewTemp
  real, dimension(-1:nAlts+2,3)           :: NewVel_GD
  real, dimension(-1:nAlts+2,nSpecies) :: NewLogNS, NewVertVel
!  real, dimension(-1:nAlts+2,nIonsAdvect) :: NewLogINS
!  real, dimension(-1:nAlts+2,nIonsAdvect) :: LogINS

  real, dimension(-1:nAlts+2,nIons) :: NewLogINS
  real, dimension(-1:nAlts+2,nIons) :: LogINS
!  real, dimension(1:2,nIons) :: OldLogINSc


  real, dimension(-1:nAlts+2, nSpecies) :: LogNS, VertVel
  real, dimension(nAlts, nSpecies) :: NDensityS_1D

  real, dimension(0:nAlts+1) :: cMax
  
  real :: SZAVertical
  real :: MLatVertical
 
  integer :: iLon1D, iLat1D, iBlock1D

  real :: Heating(nAlts)
  real, dimension(-1:nAlts+2) :: EddyCoef_1d
  real, dimension( 0:nAlts+1) :: ViscCoef_1d
  real, dimension( 0:nAlts+1,1:nSpecies) :: ViscCoefS_1d
  real :: KappaTemp_1d( 0:nAlts+1)

  real :: Centrifugal, Coriolis, Lat, Lon, dAltdlon_1D, dAltdLat_1D

end module ModVertical

