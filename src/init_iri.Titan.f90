!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine init_iri

  use EUA_ModIri90, ONLY: iri90

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModTime

  implicit none

  ! iri variables

  integer :: jmag,nzkm
  real, dimension(1:11, 1:nAlts) :: outf
  real, dimension(1:30) :: oarr
  logical, dimension(1:12) :: jf
  !  character, dimension(1:80) :: ccirnm, ursinm

  integer :: iBlock, iAlt, iLat, iLon, iIon
  real :: geo_lat, geo_lon, geo_alt(1), geo_lst

  call report("init_iri",0)

! There is no equivalent empirical model for Titan's ionosphere

end subroutine init_iri
