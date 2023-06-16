
!
! It seems that NCAR stopped supporting the old format of the IGRF
! coefs, so I hacked things together to work with the new coefs.
!
!   https://github.com/NCAR/apex_fortran/blob/master/apex.f90
! Certain variables from ModApex were grabbed from this site.
!

module ModApex

  real,parameter :: re = 6371.2,            & ! Reference radius of IGRF
                    eps = 1.e-5,            & ! Small number
                    fltnvrs = 298.257223563   ! Inverse flatness of geoid
  
  integer,parameter :: nmax=13 ! maximum degree of IGRF coefficients
  integer,parameter :: ncoef = nmax*nmax + 2*nmax + 1 ! 196
  real,dimension(ncoef) :: &
       gb, & ! Coefficients for magnetic field calculation
       gv    ! Coefficients for magnetic potential calculation

  real, parameter :: &
       dtr=0.0174532925199432957692369076847, & ! degrees to radians
       rtd=57.2957795130823208767981548147      ! radians to degrees

  real ::  & ! Formerly /APXDIPL/ and /DIPOLE/
       colat, & ! Geocentric colatitude of geomagnetic dipole north pole (deg)
       elon,  & ! East longitude of geomagnetic dipole north pole (deg)
       vp,    & ! Magnitude, in T.m, of dipole component of magnetic
       ! potential at geomagnetic pole and geocentric radius re
       ctp,stp  ! cosine and sine of colat
 
end module ModApex

