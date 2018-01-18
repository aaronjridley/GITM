!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModNumConst
  use ModKind
  implicit none

  real, parameter:: &
       cZero      =  0.0, & 
       cOne       =  1.0, & 
       cHalf      =  0.5, & 
       cTiny      =  1.0E-6,                                &
       cHuge      =  1.0E18,                                &
       cSqrtTwo   =  1.4142135623730951,                    &
       cSqrtHalf  =  0.5*cSqrtTwo,                          &
       cPi        =  3.1415926535897932384626433832795,     &
       cTwoPi     =  2*cPi,                                 &
       cHalfPi    =  0.5*cPi,                               &
       cRadToDeg  =  180.0/cPi,                             &
       cDegToRad  =  cPi/180.0,                             &
       cTolerance =  1e-10

  real(Real8_), parameter:: &
       cTiny8     =  1e-10,                              &
       cPi8       =  3.1415926535897932384626433832795,  &
       cTwoPi8    =  2*cPi8

  ! integer unit matrix
  integer, parameter, dimension(3,3) :: i_DD = reshape( &
       (/1,0,0, 0,1,0, 0,0,1/), (/3,3/))

  ! real unit matrix (also Kronecker delta)
  real, parameter, dimension(3,3) :: cUnit_DD = reshape( &
       (/1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0/),&
       (/3,3/))

  ! Levi-Civita tensor
  integer, parameter :: iLeviCivita_III(3,3,3) = reshape( (/ &
       0, 0, 0,   0, 0,-1,  0, 1, 0,   &
       0, 0, 1,   0, 0, 0, -1, 0, 0,   &
       0,-1, 0,   1, 0, 0,  0, 0, 0 /), (/3,3,3/) )
  
end module ModNumConst
