!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModAtomicMass
 implicit none
 integer,parameter::nZMax=79
 integer, parameter :: nMixMax = 6

 !Atomic masses of the first 54 elements
 !Given in atomic mass units (1 atomic mass unit = 1.66053886 × 10e-27 kilograms)
 real,parameter,dimension(nZMax) :: cAtomicMass_I=(/&
  !    1    !     2    !      3    !     4    !     5     !    6   !     7    !    8    !      9     !    10   !
  1.007947  , 4.0026022, 6.9412    , 9.0121823, 10.8115   , 12.0111, 14.006747, 15.99943, 18.99840329, 20.17976, &
  22.9897686, 24.30506 , 26.9815395, 28.08553 , 30.9737624, 32.0666, 35.45279 , 39.9481 , 39.09831   , 40.0789 , & ! 10-s
  44.9559109, 47.883   , 50.94151  , 51.99616 , 54.938051 , 55.8473, 58.933201, 58.6934 , 63.5463    , 65.392  , & ! 20-s
  69.7231   , 72.612   , 74.921592 , 78.963   , 79.904    , 83.801 , 85.46783 , 87.621  , 88.905852  , 91.2242 , & ! 30-s
  92.906382 , 95.941   , 98.0      , 101.072  , 102.905503, 106.421, 107.86822, 112.4118, 114.821    , 118.7107, & ! 40-s
  121.757   , 127.603  , 126.904473, 131.292  , -1.0      , -1.0   , -1.0     ,  -1.0   , -1.0       , -1.0    , & ! 50-s
  -1.0      , -1.0     , -1.0      , -1.0     , -1.0      , -1.0   , -1.0     ,  -1.0   , -1.0       , -1.0    , & ! 60-s
  -1.0      , -1.0     , -1.0      , -1.0     , -1.0      , -1.0   , -1.0     ,  -1.0   , 196.96656  /) ! 70-s
end module CRASH_ModAtomicMass
