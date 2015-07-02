!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! Module with the data for all plastics, used in the CRASH targets or used 
! to simulate them.
module CRASH_ModPolyimide

  ! Constants for polyimide

  use CRASH_ModAtomicMass

  implicit none
  save
  
  !Number of elements in C_22 H_10 N_2 O_5
  integer,parameter :: nPolyimide = 4 
  
  !Relative atomic concentrations in C_22 H_10 N_2 O_5
  real,dimension(nPolyimide),parameter :: CPolyimide_I = &
         (/22.0, 10.0, 2.0, 5.0/)/(22.0 + 10.0 + 2.0 +5.0)

  !Atomic number for elements:
  integer,dimension(nPolyimide),parameter :: nZPolyimide_I = &
         (/6, 1, 7, 8/)
   
  ! Averaged atomic mass
  real,parameter :: cAPolyimide = &
         (cAtomicMass_I(6) * 22.0 + &
         cAtomicMass_I(1) * 10.0 + &
         cAtomicMass_I(7) *  2.0 + &
         cAtomicMass_I(8) *  5.0) / &
         (22.0 + 10.0 + 2.0 +5.0)

end module CRASH_ModPolyimide

!=======================
module CRASH_ModAcrylic
  use CRASH_ModAtomicMass
  implicit none
  SAVE
  
  !Number of elements in C_5 O_2 H_8
  integer,parameter :: nAcrylic = 3 
  
  !Relative atomic concentrations in C_5 O_2 H_8
  real,dimension(nAcrylic),parameter :: cAcrylic_I = &
         (/5.0/15, 2.0/15, 8.0/15/)

  !Atomic number for elements:
  integer,dimension(nAcrylic),parameter :: nZAcrylic_I = &
         (/6, 8, 1/)
   
  ! Averaged atomic mass
  real,parameter :: cAAcrylic = &
         (cAtomicMass_I(6) * 5.0 + &
          cAtomicMass_I(8) * 2.0 + &
          cAtomicMass_I(1) * 8.0 ) / 15
end module CRASH_ModAcrylic
!==============================

module CRASH_ModPolysterene
  use CRASH_ModAtomicMass
  implicit none
  SAVE
  
  !Number of elements in C_1 H_1
  integer,parameter :: nPolysterene = 2 
  
  !Relative atomic concentrations in C_1 H_1
  real,dimension(nPolysterene),parameter :: cPolysterene_I = &
         (/0.5, 0.5/)

  !Atomic number for elements:
  integer,dimension(nPolysterene),parameter :: nZPolysterene_I = &
         (/6, 1/)
   
  ! Averaged atomic mass
  real,parameter :: cAPolysterene = &
         (cAtomicMass_I(6) + &
         cAtomicMass_I(1)) / 2.0
end module CRASH_ModPolysterene
!===============================
