!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-*- F90 -*- so emacs thinks this is an f90 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NAME: AB_ARRAY_module
!
! PURPOSE: a module which implements some array wrapper structures for 
!          AB2D.
!
! HISTORY:
!  1/05/01 Robert Oehmke: created
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module AB_ARRAY_module
  implicit none

  type :: AB_ARRAY_REAL
     real, pointer :: array(:)
  end type AB_ARRAY_REAL

  type :: AB_ARRAY_INT_2D
     integer, pointer :: array(:,:)
  end type AB_ARRAY_INT_2D

  type :: AB_ARRAY_INT
     integer, pointer :: array(:)
  end type AB_ARRAY_INT


end module AB_ARRAY_module


