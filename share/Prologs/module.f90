!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! Use this for your f90 module definition.

!BOP
!
!MODULE: ModExample - short description of ModExample
!
!DESCRIPTION:
! Long description of ModExample using LaTex syntax.
! Remember that underscore, per cent and dollar signs have to be escaped
! like example\_routine and 10\%, use {\bf bold face} for emphasis, etc.

!INTERFACE:
module ModExample

  !USES:
  !use ModOther ! This is commented out so that we can compile the code

  implicit none

  save

  private ! except

  !PUBLIC TYPES:
  public :: ExampleType ! short description of ExampleType
  type ExampleType
    real :: ExampleField     
  end type ExampleType

  !PUBLIC MEMBER FUNCTIONS:

  public :: example_iroutine ! short description of example_iroutine

  !PUBLIC DATA MEMBERS:

  real, public, parameter :: ExampleParameter=3.1415

  !LOCAL VARIABLES:
  real :: VeryImportantVariable

  !REVISION HISTORY: 
  !
  !  04/04/04 My Name <myname@umich.edu> - description of change
  !
  !EOP ----------------------------------------------------------------------

  character(len=*),parameter :: NameMod='ModExample'

contains

  subroutine example_iroutine

  end subroutine example_iroutine

end module ModExample
