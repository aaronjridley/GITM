!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!
!MODULE: example_program - the main executable
!
!DESCRIPTION:
!
! The main executable is documented as if it was a module, because
! Protex does not provide a !PROGRAM: tag. Use LaTex syntax to 
! describe the purpose of this program.
!
!INTERFACE:
program example_program

  !USES:
  use ModExample, ONLY: ExampleParameter
  use ModExample2, ONLY: example_iroutine

  implicit none

  !LOCAL VARIABLES:
  character (len=*), parameter :: NameProgram='example_program'
  real :: ImportantVar1, ImportantVar2

  !REVISION HISTORY:
  ! 04/27/2004 My Name <myemail@umich.edu> - initial version
  !EOP
  ! other local variables not worth documenting
  !...
  !---------------------------------------------------------------------------
  !BOC
  write(*,*)NameProgram,': documented executable statements come here'
  call example_routine(1.0, ImportantVar1, ImportantVar2)
  call example_iroutine(2.0, ImportantVar1, ImportantVar2)
  !EOC
  write(*,*)NameProgram,': this part should not appear in the documentation'

end program example_program
