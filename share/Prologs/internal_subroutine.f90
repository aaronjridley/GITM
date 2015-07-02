!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModExample2

  implicit none
  private ! except
  public :: example_iroutine
  character(len=*), parameter:: NameMod = 'ModExample2'

contains

  !BOP ======================================================================
  !IROUTINE: example_iroutine - put short description here

  !INTERFACE:
  subroutine example_iroutine(InputVar, OutputVar, InputOutputVar)

    !USES:
    use ModExample, ONLY: ExampleParameter

    !INPUT ARGUMENTS: 
    real, intent(in) :: InputVar          ! short description of InputVar

    !OUTPUT ARGUMENTS:
    real, intent(out) :: OutputVar        ! short description of OutputVar

    !INPUT/OUTPUT ARGUMENTS: 
    real, intent(inout) :: InputOutputVar ! short description of InputOutputVar

    !DESCRIPTION: 
    ! This is a template for a subroutine that comes after a CONTAINS 
    ! statement, i.e. a routine internal to a module or another 
    ! subroutine/function.
    !
    ! Provide here long description of subroutine example\_routine 
    ! in Latex format.
    ! Do not repeat information already provided by the above tags.

    !LOCAL VARIABLES:
    real :: AnImportantLocalVariable

    !REVISION HISTORY: 
    ! 04/27/2004 My Name <myemail@umich.edu> - initial version
    ! 04/28/2004 My Name <myemail@umich.edu> - fixed some typos
    !EOP  

    ! local variables not worth of documenting come here

    character(len=*), parameter:: NameSub=NameMod//'::example_iroutine'

    !------------------------------------------------------------------------
    !BOC
    write(*,*) NameSub,': documented executable statements come here'
    !EOC

    write(*,*) NameSub,': this part should not appear in the documentation'

  end subroutine example_iroutine

end module ModExample2
