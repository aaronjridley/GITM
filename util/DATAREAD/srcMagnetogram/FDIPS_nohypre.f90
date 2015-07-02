!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHypre

  implicit none

  public:: hypre_initialize
  public:: hypre_solver
  public:: hypre_preconditioner

contains

  !==========================================================================

  subroutine hypre_initialize

    call CON_stop('FDIPS is compiled without HYPRE')

  end subroutine hypre_initialize

  !==========================================================================

  subroutine hypre_solver

    call CON_stop('FDIPS is compiled without HYPRE')

  end subroutine hypre_solver

  !===========================================================================

  subroutine hypre_preconditioner(n, y_I)

    integer, intent(in):: n
    real, intent(inout):: y_I(n)

    call CON_stop('FDIPS is compiled without HYPRE')

  end subroutine hypre_preconditioner

end module ModHypre
!============================================================
subroutine read_hypre_param

  call CON_stop('FDIPS is compiled without HYPRE')

end subroutine read_hypre_param
