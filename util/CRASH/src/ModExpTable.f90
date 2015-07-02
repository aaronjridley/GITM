!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
 !-------
 ! \
module CRASH_M_expTab
  ! /
  !  tabulated exponential, 10 times faster that hard wired exp()
  !-------
  implicit none
  integer,private,save :: nb_Exp=-1	! flag of preparation

  !Number of tabulated points
  integer,parameter :: ex_nb=15000

  !Step in the argument; maximal argument!
  real,parameter :: ex_du=1e-2, ex_max=ex_nb*ex_du

  !Inverse of \Delta u
  real,parameter :: ex_sdu = 1.0/ex_du

  !To be used in interpolating the function
  real ::ex_c,ex_cm,ex_tab(0:ex_nb)

  !Auxiliary reals
  real :: ex_u,ex_y

  !Auxiliary integer
  integer :: ex_i
  !-------
contains
  !===================
  subroutine exp_tab8()
    !Initializes the table for exponential function 
    if(nb_exp.eq.ex_nb) return
    nb_exp = ex_nb

    do ex_i = 0,ex_nb
       ex_u = ex_i*ex_du
       ex_tab(ex_i) = exp(-ex_u)
    end do
    !Combination: (1 - exp(-\Delta u) )/\Delta u 
    ex_c=(1.0 - ex_tab(1))*ex_sdu
    !Combination: (exp(\Delta u) - 1  )/\Delta u
    ex_cm=(1.0/ex_tab(1) - 1.0)*ex_sdu
  end subroutine exp_tab8
 !-------
 ! \
end module CRASH_M_expTab
 ! /
 !-------
 
