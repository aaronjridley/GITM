!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!QUOTE: \clearpage
!
!BOP
!
!QUOTE: \section{util/TIMING: Library for Timing and Profiling}
!
!MODULE: ModTiming - the variables used by the TIMING utility
!
!INTERFACE:
module ModTiming
  !DESCRIPTION:
  !
  ! This module contains the variables shared by a collection of subroutines 
  ! and functions which together form the TIMING library.
  ! This utility facilitates timing and profiling of Fortran 90 codes.
  ! See the user manual for usage.
  !
  !REVISION HISTORY:
  ! 05/11/2001 G. Toth <gtoth@umich.edu> - initial version of TIMING
  !            several bug fixes and improvements later on
  !
  !EOP
  implicit none
  save

  integer, parameter :: Real8_ = selected_real_kind(12)

  logical :: UseTiming=.false.

  integer, parameter :: maxtiming=1000, maxclock=3

  character (LEN=40), dimension(maxtiming) :: sa_name

  ! These two declarations were replaced to avoid a pgf90 compiler error:
  ! character (LEN=40), parameter :: spaces = repeat(' ',40)
  ! character (LEN=79), parameter :: sepline = repeat('-',79)
  character (LEN=40), parameter :: spaces = &
       '                                        '
  character (LEN=79), parameter :: sepline = &
       '-------------------------------------------------------------------------------'

  real(Real8_), dimension(maxtiming):: da_start=0.0, da_last=0.0, da_sum_other
  real(Real8_), dimension(maxtiming,maxclock) :: da_sum=0.0

  integer, dimension(maxtiming) :: ia_step=-1, ia_depth=0, ia_parent=1
  integer, dimension(maxtiming,maxclock)   :: ia_call=0,  ia_iter=0

  logical, dimension(maxtiming) :: la_active=.false.

  integer :: ntiming=0, i_last=1, step_reset(maxclock)=-1, &
             current_depth=0, max_depth=-1

  character (LEN=10) :: report_style='cumu'

  integer :: step=-1

  integer :: lVerbose=-1

  integer :: iProc=0

  integer :: iUnit=6

  character (LEN=2) :: NameComp='  '

  interface
     function timing_cpu()
       implicit none
       real(selected_real_kind(12)) :: timing_cpu
     end function timing_cpu
  end interface

end module ModTiming
