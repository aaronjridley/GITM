!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
program timing_test

  use ModMpi
  implicit none

  integer:: i, n_iter=0, n_step=10, me_world, numprocs, error

  logical:: timing_version_on
  character (len=40) :: timing_version_name
  real   :: timing_version_number

  integer, parameter :: Real8_= selected_real_kind(12)
  real(Real8_), external:: timing_func_d

  !------------------------------------------------------------------------

  call MPI_INIT(error)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me_world, error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, error)

  if(me_world==0)then
     write(*,'(a)')'========================================================='
     write(*,'(a,i3,a)')&
          'TIMING_TEST is running on',numprocs,' processors'
     write(*,*)
     call timing_version(timing_version_on,timing_version_name,&
          timing_version_number)
     write(*,'(2a,f5.2,a,l2)') 'Testing module ',&
          timing_version_name(1:len_trim(timing_version_name))//' version ',&
          timing_version_number,' functional=',timing_version_on
     write(*,'(a)')'========================================================='
  end if

  ! Time processor 0 only
  ! Other processors can also call timing routines, but it has no effect
  if(me_world==0)then
     ! Test behaviour when using an incorrect parameter name
     call timing_param_put_i('TestWrongName',13,error)
     if(error/=0)write(*,*)'timing_param_put_i("TestWrongName",13,error)',&
          ' resulted in error with error=',error

     call timing_active(.true.)
     call timing_depth(3)
     call timing_step(0)
  end if

  call timing_start('timing_test')

  call sleep(0.1)
  call timing_start('initialize')
  call sleep(0.3)
  call timing_stop('initialize')

  call timing_show('initialize',1)

  ! Test calling a function from a processor with inactive timing
  if(me_world==1)write(*,*)'Inactive timing should return -1:',&
       timing_func_d('sum',1,'initialize','timing_test')

  if(me_world==0)then
     write(*,*)'============================================================'
     write(*,*)'               STARTING ITERATIONS'
     write(*,*)'============================================================'
  end if
  do
     if(n_iter>=10)EXIT

     n_iter=n_iter+1; n_step=n_step+1
     call timing_step(n_step)

     call timing_start('advance_expl')
     do i=1,10
        call timing_start('calc_gradients')
        call timing_start('apply_limiters')
        call sleep(0.01)
        call timing_stop('apply_limiters')
        call timing_stop('calc_gradients')
        call timing_start('calc_facevalues')
        call sleep(0.02)
        call timing_stop('calc_facevalues')
        call sleep(0.01) ! other stuff
     end do
     call timing_stop('advance_expl')
     if(n_iter==5)then
        call timing_tree(2,2)
        call timing_reset('#all',2)
     end if

     if(me_world==0)write(*,*)'speed:',&
          100/max(1.D-10,timing_func_d('sum',1,'advance_expl','timing_test')),&
          ' at n_step=',n_step
  end do
  if(me_world==0)then
     write(*,*)'============================================================'
     write(*,*)'               STOPPING ITERATIONS'
     write(*,*)'============================================================'
  end if
  call timing_start('calc_gradients')
  call timing_start('apply_limiters')
  call sleep(0.01)
  call timing_stop('apply_limiters')
  call timing_stop('calc_gradients')
  call timing_start('save_output')
  call sleep(0.2)
  call timing_stop('save_output')
  call sleep(0.3) ! other stuff
  call timing_stop('timing_test')

  if(me_world==0)then
     write(*,*)'============================================================'
     write(*,*)'               STOPPING CALCULATIONS'
     write(*,*)'============================================================'
  end if

  call timing_report_style('tree')
  call timing_report
  call timing_report_total

  call timing_tree(2,2)
  call timing_tree(3,2)

  call timing_show('calc_gradients',1)
  call timing_show('calc_gradients',2)
  call timing_show('calc_gradients',3)

  call timing_sort(1,-1,.true.)
  call timing_sort(2,4,.true.)
  call timing_sort(3,-1,.true.)
  call timing_sort(3,-1,.false.)

  call MPI_FINALIZE(error)

contains

  !=======================================================================
  subroutine sleep(len)
    implicit none
    
    real, intent(in):: len
    real(Real8_) :: time_before
    !---------------------------------------------------------------------
    time_before=MPI_WTIME()
    do
       if(MPI_WTIME()>time_before+len) EXIT
    end do
  end subroutine sleep
  
end program timing_test

