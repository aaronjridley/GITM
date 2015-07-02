!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP -------------------------------------------------------------------
!
!MODULE: ModIoUnit - general utilities for Fortran I/O units.
!
!DESCRIPTION:
!
! This module provides various utilities related to Fortran I/O units.
! In particular independently developped components can use the 
! io\_unit\_new() function to obtain an unused IO unit for extended use.
! 
! The unit number in UNITTMP\_ is a safe unit number to open and close a file 
! if no other file is opened between the open and close and all programs
! use ModIoUnit to obtain unit numbers.
!
! Standard output has unit number STDOUT\_=6. This constant is easier to read.
!
! The io\_unit\_clean subroutine closes all open IO units and deletes the
! empty files.
!
! The methods in this module can be tested by running the 
! io\_unit\_test subroutine.
!
!INTERFACE:

module ModIoUnit

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: io_unit_new    ! Return an unused unit number for extended use
  public :: io_unit_clean  ! Close open units, delete empty files
  public :: io_unit_test   ! Test the functionality of this module

  !PUBLIC DATA MEMBERS:

  public :: UNITTMP_       ! For open read/write close without intervening open
  public :: STDOUT_        ! Fortran unit number for standard output

  integer, parameter :: STDOUT_       = 6     ! Standard output
  integer, parameter :: UNITTMP_      = 9     ! Temporary unit number

  !LOCAL VARIABLES:

  integer, parameter :: MinUnitNumber = 20    ! Smallest allowed unit number
  integer, parameter :: MaxUnitNumber = 1000  ! Largest allowed unit number

  integer :: iUnitMax = UNITTMP_              ! The largest unit number used

  !REVISION HISTORY:
  ! 01Aug03  Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  ! 20Aug04  Gabor Toth                     added debugging for io_unit_new
  !EOP ___________________________________________________________________

  character (len=*), parameter :: NameMod = 'ModIoUnit'

contains

  function io_unit_new()  result(iUnit)

    !  Returns a unit number of a unit that exists and is not connected
    integer :: iUnit
    logical :: IsExisting, IsOpened
    integer :: iError

    character (len=*), parameter :: NameSub = NameMod//'::io_unit_new'
    !--------------------------------------------------------------------

    do iUnit = MinUnitNumber, MaxUnitNumber
       inquire (&
            unit   = iUnit,       &
            exist  = IsExisting,  &
            opened = IsOpened,    &
            iostat = iError)
       if (IsExisting .and. .not. IsOpened .and. iError == 0) then
          iUnitMax = max(iUnitMax, iUnit)
          return
       end if
    end do

    iUnit = -1

  end function io_unit_new
  !===========================================================================
  subroutine io_unit_clean

    ! Close all open units for this processor
    integer :: iUnit, iError
    logical :: IsOpen
    character(len=100) :: Name
    character :: String
    !------------------------------------------------------------------------
    do iUnit = UNITTMP_,iUnitMax

       inquire(iUnit,OPENED=IsOpen,NAME=Name)
       if(IsOpen)then
          ! Clos file so that output is flushed
          close(iUnit)
          ! Try to open file and read 1 character
          open(iUnit,FILE=Name,STATUS='old',IOSTAT=iError)
          if(iError/=0) CYCLE
          read(iUnit,'(a1)',IOSTAT=iError)String
          if(iError<0)then
             ! Delete empty files
             close(iUnit,STATUS='delete')
          else
             ! Close file again
             close(iUnit)
          end if
       end if
    end do

  end subroutine io_unit_clean
  !==========================================================================
  subroutine io_unit_test

    integer :: iUnit1, iUnit2, iUnit3, iUnit4
    logical :: IsExisting
    !---------------------------------------------------------------------

    write(*,'(a)')'Testing io_unit_new()'
    iUnit1 = io_unit_new()
    if(iUnit1/=MinUnitNumber)write(*,*)'test io_unit_new() failed: ',&
         'iUnit1=',iUnit1,' should be MinUnitNumber=',MinUnitNumber
    open(iUnit1,file='ascii',status='unknown',form='formatted')
    write(iUnit1,*)1

    iUnit2 = io_unit_new()
    if(iUnit2/=MinUnitNumber+1)write(*,*)'test io_unit_new() failed: ',&
         'iUnit2=',iUnit2,' should be MinUnitNumber+1=',MinUnitNumber+1
    open(iUnit2,file='binary',status='unknown',form='unformatted')
    write(iUnit2)1

    iUnit3 = io_unit_new()
    if(iUnit3/=MinUnitNumber+2)write(*,*)'test io_unit_new() failed: ',&
         'iUnit3=',iUnit3,' should be MinUnitNumber+2=',MinUnitNumber+2
    open(iUnit3,file='empty_ascii',status='unknown',form='formatted')

    iUnit4 = io_unit_new()
    if(iUnit4/=MinUnitNumber+3)write(*,*)'test io_unit_new() failed: ',&
         'iUnit4=',iUnit4,' should be MinUnitNumber+3=',MinUnitNumber+3
    open(iUnit4,file='empty_binary',status='unknown',form='unformatted')

    write(*,'(a)')'Testing io_unit_clen'
    call io_unit_clean

    inquire(file='ascii',exist=IsExisting)
    if(.not.IsExisting)then
       write(*,*)'test io_unit_clean failed: ',&
            'file "ascii" should not have been deleted'
    else
       open(iUnit1,file='ascii',status='unknown',form='formatted')
       close(iUnit1,STATUS='delete')
       inquire(file='ascii',exist=IsExisting)
       if(IsExisting)write(*,*)'failed to delete file "ascii"'
    end if

    inquire(file='binary',exist=IsExisting)
    if(.not.IsExisting)then
       write(*,*)'test io_unit_clean failed: ',&
            'file "binary" should not have been deleted'
    else
       open(iUnit2,file='binary',status='unknown',form='unformatted')
       close(iUnit2,STATUS='delete')
       inquire(file='binary',exist=IsExisting)
       if(IsExisting)write(*,*)'failed to delete file "binary"'
    end if

    inquire(file='empty_ascii',exist=IsExisting)
    if(IsExisting)then
       write(*,*)'test io_unit_clean failed: ',&
            'file "empty_ascii" should have been deleted'
       open(iUnit3,file='empty_ascii',status='unknown',form='formatted')
       close(iUnit3,STATUS='delete')
       inquire(file='empty_ascii',exist=IsExisting)
       if(IsExisting)write(*,*)'failed to delete file "empty_ascii"'
    end if

    inquire(file='empty_binary',exist=IsExisting)
    if(IsExisting)then
       write(*,*)'test io_unit_clean failed: ',&
            'file "empty_binary" should have been deleted'
       open(iUnit4,file='empty_binary',status='unknown',form='unformatted')
       close(iUnit4,STATUS='delete')
       inquire(file='empty_binary',exist=IsExisting)
       if(IsExisting)write(*,*)'failed to delete file "empty_binary"'
    end if

  end subroutine io_unit_test

end module ModIoUnit
