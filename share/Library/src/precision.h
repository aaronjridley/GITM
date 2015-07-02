!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!  The integer iRealPrec is set to unity while compiling, if the length
!  of 1.0 is 8 bytes or longer (=2 or more words). 
!  It is set to zero if the length of 1.0 is 4 bytes (=1 word)
!
!  If the header precision.h is included into f90 or f77 file, then the
!  value of iRealPrec is set while compiling this file. It can be 
!  different in the different files if they are compiled with different
!  compiler flags.
! 
!  If the header precision.h is uncluded into the f90 module, then the 
!  value of iRealPrec is set while compiling this module. The value will
!  be the same in all files which use this module 
      INTEGER iRealPrec
      PARAMETER (iRealPrec=(1.00000000011 - 1.0)*10000000000.0)
