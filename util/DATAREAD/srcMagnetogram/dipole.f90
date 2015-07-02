!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program dipole

  ! Create a magnetogram with spherical harmonics (1,0)
  ! proportional to cos(Theta), a simple dipole field

  implicit none

  real,parameter::Amplitude = 5.0
  integer, parameter:: nPhi = 360, nTheta = 180
  real, parameter:: dSinTheta =  2.0/nTheta

  real :: Theta
  integer:: iPhi, iTheta
  !--------------------------------------------------------------

  open(11,file='fitsfile.dat',status='replace')
  write(11,'(a)')'#CR'
  write(11,'(a)')'1200'
  write(11,'(a)')'#nMax'
  write(11,'(a)')'90'
  write(11,'(a)')'#ARRAYSIZE'
  write(11,*)nPhi
  write(11,*)nTheta
  write(11,'(a)')'#START'

  do iTheta = 0, nTheta-1
     Theta = acos( (iTheta+0.5)*dSinTheta - 1.0 )
     do iPhi = 0, nPhi-1
        write(11,'(e15.6)')Amplitude * cos(Theta)
     end do
  end do
  close(11)
  stop

end program dipole
!=================================================================
subroutine CON_stop(String)
  character(LEN=*),intent(in):: String
  write(*,*) String
  stop
end subroutine CON_stop
