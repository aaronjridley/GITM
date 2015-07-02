!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program dipole11

  ! Create a magnetogram with spherical harmonics (1,1)
  ! proportional to sin(Theta)*cos(Phi)

  use ModNumConst, ONLY: cPi, cTwoPi, cRadToDeg
  use ModPlotFile, ONLY: save_plot_file

  implicit none

  real,parameter::Amplitude = 5.0
  integer, parameter:: nPhi = 360, nTheta = 180
  real, parameter:: dPhi= cTwoPi/nPhi, dSinTheta = 2.0/nTheta

  real:: Phi_I(nPhi), Theta_I(nTheta), Br_II(nPhi, nTheta)

  real :: Theta
  integer:: iPhi, iTheta
  !---------------------------------------------------------------------------

  ! Uniform in Phi
  do iPhi = 1, nPhi
     Phi_I(iPhi) = dPhi*(iPhi - 1)
  end do

  ! Uniform in sin(theta) and start from Theta = Pi (!) and finish at 0
  do iTheta = 1, nTheta
     Theta_I(iTheta) = acos( (iTheta-0.5)*dSinTheta - 1.0 )
  end do

  ! Br = A*sin(theta)*cos(phi)
  do iTheta = 1, nTheta
     do iPhi = 1, nPhi
        ! NOTE: this is not optimized at all
        Br_II(iPhi,iTheta) = Amplitude * sin(Theta_I(iTheta)) * cos(Phi_I(iPhi))
     end do
  end do

  open(11,file='fitsfile.dat',status='replace')
  write(11,'(a)')'#CR'
  write(11,'(a)')'1200'
  write(11,'(a)')'#nMax'
  write(11,'(a)')'90'
  write(11,'(a)')'#ARRAYSIZE'
  write(11,*)nPhi
  write(11,*)nTheta
  write(11,'(a)')'#START'

  do iTheta = 1, nTheta
     do iPhi = 1, nPhi
        write(11,'(e15.6)') Br_II(iPhi,iTheta)
     end do
  end do
  close(11)

  call save_plot_file('dipole11.out', &
       StringHeaderIn = 'DIPOLE11 output: [deg] [G]', &
       NameVarIn = 'Longitude Latitude Br LongitudeShift', &
       ParamIn_I = (/ 0.0 /), &
       Coord1In_I = cRadToDeg*Phi_I, &
       Coord2In_I = 90.0 - cRadToDeg*Theta_I, &
       VarIn_II  = Br_II)

  ! Uniform latitude distribution and start from Theta = Pi (!) and finish at 0
  do iTheta = 1, nTheta
     Theta_I(iTheta) = cPi*(1 - (iTheta - 0.5)/nTheta)
  end do

  ! Br = A*sin(theta)*cos(phi)
  do iTheta = 1, nTheta
     do iPhi = 1, nPhi
        ! NOTE: this is not optimized at all
        Br_II(iPhi,iTheta) = Amplitude * sin(Theta_I(iTheta)) * cos(Phi_I(iPhi))
     end do
  end do

  call save_plot_file('dipole11uniform.out', &
       StringHeaderIn = 'DIPOLE11 output: [deg] [G]', &
       NameVarIn = 'Longitude Latitude Br LongitudeShift', &
       ParamIn_I = (/ 0.0 /), &
       Coord1In_I = cRadToDeg*Phi_I, &
       Coord2In_I = 90.0 - cRadToDeg*Theta_I, &
       VarIn_II  = Br_II)

end program dipole11
!=================================================================
subroutine CON_stop(String)
  character(LEN=*),intent(in):: String
  write(*,*) String
  stop
end subroutine CON_stop
