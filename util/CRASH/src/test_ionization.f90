!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!***********************************************************************
!    calculation of ionization equilibrium, single material
!    for given - concentration of  atomic particles Co+Ciz [cm^{-3}]
!              - electron temperature Te[eV]
!***********************************************************************
program saha
  use CRASH_ModStatSum
  use CRASH_ModPartition
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  implicit NONE
  real :: &
       Nao = 1.00e18,  &  ! cm-3
       vTe=5., NaTrial, vU 

  integer,parameter :: nN=5 , nT=25, nU=40
  real    :: dTe, dLogN, dU
  integer :: iT, nT1=1000000, iN, iU, iLoop

  real    :: Z_I(0:nN),Z2Out_I(0:nN),Uav_I(0:nN),Cv_I(0:nN), Te_I(0:nN)
  real    :: ZNoCoulomb_II(1:nT, 0:nN)
  integer :: iIter_I(0:nN)
  !character(LEN=*),parameter,dimension(0:nN) :: Separator_I='|'
  !character(LEN=*),parameter,dimension(0:nN) :: Separator1_I='/'
  integer::iError


  !-------------------------------------------
  LogGeMinFermi = -4.0

  dTe = 5.0; dLogN=log(10.0); dU=100.0

  UsePreviousTe = .true.

  call set_element( 54 )


  open(24,file='../doc/Table1.tex',status='replace')
  write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Te[eV]\textbackslash \textbackslash Na[$1/cm^3$] & $10^{18}$'//&
       ' & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(24,'(a)')'\hline' 
  write(24,'(a)')'\hline'


  do iT  = 1,nT
     if (((iT-1)/25)*25==(iT-1).and.iT>25) then

        write(24,'(a)')'\end{tabular}', char(10)
        !--------------
        write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Te[eV]\textbackslash \textbackslash Na[$1/cm^3$] & $10^{18}$ & $10^{19}$'//&
             ' & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(24,'(a)')'\hline' 
        write(24,'(a)')'\hline'
     end if
     vTe = dTe * iT

     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
        Z_I(iN) = z_averaged() 
        Z2Out_I(iN)= z2_averaged()/Z_I(iN)
        Uav_I(iN)=internal_energy()
        Cv_I(iN)=heat_capacity()
        if(iError/=0)then
           Z_I(iN) = -  Z_I(iN)
           Z2Out_I(iN)= - Z2Out_I(iN)
           Uav_I(iN)= -Uav_I(iN)
           Cv_I(iN)= - Cv_I(iN)
           write(*,*)'Error=',iError,vTe,NaTrial*1000000.0
        end if
     end do
     write(24,'(f5.0,6(a,f7.1),a)') vTe,&
          (' & ', Z_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(a)')'\hline'

     ZNoCoulomb_II(iT,:) = Z_I
  end do

  write(24,'(a)')'\end{tabular}'

  close(24)

!_______________________________________________
  open(24,file='../doc/Table2.tex',status='replace')
  write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Te[eV] & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ '//&
       '& $U_{av} | C_v$ & $U_{av} | C_v$\tabularnewline'
  write(24,'(a)')'\hline' 
  write(24,'(a)')'\hline'


  do iT  = 1,nT
     if (((iT-1)/25)*25==(iT-1).and.iT>25) then

        write(24,'(a)')'\end{tabular}', char(10)
        !--------------
        write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Te[eV] & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ '//&
             '& $U_{av} | C_v$ & $U_{av} | C_v$\tabularnewline'
        write(24,'(a)')'\hline' 
	write(24,'(a)')'\hline' 
     end if
     vTe = dTe * iT

     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
        Z_I(iN) = z_averaged() 
        Z2Out_I(iN)= z2_averaged()/Z_I(iN)
        Uav_I(iN)=internal_energy()
        Cv_I(iN)=heat_capacity()
        if(iError/=0)then
           Z_I(iN) = -  Z_I(iN)
           Z2Out_I(iN)= - Z2Out_I(iN)
           Uav_I(iN)= -Uav_I(iN)
           Cv_I(iN)= - Cv_I(iN)
           write(*,*)'Error=',iError,vTe,NaTrial*1000000.0
        end if
     end do
     write(24,'(f5.0,6(a,f8.1,a,f7.1),a)') vTe,&
          (' & ', Uav_I(iN), ' | ', Cv_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(a)')'\hline'


  end do

  write(24,'(a)')'\end{tabular}'

  close(24)

  !_____________________________________
  open(25,file='../doc/Table3.tex',status='replace')
  write(25,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(25,'(a)')'\hline'
  write(25,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(25,'(a)')'\hline'
  write(25,'(a)')'U[eV] & Te (Iterations) &  Te (Iterations) &  Te (Iterations) & '//&
       ' Te (Iterations) &  Te (Iterations) &  Te (Iterations) \tabularnewline'
  write(25,'(a)')'\hline' 
  write(25,'(a)')'\hline'
  do iU  = 1,nU
     if (((iU-1)/50)*50==(iU-1).and.iU>50) then
        write(25,'(a)')'\end{tabular}', char(10)
        !------------
        write(25,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(25,'(a)')'\hline'
        write(25,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(25,'(a)')'\hline'
        write(25,'(a)')'U[eV] & Te (Iterations) &  Te (Iterations) &  Te (Iterations) & '//&
             ' Te (Iterations) &  Te (Iterations) &  Te (Iterations) \tabularnewline'
        write(25,'(a)')'\hline' 
        write(25,'(a)')'\hline'
     end if
     vU = dU * iU 
     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_temperature(vU,NaTrial*1000000.0,iError)
        Te_I(iN) = Te 
        iIter_I(iN) = iIterTe
        if(iError/=0)then
           write(*,*)'Error=',iError,Te_I(iN),NaTrial*1000000.0
           Te_I(iN)=-Te_I(iN)
        end if
     end do

     write(25,'(f6.0,6(a,f7.1,a,i7,a),a)') vU,&
          (' & ', Te_I(iN), ' (', iIter_I(iN),')', iN=0,nN ),'\tabularnewline'
     write(25,'(a)')'\hline'

  end do


  write(25,'(a)')'\end{tabular}'
  close(25)
  !_____________________________________

  open(25,file='../doc/Table4.tex',status='replace')
  write(25,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(25,'(a)')'\hline'
  write(25,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(25,'(a)')'\hline'
  write(25,'(a)')'P/Na[eV] & Te (Iterations) &  Te (Iterations) &  Te (Iterations) & '//&
       ' Te (Iterations) &  Te (Iterations) &  Te (Iterations) \tabularnewline'
  write(25,'(a)')'\hline' 
  write(25,'(a)')'\hline'

  dU = 50.0
  do iU  = 1,nU
     if (((iU-1)/50)*50==(iU-1).and.iU>50) then
        write(25,'(a)')'\end{tabular}', char(10)
        !------------
        write(25,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(25,'(a)')'\hline'
        write(25,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(25,'(a)')'\hline'
        write(25,'(a)')'U[eV] & Te (Iterations) &  Te (Iterations) &  Te (Iterations) & '//&
             ' Te (Iterations) &  Te (Iterations) &  Te (Iterations) \tabularnewline'
        write(25,'(a)')'\hline' 
        write(25,'(a)')'\hline'
     end if
     vU = dU * iU 
     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call pressure_to_temperature(vU,NaTrial*1000000.0,iError)
        Te_I(iN) = Te 
        iIter_I(iN) = iIterTe
        if(iError/=0)Te_I(iN)=-Te_I(iN)
     end do

     write(25,'(f6.0,6(a,f7.1,a,i7,a),a)') vU,&
          (' & ', Te_I(iN), ' (', iIter_I(iN),')', iN=0,nN ),'\tabularnewline'
     write(25,'(a)')'\hline'

  end do


  write(25,'(a)')'\end{tabular}'
  close(25)

  !_______________________________________________
  open(24,file='../doc/Table5.tex',status='replace')
  write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Te[eV] & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ '//&
       '& $U_{av} | C_v$ & $U_{av} | C_v$\tabularnewline'
  write(24,'(a)')'\hline' 
  write(24,'(a)')'\hline'

  call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)
  UsePreviousTe = .false.
  do iT  = 1,nT
     if (((iT-1)/25)*25==(iT-1).and.iT>25) then

        write(24,'(a)')'\end{tabular}', char(10)
        !--------------
        write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Te[eV] & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ '//&
             '& $U_{av} | C_v$ & $U_{av} | C_v$\tabularnewline'
        write(24,'(a)')'\hline' 
	write(24,'(a)')'\hline' 
     end if
     vTe = dTe * iT

     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
        Uav_I(iN)= internal_energy()
        Cv_I(iN) = heat_capacity()
        if(iError/=0)then
           Z_I(iN) = -  Z_I(iN)
           Z2Out_I(iN)= - Z2Out_I(iN)
           Uav_I(iN)= -Uav_I(iN)
           Cv_I(iN)= - Cv_I(iN)
           write(*,*)'Error=',iError,vTe,NaTrial*1000000.0
        end if
     end do
     write(24,'(f5.0,6(a,f8.1,a,f7.1),a)') vTe,&
          (' & ', Uav_I(iN), ' | ', Cv_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(a)')'\hline'
  end do
  write(24,'(a)')'\end{tabular}'

  close(24)

  open(24,file='../doc/Table6.tex',status='replace')
  write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(24,'(a)')'\hline'

  write(24,'(a)')'Te[eV] & $DH | Mad$ & $DH | Mad$ & $DH | Mad$ & $DH | Mad$ '//&
       '& $DH | Mad$ & $DH | Mad$\tabularnewline'

  write(24,'(a)')'\hline' 
  write(24,'(a)')'\hline'

  call set_element(54)
  UsePreviousTe = .true.
  do iT  = 1,nT
     if (((iT-1)/25)*25==(iT-1).and.iT>25) then

        write(24,'(a)')'\end{tabular}', char(10)
        !--------------
        write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Te[eV] & $DH | Mad$ & $DH | Mad$ & $DH | Mad$ & $DH | Mad$ '//&
             '& $DH | Mad$ & $DH | Mad$\tabularnewline'
        write(24,'(a)')'\hline' 
        write(24,'(a)')'\hline' 
     end if
     vTe = dTe * iT
     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
        Uav_I(iN)= eDebyeHuekel 
        Cv_I(iN) = eMadelungPerTe * 3.0 * Te * Z2
        if(iError/=0)then
           Z_I(iN) = -  Z_I(iN)
           Z2Out_I(iN)= - Z2Out_I(iN)
           Uav_I(iN)= -Uav_I(iN)
           Cv_I(iN)= - Cv_I(iN)
           write(*,*)'Error=',iError,vTe,NaTrial*1000000.0
        end if
     end do
     write(24,'(f5.0,6(a,f8.1,a,f7.1),a)') vTe,&
          (' & ', Uav_I(iN), ' | ', Cv_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(a)')'\hline'


  end do

  write(24,'(a)')'\end{tabular}'

  close(24)

  UseCoulombCorrection = .true.
  open(24,file='../doc/Table7.tex')
  write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|c|c|c|c|c|c|}'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Na[$1/cm^3$] & '//&
       '\multicolumn{2}{|c|}{$10^{18}$} & \multicolumn{2}{|c|}{$10^{19}$} & '//&
       '\multicolumn{2}{|c|}{$10^{20}$} & \multicolumn{2}{|c|}{$10^{21}$} & '//&
       '\multicolumn{2}{|c|}{$10^{22}$} & \multicolumn{2}{|c|}{$10^{23}$}\tabularnewline'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Te[eV] &'//&
       ' no & Mad & no & Mad & no & Mad & no & Mad & no & Mad & no & Mad\tabularnewline'
  write(24,'(a)')'\hline' 
  write(24,'(a)')'\hline'


  do iT  = 1,nT
     if (((iT-1)/25)*25==(iT-1).and.iT>25) then

        write(24,'(a)')'\end{tabular}', char(10)
        !--------------
        write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|c|c|c|c|c|c|}'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Na[$1/cm^3$] & '//&
             '\multicolumn{2}{|c|}{$10^{18}$} & \multicolumn{2}{|c|}{$10^{19}$} & '//&
             '\multicolumn{2}{|c|}{$10^{20}$} & \multicolumn{2}{|c|}{$10^{21}$} & '//&
             '\multicolumn{2}{|c|}{$10^{22}$} & \multicolumn{2}{|c|}{$10^{23}$}\tabularnewline'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Te[eV] &'//&
             ' no & Mad & no & Mad & no & Mad & no & Mad & no & Mad & no & Mad\tabularnewline'
        write(24,'(a)')'\hline' 
        write(24,'(a)')'\hline'
     end if
     vTe = dTe * iT

     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
        Z_I(iN) = z_averaged() 
        Z2Out_I(iN)= z2_averaged()/Z_I(iN)
        Uav_I(iN)=internal_energy()
        Cv_I(iN)=heat_capacity()
        if(iError/=0)then
           Z_I(iN) = -  Z_I(iN)
           Z2Out_I(iN)= - Z2Out_I(iN)
           Uav_I(iN)= -Uav_I(iN)
           Cv_I(iN)= - Cv_I(iN)
           write(*,*)'Error=',iError,vTe,NaTrial*1000000.0
        end if
     end do
     write(24,'(f5.0,6(a,f7.1,a,f7.1),a)') vTe,&
          (' & ', ZNoCoulomb_II(iT,iN),' & ', Z_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(a)')'\hline'


  end do

  write(24,'(a)')'\end{tabular}'

  close(24)
  !________


end program saha

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  write(*,*)StringError
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test


