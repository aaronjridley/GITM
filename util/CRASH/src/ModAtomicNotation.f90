!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModAtomicNotation
  use CRASH_ModAtomicMass, ONLY: nMixMax

  implicit none

  !Chemistry and spectroscopy
  !\
  ! The Mendeleev table
  !/
  integer, parameter :: MaterialMin_=-1, MaterialMax_ = 79
  character(LEN=2),parameter:: NameElement_I(MaterialMin_:MaterialMax_)=(/&
     'Ay' ,'Pl',                                                                               & !-1:0
     'H_' ,                                                                               'He',& ! 1:2
     'Li','Be',                                                  'B_','C_','N_','O_','F_','Ne',& ! 3:10
     'Na','Mg',                                                  'Al','Si','P_','S_','Cl','Ar',& !11:18
     'K_','Ca','Sc','Ti','V_','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',& !19:36
     'Rb','Sr','Y_','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I_','Xe',& !37:54
     'Cs','Ba',                                                                                & !55:56
               'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Th','Dy','Ho','Er','Tm','Yb','Lu',     & !57:71
                    'Hf','Ta','W_','Re','Os','Ir','Pt',                                        & !72:78
     'Au'/)                                                                                      !79

  integer, parameter, dimension(1:nMixMax,MaterialMin_:0) :: nZMixStored_II = reshape( (/&
                                                  6, 8, 1, 0, 0, 0, & ! Ay (acrylic)
                                                  6, 1, 7, 8, 0, 0  & ! Pl (polyimide)
                                                  /),(/nMixMax,1-MaterialMin_/) )

  real,    parameter, dimension(1:nMixMax,MaterialMin_:0) :: cMixStored_II = reshape( (/& 
                                             5.0/15,  2.0/15, 8.0/15, 0.0   , 0.0, 0.0,& ! Ay (acrylic)
                                            22.0/39, 10.0/39, 2.0/39, 5.0/39, 0.0, 0.0 & ! Pl (polyimide)
                                            /),(/nMixMax,1-MaterialMin_/) )
       
  
  !Orbital quantum number: transform a numerical value to a symbol:
  character(LEN=1),parameter::TypeL_I(0:9) = (/'s','p','d','f','g','h','i','k','l','m'/)
contains
    integer function l_orbital(TypeL)
    character(LEN=1),intent(in)::TypeL
    !--------------------------------!
    l_orbital=-1
    select case(TypeL)
    case('s','S')
       l_orbital = 0
    case('p','P')
       l_orbital = 1
    case('d','D')
       l_orbital = 2
    case('f','F')
       l_orbital = 3
    case('g','G')
       l_orbital = 4
    case('h')
       l_orbital = 5
    case('i')
       l_orbital = 6
    case('k')
       l_orbital = 7
    case('l')
       l_orbital = 8
    case('m')
       l_orbital = 9
    case default
       call CON_stop('The spectroscopy symbol '//TypeL//' is not implemented')
    end select
  end function l_orbital
  !===========================
  integer function i_material(TypeMaterial)
    character(LEN=2),intent(in)::TypeMaterial
    integer:: iLoop
    !------------
    i_material = MaterialMin_-1
    do iLoop = MaterialMin_,MaterialMax_
       if(TypeMaterial==NameElement_I(iLoop))then
          i_material = iLoop
          return
       end if
    end do
  end function i_material
  !====================================================================

end module CRASH_ModAtomicNotation
