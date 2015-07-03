!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModCIMPIT

  use Mod_Glow, only: ALTSMAX
  implicit none


  real  :: ALPHA(ALTSMAX), BETA(ALTSMAX), GAMMA2(ALTSMAX), PSI(ALTSMAX),&
       DELZ(ALTSMAX), DEL2(ALTSMAX), DELA(ALTSMAX), DELP(ALTSMAX), &
       DELM(ALTSMAX), DELS(ALTSMAX), DEN(ALTSMAX), FAC

end module ModCIMPIT
