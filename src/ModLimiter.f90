!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModLimiterGitm
  implicit none
contains
  !============================================================================
  real function Limiter_minmod(dUp, dDown)

    real :: dUp, dDown

    Limiter_minmod = (sign(0.5,dUp) + sign(0.5,dDown))*min(abs(dUp),abs(dDown))

  end function Limiter_minmod
  !============================================================================
  real function Limiter_mc(dUp, dDown)

    use ModInputs, only:BetaLimiter

    real :: dUp, dDown

    if (dUp > 0.0) then
       if (dDown > 0.0) then
          Limiter_mc = min(BetaLimiter*dUp,BetaLimiter*dDown,(dUp+dDown)*0.5)
       else
          Limiter_mc = 0.0
       endif
    else
       if (dDown < 0.0) then
          Limiter_mc = max(BetaLimiter*dUp,BetaLimiter*dDown,(dUp+dDown)*0.5)
       else
          Limiter_mc = 0.0
       endif
    endif

  end function Limiter_mc

end module ModLimiterGitm
