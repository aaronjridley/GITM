;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
filename='potentialtest.out'
.r getpict
cut=grid(0,*,*)
func='{pot}-(x-2.5^3/x^2)/(1+2*2.5^3)*5*sin(y)*cos(z)
plotmode='contbar
.r plotfunc
print,'L1 error for slice =',total(abs(f))/n_elements(f)
r     = x(*,*,*,0)
theta = x(*,*,*,1)
phi   = x(*,*,*,2)
print,'L1 error for volume=', $
  total(abs(w-(r-2.5^3/r^2)/(1+2*2.5^3)*5*sin(theta)*cos(phi)))/n_elements(w)
