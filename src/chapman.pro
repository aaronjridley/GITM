;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf

nLats = 36
nLons = 72
nAlts = 20

r = fltarr(72,36,20)
t = fltarr(72,36,20)
p = fltarr(72,36,20)

for i=0,nLats-1 do t(*,i,*)  = -90.0 + (float(i)+0.5)/nLats * 180.0
for i=0,nLons-1 do p(i,*,*)  = 180.0 + (float(i)+0.5)/nLons * 360.0
for i=0,nAlts-1 do r(*,*,i)  = 6372.0 + 95.0 + (float(i)+0.5) * 20.0

t = t * !pi / 180.0
p = (p mod 360.0) * !pi / 180.0
r = r * 1000.0

h = 5.0 * 1000.0

x = r * cos(t) * cos(p)
y = r * cos(t) * sin(p)
z = r * sin(t)

yz = sqrt(y^2 + z^2)

sza = acos(x / r)

xp = r / h

y = (0.5 * !pi * xp)^0.5 * abs(x/r)

ch = fltarr(72,36,20)

erfcy = fltarr(72,36,20)

a = 1.06069630
b = 0.55643831
c = 1.06198960
d = 1.72456090
f = 0.56498823
g = 0.06651874

loc = where(y lt 8)
erfcy(loc) = (a+b*y(loc)) / (c+d*y(loc)+y(loc)^2)
loc = where(y gt 8)
erfcy(loc) = f / (g + y(loc))

loc = where(sza lt !pi/2)

ch(loc) = h * (!pi * 0.5 * xp(loc))^0.5 * erfcy(loc)

rg = r * cos(sza-!pi/2)
loc = where((r - rg)/h lt 30.0)

eg = fltarr(72,36,20) + exp(30.0)
eg(loc) = exp((r(loc) - rg(loc))/h)

xg = rg/h

loc = where(sza gt !pi/2)

ch(loc) = (!pi * 0.5 * xg(loc))^0.5 * h * (2 * eg(loc) - erfcy(loc))

loc = where((r - rg)/h lt 10.0 and sza gt !pi/2)
maxi = max(ch(loc))

print, maxi

loc = where((r - rg)/h ge 10.0 and sza gt !pi/2)
ch(loc) = maxi

makect, 'mid'
contour, alog10(ch(*,*,10)+0.001), nlevels = 30, /cell_fill
contour, alog10(ch(*,*,10)+0.001), /follow, nlevels = 10, /noerase



end

