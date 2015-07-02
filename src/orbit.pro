;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf

dpy = 670

orbitangle = findgen(dpy+1) * 2.0 * !pi / dpy

; Earth
a = 1.000110
b = 0.034221
c = 0.001280
d = 0.000719
e = 0.000077

; Mars
a =  1.52
b =  0.04
c =  0.15
d = 0.0
e = 0.0

SunOrbitEccentricity = a                     + $
                       b*cos(OrbitAngle)    + $
                       c*sin(OrbitAngle)    + $
                       d*cos(2.*OrbitAngle) + $
                       e*sin(2.*OrbitAngle)

; 1.56
; 1.67
; 1.47
; 1.38

print, SunOrbitEccentricity(0)
print, SunOrbitEccentricity(dpy/4)
print, SunOrbitEccentricity(dpy/2)
print, SunOrbitEccentricity(3*dpy/4)

x = SunOrbitEccentricity * cos(orbitangle)
y = SunOrbitEccentricity * sin(orbitangle)

;plot, orbitangle*dpy/(2*!pi), SunOrbitEccentricity, xstyle = 1, ystyle = 1

setdevice, 'orbit.ps', 'l', 5, 0.95

ppp = 1
space = 0.05
pos_space, ppp, space, sizes

get_position, ppp, space, sizes, 0, pos

plot, x,y, pos = pos, thick = 2.0, $
  xtitle = "(AU)", ytitle = "(AU)"

xyouts, 0.15, 0.01, "Equinox", charsize=0.8

xyouts, -0.15, 0.01, "Equinox", alignment = 1.0, charsize=0.8

xyouts, -0.01, 0.15, "Northern Summer", orient = 90, charsize=0.8
xyouts, -0.01, -0.15, "Southern Summer", orient = 90, $
  alignment = 1.0, charsize=0.8


xyouts, x(dpy/4), y(dpy/4)+0.05, "Mars", alignment = 0.5, charsize = 1.3

; Earth
a = 1.000110
b = 0.034221
c = 0.001280
d = 0.000719
e = 0.000077

dpy = 365

orbitangle = findgen(dpy+1) * 2.0 * !pi / dpy

SunOrbitEccentricity = a                     + $
                       b*cos(OrbitAngle)    + $
                       c*sin(OrbitAngle)    + $
                       d*cos(2.*OrbitAngle) + $
                       e*sin(2.*OrbitAngle)

x = SunOrbitEccentricity * cos(orbitangle)
y = SunOrbitEccentricity * sin(orbitangle)

oplot, x,y, thick = 2.0

xyouts, x(dpy/4), y(dpy/4)+0.05, "Earth", alignment = 0.5, charsize = 1.3

SunOrbitEccentricity = 0.05

x = SunOrbitEccentricity * cos(orbitangle)
y = SunOrbitEccentricity * sin(orbitangle)

oplot, x,y

oplot, [-2.0,2.0], [0.0,0.0], linestyle = 1
oplot,  [0.0,0.0], [-2.0,2.0],linestyle = 1

closedevice

end
