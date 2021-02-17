
IsNorth = fix(ask('North (1) or South (1)','1'))
IncludeCusp = 0

IncludeIons = 1

max_efield_var = 0.0           ; mV/m
max_pot_var = max_efield_var / 1000.0 * 111.0 * 1000.0 /1000.0 

variability = max_pot_var
cc = 0.0

bz = [-2.0, -2.0, -2.0, -2.0]
by = [-5.0, -5.0,  5.0,  5.0]

if (IsNorth) then sHem = 'n' else sHem = 's'

filename = 'amie_test_'+sHem+'.bin'

iTime = [2000,3,21,0,0,0]
c_a_to_r, iTime, firsttime

time = firsttime + [0.0, 12.0, 12.1667, 24.0] * 3600.0D


nTimes = n_elements(time)

imf = fltarr(nTimes, 4)
ae = fltarr(nTimes, 4)
dst = fltarr(nTimes, 2)
hpi = fltarr(nTimes, 2)
cpcp = fltarr(nTimes)

by = by * (IsNorth*2 - 1)

; This sets the resolution:

dlat = 0.5
nLats = 45.0/dlat + 1

dmlt = 0.25
nMlts = 24.0/dmlt + 1

lats = 90.0-findgen(nLats)*dlat
mlts = findgen(nMlts)*dmlt

; This sets where the potential and aurora maximizes (potential has
; max/min at the ocflb, while the aurora is equatorward of this.)

ocflbBase = 70.0 - 3.0*cos(mlts*!pi/12.0)

nVars = 5

data = fltarr(nTimes, nVars, nMlts, nLats)

potential = fltarr(nMlts, nLats)
eflux     = fltarr(nMlts, nLats)
avee      = fltarr(nMlts, nLats)
ion_eflux     = fltarr(nMlts, nLats)
ion_avee      = fltarr(nMlts, nLats)
; area is used to calculate the hemispheric power
area      = fltarr(nMlts, nLats)
for iLat = 0,nLats-1 do $
   area(*,iLat) = 111.0*dlat*(111.0*360.0/nMlts) * cos(lats(iLat)*!dtor)*1000.0*1000.0


; Cusp, LLBL, and Mantle (HALF widths!!!)
cusp_ave = 0.1
cusp_eflux = 2.0
cusp_center_mlt = 12.0 + by/10.0
cusp_width_mlt = 1.0 + abs(by)/20.0
cusp_width_lat = 0.5

llbl_ave = 0.5
llbl_eflux = 1.0
llbl_center_mlt = 12.0 - by/5.0
llbl_width_mlt = 3.0 + abs(by)/10.0
llbl_width_lat = 1.0

mantle_ave = cusp_ave
mantle_eflux = cusp_eflux/10.0
mantle_center_mlt = llbl_center_mlt
mantle_width_mlt = llbl_width_mlt
mantle_width_lat = 3.0*cusp_width_lat

for i=0,nTimes-1 do begin

   ocflb = ocflbBase

   by_sharpen = (25.0-abs(by(i)))/25.0
   fac_pt_base = 7.0 * by_sharpen^0.5
   fac_ef_base = 1.0
   fac_ae_base = 3.0

   ; This is the amplitude of the potential in the harang
   ; discontinuity.  I would set this to zero, or you could play with
   ; it for a strong jet....

   amp_h = 10.0
   amp_h = amp_h(0)

   ; This is the amplitude, multiplied by Bz for the strength of the CPCP
   amp    = 10.0*1000.0
   ; Amplitude of the E-flux
   amp_ef = 2.0
   ; Amplitude of the average energy
   amp_ae = 2.5

   amp_h = amp_h * 1000.0

   m = (90.0-mean(ocflb))/1.5

   for iMlt = 0,nMlts-1 do begin

      fac = fltarr(nLats) + fac_pt_base/1.5
      fac_ef = fltarr(nLats) + fac_ef_base * (0.5+1.5*(cos(mlts(iMlt)*!pi/12.0)+1.0))/3.5
      ion_ef = fltarr(nLats) + 0.33 * fac_ef_base * $
               exp((cos((mlts(iMlt)-21.0)*15.0*!dtor)+1.0)/0.5)/exp(4)

      fac_ae = fltarr(nLats) + fac_ae_base
      ion_ae = fltarr(nLats) + fac_ae_base

      d = lats-ocflb(iMlt)

      ; build the potential pattern

      l = where(d ge 0)
      fac(l) = fac_pt_base
      line = exp(-abs(lats-ocflb(iMlt))/fac)

      a = amp
      if (by(i) gt 0 and mlts(iMlt) gt 12.0) then a = a * by_sharpen
      if (by(i) lt 0 and mlts(iMlt) lt 12.0) then a = a * by_sharpen

      potential(iMlt,*) = -bz(i) * a * line * sin(mlts(iMlt)*!pi/12.0)

      t = mlts(iMlt)*!pi/12.0 - !pi/2
      x = (90.0-lats) * cos(t)
      y = (90.0-lats) * sin(t)

      t0 = (12.0 + by(i)/4.0)*!pi/12.0 - !pi/2
      l0 = 85.0
      x0 = (90.0-l0) * cos(t0)
      y0 = (90.0-l0) * sin(t0)

      r = sqrt((x-x0)^2 + (y-y0)^2)
      
      ; Some By
      line = -amp/2 * by(i) * exp(-r/(m/1.5))
      potential(iMlt,*) = potential(iMlt,*) + line

      ; Harang Discontinuity - centered at midnight
      line = -amp_h * exp(-2*(lats-ocflb(iMlt))^2/fac^2) * $
             ((cos((mlts(iMlt)-1)*!pi/12.0)+1.0)/2.0)^3.0
      l = where(line gt 0,c)
      if (c gt 0) then line(l) = 0.0
      potential(iMlt,*) = potential(iMlt,*) + line

      ; build the energy flux pattern
      l = where(d lt 0)
      fac_ef(l) = fac_ef(l)*2
      line = exp(-abs(lats-ocflb(iMlt))^4/fac_ef^4)
      eflux(iMlt,*) = amp_ef * line * (cos(mlts(iMlt)*!pi/12.0)+2.0)/3.0

      ; ion eflux
      line = exp(-abs(lats-(ocflb(iMlt)-2.5*mean(fac_ef)))^4/fac_ef^4)
      ion_eflux(iMlt,*) = ion_ef * line

      ; build the averaged energy pattern
      line = exp(-(lats-ocflb(iMlt))^2/fac_ae^2)
      avee(iMlt,*) = amp_ae * line * (cos(mlts(iMlt)*!pi/12.0)+5.0)/6

      ; ion average energy:
      line = exp(-(lats-(ocflb(iMlt)-2.5*mean(fac_ef)))^2/fac_ae^2)
      ion_avee(iMlt,*) = 7.0 * amp_ae * line * (cos((mlts(iMlt)-21.0)*!pi/12.0)+5.0)/6

;      eflux(iMlt,*) = eflux(iMlt,*) + $
;                      amp_ef*10 * line * exp(-(mlts(iMlt)-12-by(i)*0.05)^2/2.0)

;      eflux(iMlt,*) = eflux(iMlt,*) + $
;                      cc * vari_line * eflux(iMlt,*) + $
;                      sign(cc) * (1.0 - abs(cc)) * vari_line2 * eflux(iMlt,*)

      ; Add LLBL:
      llbl_lat_shift = -llbl_width_lat/2.5 * by_sharpen
      p_mlt = ((llbl_center_mlt(i) - mlts(iMlt))^4)/llbl_width_mlt(i)^4
      p_lat = ((lats - (ocflb(iMlt)+llbl_lat_shift))^4)/llbl_width_lat^4
      llbl_eflux_line = llbl_eflux * exp(-p_mlt) * exp(-p_lat)
      eflux(iMlt,*) = eflux(iMlt,*) + llbl_eflux_line
      
      ; Add Mantle:
      p_mlt = ((mantle_center_mlt(i) - mlts(iMlt))^4)/mantle_width_mlt(i)^4
      p_lat = ((lats - (ocflb(iMlt) - llbl_lat_shift + mantle_width_lat/2))^4)/mantle_width_lat^4
      mantle_eflux_line = mantle_eflux * exp(-p_mlt) * exp(-p_lat)
      eflux(iMlt,*) = eflux(iMlt,*) + mantle_eflux_line
      
      ; Add Cusp:
      p_mlt = ((cusp_center_mlt(i) - mlts(iMlt))^4)/cusp_width_mlt(i)^4
      p_lat = ((lats - (ocflb(iMlt)+cusp_width_lat/2))^4)/cusp_width_lat^4
      cusp_eflux_line = cusp_eflux * exp(-p_mlt) * exp(-p_lat)
      eflux(iMlt,*) = eflux(iMlt,*) + cusp_eflux_line
      
      ; Add Mantle to avee:
      l = where(mantle_eflux_line gt 0.5*mantle_eflux,c)
      if (c gt 0) then avee(iMlt,l) = mantle_ave

      ; Add LLBL to avee:
      l = where(llbl_eflux_line gt 0.5*llbl_eflux,c)
      if (c gt 0) then avee(iMlt,l) = llbl_ave

      ; Add cusp to avee:
      l = where(cusp_eflux_line gt 0.5*cusp_eflux,c)
      if (c gt 0) then avee(iMlt,l) = cusp_ave

      
   endfor
;stop
   
   ; set base values for aurora

   l = where(avee lt 0.05)
   avee(l) = 0.05

   l = where(eflux lt 0.02)
   eflux(l) = 0.002

   l = where(ion_avee lt 0.05)
   ion_avee(l) = 0.05

   l = where(ion_eflux lt 0.02)
   ion_eflux(l) = 0.002

   ; fill the data array

   data(i,0,*,*) = potential
   data(i,1,*,*) = eflux
   data(i,2,*,*) = avee
   data(i,3,*,*) = ion_eflux
   data(i,4,*,*) = ion_avee

   ; calculate hemispheric power

   hpi(i,0) = total(area*eflux*0.001)/1.0e9
   hpi(i,1) = total(area*ion_eflux*0.001)/1.0e9

   ; calculate cpcp

   cpcp(i) = max(potential)-min(potential)

   print, 'hp (e/i): ', i, hpi(i,0), hpi(i,1), $
          amp_h/1000.0, by(i), cpcp(i)/1000.0

endfor

imf(*,0) = 400.0
imf(*,1) = 0.0
imf(*,2) = by
imf(*,3) = bz

Version = 1.0

vars = ['Potential (kV)', $
        'Total Energy Flux (ergs/cm2/s)', $
        'Mean Energy (ergs)', $
        'Ion Energy Flux (ergs/cm2/s)', $
        'Ion Mean Energy (ergs)']

; write amie file

amie_write_binary, filename, Vars, lats, mlts, time, data, $
                   imf = imf, ae = ae, dst = dst, hpi = hpi, cpcp = cpcp, $
                   Version = Version

lats2d = fltarr(nmlts,nlats)
mlts2d = fltarr(nmlts,nlats)

for i=0,nLats-1 do mlts2d(*,i) = mlts
for i=0,nmlts-1 do lats2d(i,*) = lats

x = (90.0-lats2d)*cos(mlts2d/12.0*!pi-!pi/2)
y = (90.0-lats2d)*sin(mlts2d/12.0*!pi-!pi/2)

end

