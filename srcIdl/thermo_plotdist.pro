
filename = '2DUSR_t000320_000000.bin'

read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                        vars, data, rb, cb, bl_cnt, iTime, Version

alt = reform(data(2,*,*,*)) / 1000.0
lat = reform(data(1,*,*,*)) / !dtor
lon = reform(data(0,*,*,*)) / !dtor

nEnergies = nVars-10
Energies = fltarr(nEnergies)
for i=0,nEnergies-1 do begin
   Energies(i) = float(strmid(vars(i+10),5,9))
endfor

setdevice, 'aurora2.ps','p',5

plot_oi, energies,data(10:59,20,12), psym =-4, pos = [0.15,0.3,0.9,0.7], $
         xtitle = 'Energy (eV)', ytitle = 'Differential Flux (#/cm2/s/eV)',$
         thick = 3, charsize = 1.2

oplot, data(4,20,12)*[1000.0,1000.0], [0.0,1.0e7], thick = 3
oplot, data(6,20,12)*[1000.0,1000.0], [0.0,1.0e7], linestyle = 1, thick = 3
oplot, data(8,20,12)*[1000.0,1000.0], [0.0,1.0e7], linestyle = 2, thick = 3

closedevice

end

