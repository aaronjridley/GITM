
filelist = findfile('log*.dat')

thermo_readlog, filelist, data, time

setdevice, 'log.ps','p',5

stime = time(0)
c_r_to_a, itime, stime
c_a_to_s, itime, strtime

t = time-time(0)
cTime = 'Simulation Time'

if (max(t) gt 7200.0) then begin
    t = t/3600.0
    cTime = cTime+' (Hours)'
    nDays = max(t)/(24.0)
endif else begin
    cTime = cTime+' (Seconds)'
    nDays = 0
endelse

plot, t, data(10,*), linestyle = 2, ytitle = 'Temperature (K)', $
  xtitle = cTime, xstyle = 1, pos = [0.1,0.3,0.9,0.7]
oplot, t, data(9,*), linestyle = 2
oplot, t, data(11,*), thick = 3

xyouts, 0.9, 0.705, /norm, 'Start Time : '+strtime, alignment = 1

if (nDays gt 1) then begin
   for i=1, nDays do begin
      oplot, [i,i]*24.0, [0, max(t)*10.0], linestyle = 1
   endfor
endif

closedevice

end


