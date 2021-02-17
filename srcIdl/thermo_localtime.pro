
filelist_new = [findfile('1DALL_t120[12]*bin'),findfile('1DALL_t120[34]*bin')];,findfile('1DALL_t120[56]*bin')]
nfiles_new = n_elements(filelist_new)
getnewdata = 1

if n_elements(nfiles) gt 0 then begin
    if (nfiles_new eq nfiles) then default = 'n' else default='y'
    if (strpos(old_start,start) eq -1) then default = 'n'
    GetNewData = mklower(strmid(ask('whether to reread data',default),0,1))
    if (GetNewData eq 'n') then GetNewData = 0 else GetNewData = 1
endif

if (GetNewData) then begin

;    thermo_readsat, filelist_new, data, time, nTimes, Vars, nAlts, nSats, Files
;    nFiles = n_elements(filelist_new)

    gitm_read_bin, filelist_new, data, time, nVars, Vars, version
    nTimes = n_elements(time)
    nAlts = n_elements(data(0,0,0,0,*))
    
    ; Need to reform the data
    newdata = dblarr(1,nTimes,nVars,nAlts)
    newdata(0,*,*,*) = data(*,*,0,0,*)
    data = newdata

    nSats = 1

endif

iAlt = 5

iVar = 8

ntimes = n_elements(time)
dt = round(mean(time(1:ntimes-1)-time(0:ntimes-2)))

x = fix((time-time(0))/86400.0)
y = ((time-time(0)) mod 86400.0)/dt

ny = max(y)+1
nx = max(x)+1

v = fltarr(nx,ny)
x2d = fltarr(nx,ny)
y2d = fltarr(nx,ny)
x2d(x,y) = x
x2d(nx-1,*) = x2d(nx-1,0)
y2d(x,y) = y*dt/3600.0
y2d(nx-1,*) = y2d(nx-2,*)

setdevice, 'no_low.ps','p',4
makect,'all'

v(x,y) = newdata(0,*,iVar,iAlt)
v(nx-1,*)=v(nx-2,*)

contour, v, x2d, y2d, xtitle = 'Days', ytitle = 'hours', $
         /fill, nlevels = 31, title = vars(iVar),$
         pos =[0.1,0.3,1.0,0.8], xstyle = 1, ystyle = 1
closedevice

setdevice, 'o_low.ps','p',4
makect,'all'

iVar = 4
v(x,y) = newdata(0,*,iVar,iAlt)
v(nx-1,*)=v(nx-2,*)

contour, v, x2d, y2d, xtitle = 'Days', ytitle = 'hours', $
         /fill, nlevels = 31, title = vars(iVar),$
         pos =[0.1,0.3,1.0,0.8], xstyle = 1, ystyle = 1
closedevice

setdevice, 'temp.ps','p',4
makect,'all'

iVar = 15
iAlt = 30
v(x,y) = newdata(0,*,iVar,iAlt)
;v(nx-1,*)=v(nx-2,*)

contour, v, x2d, y2d, xtitle = 'Days', ytitle = 'hours', $
         /fill, nlevels = 31, title = vars(iVar),$
         pos =[0.1,0.3,1.0,0.8], xstyle = 1, ystyle = 1
closedevice

last = reform(v(nx-1,*))
mini = min(last)
maxi = max(last)
l = where(last eq mini)

print, mini,maxi,mini/maxi,y2d(nx-1,l(0)), mean(v(nx-1,*))

iVar = 34
iAlt = 30
v(x,y) = newdata(0,*,iVar,iAlt)
;v(nx-1,*)=v(nx-2,*)

contour, v, x2d, y2d, xtitle = 'Days', ytitle = 'hours', $
         /fill, nlevels = 31, title = vars(iVar),$
         pos =[0.1,0.3,1.0,0.8], xstyle = 1, ystyle = 1
closedevice

last = reform(v(nx-1,*))
mini = min(last)
maxi = max(last)
l = where(last eq mini)

print, mini,maxi,mini/maxi,y2d(nx-1,l(0)), mean(v(nx-1,*))

end



