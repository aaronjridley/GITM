
filelist = findfile('3DALL*.bin')
if (n_elements(filelist) eq 1) then filelist = findfile('3DLST*.bin')

if (n_elements(filelist) eq 1) then begin
   print, 'Something is wrong here. Only one file found.'
   stop
endif

file = filelist(0)

gitm_read_header, file, StartTime, nVars, Vars, $
                  nLons, nLats, nAlts, version

file = filelist(n_elements(filelist)-1)

gitm_read_header, file, EndTime, nVars, Vars, $
                  nLons, nLats, nAlts, version

; In 3DLST, east velocity is  Vn(east)(m/s)
; In 3DALL, east is V!Dn!N(east)

sEast = ['V!Dn!N(east)','Vn(east)(m/s)']
iEast = -1
for iVar = 0, nVars-1 do begin
   for iE = 0,n_elements(sEast)-1 do begin
      if (iEast eq -1 and $
          strpos(Vars(iVar),sEast(iE)) gt -1) then iEast = iVar
   endfor
endfor

if (iEast eq -1) then begin
   print, 'Cant find east variable!'
   stop
endif else print, 'Eastward Velocity: ',Vars(iEast)

goce_dir = '/raid3/Data/GOCE'

read_goce, goce_dir, goce_data, starttime, endtime

if (goce_data.nPts lt 2) then begin
   print,"No GOCE data available for selected time!"
   c_r_to_a, iTime, StartTime
   c_a_to_s, iTime, sTime
   print,"Start Time : ", sTime
   c_r_to_a, iTime, EndTime
   c_a_to_s, iTime, eTime
   print,"End Time : ", eTime
endif

nPts = goce_data.nPts

VarsToGet = [3,iEast,iEast+1,iEast+2]
GoceTimes = goce_data.time
GoceLons  = goce_data.lon
GoceLats  = goce_data.lat
GoceAlts  = goce_data.alt/1000.0
; Let's adjust the Goce Altitude, since GITM assumes a
; spherical Earth.
re = 6378.137
rp = 6356.75
diff = re-rp
rgitm = 6372.0
r = rp+diff*cos(GoceLats*!dtor)
correction = r - rgitm
GoceAlts = GoceAlts + correction

GoceDensity = goce_data.density

GoceMag = sqrt(goce_data.east^2 + goce_data.north^2 + goce_data.vertical^2)
GoceEast = goce_data.east
GoceNorth = goce_data.north
GoceVertical = goce_data.vertical

UnitE = GoceEast/GoceMag
UnitN = GoceNorth/GoceMag
UnitV = GoceVertical/GoceMag

get_gitm_points, GoceTimes, GoceLons, GoceLats, GoceAlts, VarsToGet, GitmData

GitmDensity = reform(GitmData(*,0))
GitmEast     = reform(GitmData(*,1))
GitmNorth    = reform(GitmData(*,2))
GitmVertical = reform(GitmData(*,3))

GitmCross = GitmEast*UnitE + GitmNorth*UnitN + GitmVertical*UnitV

GitmCrossE = GitmCross * UnitE
GitmCrossN = GitmCross * UnitN
GitmCrossV = GitmCross * UnitV

LocalTime = (GoceLons/15.0 + (GoceTimes/3600.0)) mod 24.0

dl           = fltarr(nPts)
dl(1:nPts-2) = (GoceLats(2:nPts-1) - GoceLats(0:nPts-3))/2
dl(0)        = dl(1)
dl(nPts-1)   = dl(nPts-2)

AltMean = mean(GoceAlts)

Re = 6378.0D
r = (Re + AltMean) * 1000.0D
G = 6.67259d-11
Me = 5.9722d24                  ; reference 1
mu = G * Me
v0 = sqrt(mu/r)
period = 2*!pi*r / v0 / 60.0

dLat = 1.0
nX = fix((EndTime-Starttime)/(period*60.0))+2
nY = fix(180.0/dlat) + 1

GitmDensityMapa  = fltarr(nX, nY)
GitmEastMapa     = fltarr(nX, nY)
GitmNorthMapa    = fltarr(nX, nY)

GoceDensityMapa  = fltarr(nX, nY)
GoceEastMapa     = fltarr(nX, nY)
GoceNorthMapa    = fltarr(nX, nY)

GitmDensityMapd  = fltarr(nX, nY)
GitmEastMapd     = fltarr(nX, nY)
GitmNorthMapd    = fltarr(nX, nY)

GoceDensityMapd  = fltarr(nX, nY)
GoceEastMapd     = fltarr(nX, nY)
GoceNorthMapd    = fltarr(nX, nY)

x = round((GoceTimes-StartTime)/(period*60.0))
y = round((GoceLats+90.0)/dLat)

maptime = fltarr(nX, nY)
maplat  = fltarr(nX, nY)

for i=0,nX-1 do maplat(i,*) = findgen(nY)*dLat-90.0 + dLat/2.0

l = where(dl gt 0 and abs(GoceLats) lt 45.0, c)
if (c gt 0) then lta = mean(localtime(l))

l = where(dl lt 0 and abs(GoceLats) lt 45.0, c)
if (c gt 0) then ltd = mean(localtime(l))

for i=0L,nPts-2 do begin
   if (dl(i) gt 0) then begin

      GitmDensityMapa(x(i),y(i))  = GitmDensity(i)
      GitmEastMapa(x(i),y(i))     = GitmCrossE(i)
      GitmNorthMapa(x(i),y(i))    = GitmCrossN(i)

      GoceDensityMapa(x(i),y(i))  = GoceDensity(i)
      GoceEastMapa(x(i),y(i))     = GoceEast(i)
      GoceNorthMapa(x(i),y(i))    = GoceNorth(i)

   endif else begin

      GitmDensityMapd(x(i),y(i)) = GitmDensity(i)
      GitmEastMapd(x(i),y(i))    = GitmCrossE(i)
      GitmNorthMapd(x(i),y(i))   = GitmCrossN(i)

      GoceDensityMapd(x(i),y(i)) = GoceDensity(i)
      GoceEastMapd(x(i),y(i))    = GoceEast(i)
      GoceNorthMapd(x(i),y(i))   = GoceNorth(i)

   endelse
endfor

for i=0,nY-1 do maptime(*,i) = dindgen(nX)*period*60.0

c_r_to_a, itime, StartTime
c_a_to_y_m_d, itime, ymd

psfile = 'gitm_east_ascent_'+ymd+'.ps'
title = 'Eastward Wind'
satellite = 'GOCE'
thermo_compare_sat_plot, maplat, maptime, GitmEastMapa, GoceEastMapa, $
                         lta, StartTime, EndTime, $
                         psfile, title, satellite

psfile = 'gitm_north_ascent_'+ymd+'.ps'
title = 'Northward Wind'
thermo_compare_sat_plot, maplat, maptime, GitmNorthMapa, GoceNorthMapa, $
                         lta, StartTime, EndTime, $
                         psfile, title, satellite

psfile = 'gitm_density_ascent_'+ymd+'.ps'
title = 'Density'
thermo_compare_sat_plot, maplat, maptime, GitmDensityMapa, GoceDensityMapa, $
                         lta, StartTime, EndTime, $
                         psfile, title, satellite

psfile = 'gitm_east_descent_'+ymd+'.ps'
title = 'Eastward Wind'
thermo_compare_sat_plot, maplat, maptime, GitmEastMapd, GoceEastMapd, $
                         ltd, StartTime, EndTime, $
                         psfile, title, satellite

psfile = 'gitm_north_descent_'+ymd+'.ps'
title = 'Northward Wind'
thermo_compare_sat_plot, maplat, maptime, GitmNorthMapd, GoceNorthMapd, $
                         ltd, StartTime, EndTime, $
                         psfile, title, satellite

psfile = 'gitm_density_descent_'+ymd+'.ps'
title = 'Density'
thermo_compare_sat_plot, maplat, maptime, GitmDensityMapd, GoceDensityMapd, $
                         ltd, StartTime, EndTime, $
                         psfile, title, satellite

nTimesA = (max(GoceTimes)-min(GoceTimes)) / (period*60.0)

GoceAve = fltarr(nTimesA)
GitmAve = fltarr(nTimesA)
TimesAve = dindgen(nTimesA)*period*60.0 + min(GoceTimes) + period*60.0/2.0

for iTime = 0L, nTimesA-1 do begin
   loc = where(abs(GoceTimes-TimesAve(iTime)) lt period*60.0/2, count)
   if (count gt 0) then begin
      GitmAve(iTime) = mean(GitmDensity(loc))/1.0e-11
      GoceAve(iTime) = mean(GoceDensity(loc))/1.0e-11
   endif
endfor

ytitle = 'Mass Density (10!E-11!N kg/m!E3!N)'


psfile = 'gitm_goce_density_ave.ps'
setdevice, psfile, 'p', 5

ppp = 3
space = 0.04
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

get_position, ppp, space, sizes, 2, pos2, /rect
pos2(0) = pos2(0) + 0.05
pos2(2) = pos2(2) - 0.05

stime = min(TimesAve)
etime = max(TimesAve)
time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn

maxi = max([GoceAve,GitmAve])
yrange = [0.0, maxi]
plot, TimesAve-stime, GoceAve, yrange = yrange*1.2, pos = pos, $
      xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = ytitle,   $
      thick = 3, /noerase, min_val = 0.01
oplot, TimesAve-stime, gitmAve, linestyle = 2, thick = 3, min_val = 0.01
    
t1 = (etr-btr)*0.05
t2 = (etr-btr)*0.10
t3 = (etr-btr)*0.11
    
oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3, linestyle = 2
xyouts, t3, max(yrange)*1.05, 'GITM'
    
t1 = (etr-btr)*0.55
t2 = (etr-btr)*0.60
t3 = (etr-btr)*0.61
    
oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3
xyouts, t3, max(yrange)*1.05, Satellite

d = GoceAve-gitmAve 
yr = [-maxi,maxi]

plot, TimesAve-stime, d, $
      yrange = yr, pos = pos2, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = ytitle,   $
      thick = 3, /noerase, ystyle = 1
    
l = where(GoceAve gt 0.0 and gitmAve gt 0.0)
rmse = sqrt(mean((GoceAve(l)-gitmAve(l))^2))
rmsd = sqrt(mean((GoceAve(l))^2))
nrms = rmse/rmsd * 100.0
    
pdif = mean((GoceAve(l)-gitmAve(l))/GoceAve(l)) * 100.0
    
srms = +' (nRMS: '+string(nrms,format = '(f5.1)')+'%, '
srms = srms+string(pdif,format = '(f5.1)')+'% Difference)'
    
xyouts, pos2(2)+0.01, (pos2(1)+pos2(3))/2.0, srms, $
        alignment = 0.5, orient=270, /norm 

xyouts, pos2(2)+0.035, (pos2(1)+pos2(3))/2.0, 'GOCE - GITM', $
        alignment = 0.5, orient=270, /norm 

oplot, [btr,etr], [0.0,0.0], linestyle = 1
    
spawn, 'pwd', dir
dir = dir(0)
xyouts, 0.0, -0.05, dir, /norm, charsize = 0.8
    
closedevice

;---------------------------------------------------------------------------
; raw!
;---------------------------------------------------------------------------

psfile = 'gitm_goce_density_raw.ps'
setdevice, psfile, 'p', 5

ppp = 3
space = 0.04
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

get_position, ppp, space, sizes, 2, pos2, /rect
pos2(0) = pos2(0) + 0.05
pos2(2) = pos2(2) - 0.05

stime = min(GoceTimes)
etime = max(GoceTimes)
time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn

maxi = max([GoceDensity,GitmDensity])
yrange = [0.0, maxi]
plot, Gocetimes-stime, GoceDensity, yrange = yrange*1.2, pos = pos, $
      xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = ytitle,   $
      thick = 3, /noerase
oplot, GoceTimes-stime, gitmDensity, linestyle = 2, thick = 3
    
t1 = (etr-btr)*0.05
t2 = (etr-btr)*0.10
t3 = (etr-btr)*0.11
    
oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3, linestyle = 2
xyouts, t3, max(yrange)*1.05, 'GITM'
    
t1 = (etr-btr)*0.55
t2 = (etr-btr)*0.60
t3 = (etr-btr)*0.61
    
oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3
xyouts, t3, max(yrange)*1.05, Satellite

d = GoceDensity-gitmDensity
yr = [-maxi,maxi]

plot, GoceTimes-stime, d, $
      yrange = yr, pos = pos2, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = ytitle,   $
      thick = 3, /noerase, ystyle = 1
    
rmse = sqrt(mean((GoceDensity-gitmDensity)^2))
rmsd = sqrt(mean((GoceDensity)^2))
nrms = rmse/rmsd * 100.0
    
pdif = mean((GoceDensity-gitmDensity)/GoceDensity) * 100.0
    
srms = +' (nRMS: '+string(nrms,format = '(f5.1)')+'%, '
srms = srms+string(pdif,format = '(f5.1)')+'% Difference)'
    
xyouts, pos2(2)+0.01, (pos2(1)+pos2(3))/2.0, srms, $
        alignment = 0.5, orient=270, /norm 

xyouts, pos2(2)+0.035, (pos2(1)+pos2(3))/2.0, 'GOCE - GITM', $
        alignment = 0.5, orient=270, /norm 

oplot, [btr,etr], [0.0,0.0], linestyle = 1
    
spawn, 'pwd', dir
dir = dir(0)
xyouts, 0.0, -0.05, dir, /norm, charsize = 0.8
    
closedevice


WriteFile = 1
if (WriteFile) then begin
   openw,1,'gitm_goce_'+ymd+'.dat'
   printf,1,'year month day hour minute second longitude latitude altitude gitm_rho gitm_east gitm_north gitm_vertical gitm_east_projected gitm_north_projected goce_rho goce_east_projected goce_north_projected goce_vertical_projected'
   for i=0L,nPts-2 do begin
      c_r_to_a, iTime, GoceTimes(i)
      printf,1,iTime, GoceLons(i), GoceLats(i), GoceAlts(i), $
             GitmDensity(i),GitmEast(i),GitmNorth(i),GitmVertical(i), $
             GitmCrossE(i),GitmCrossN(i), GoceDensity(i), $
             GoceEast(i),GoceNorth(i),GoceVertical(i), $
             format = '(6i5,3f7.1,e10.3,5f7.1,e10.3,3f7.1)' 
   endfor
   close,1
endif

end
