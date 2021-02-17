
ChampDir = '/raid3/Data/CHAMP/'
GraceDir = '/raid3/Data/Grace/'

reread = 1
if (n_elements(GitmData) gt 0) then begin
    answer = ask('whether to re-read gitm data','n')
    if (strpos(mklower(answer),'n') gt -1) then reread = 0
endif

if (reread) then begin

   if (n_elements(SatType) eq 0) then SatType = 'champ'
   SatType = ask('satellite type (grace or champ)',SatType)

   if (strpos(SatType,'cha') eq 0) then begin
      Satellite = 'CHAMP'
      IsChamp = 1
   endif else begin
      Satellite = 'Grace-B'
      IsChamp = 0
   endelse

   print, 'Satellite : ',Satellite

   filelist = file_search(SatType+'_*.bin')
   gitm_read_bin, filelist, GitmData, GitmTime, nVars, Vars, version, /skip
   nTimes = n_elements(GitmTime)
   GitmAlts = reform(GitmData(0,2,0,0,*))/1000.0
   nAlts = n_elements(GitmAlts)

endif

c_r_to_a, itime_start, min(gitmtime)
c_r_to_a, itime_end, max(gitmtime)

c_a_to_y_m_d, itime_start, stime_start
c_a_to_y_m_d, itime_end, stime_end

print, 'Start Time : ', stime_start
print, 'End   Time : ', stime_end

start_time = ask('start time of plotting (yyyy-mm-dd)',stime_start)
start_time = strmid(start_time,0,4)+$
             strmid(start_time,5,2)+$
             strmid(start_time,8,2)
end_time   = ask('end time of plotting',stime_end)
end_time = strmid(end_time,0,4)+$
           strmid(end_time,5,2)+$
           strmid(end_time,8,2)

c_ymd_to_a, itime_start, start_time
c_a_to_r, itime_start, rtime_start

c_ymd_to_a, itime_end, end_time
c_a_to_r, itime_end, rtime_end

c_r_to_a, itime, rtime_start
itime(3:5) = 0
c_a_to_r,  itime, stime

c_r_to_a, ietime, rtime_end ; max(GitmTime)
ietime(3:5) = [23,59,59]
c_a_to_r, ietime, etime

nDays = round((etime-stime)/3600.0/24.0)
doy0 = jday(itime(0), itime(1), itime(2))
iYear = itime(0)
sYear = tostr(iYear)

l = where(gitmtime ge stime and gitmtime le etime, nTimes)
if (nTimes gt 0) then begin
   GitmTime = GitmTime(l)
   GitmData = GitmData(l,*,*,*,*)
endif else stop

if (IsChamp) then begin
   if (n_elements(type) eq 0) then type = '0' else type = tostr(type)
   type = fix(ask('data type to compare to (0-rho, 1-electron, 2-winds)',type))
endif else type = 0

Lats = reform(GitmData(*, 1,0,0,0))/!dtor
LocalTime = (reform(GitmData(*, 0,0,0,0))/!dtor/15.0 + $
            (GitmTime/3600.0)) mod 24.0

; mass density
if (type eq 0) then GitmRho = reform(GitmData(*, 3,0,0,*))
; electron density
if (type eq 1) then GitmRho = reform(GitmData(*,34,0,0,*))
; eastward wind
if (type eq 2) then begin
   east  = reform(GitmData(*,16,0,0,*))
   north = reform(GitmData(*,17,0,0,*))
   Lons = reform(GitmData(*, 0,0,0,0))/!dtor
   GitmRho = reform(GitmData(*,16,0,0,*))
endif
; read champ data

for cday = 0, ndays - 1 do begin

   doy = doy0 + cday

   if (type eq 0) then begin
      if (IsChamp) then $
         filehead = champdir+sYear+'/Density_3deg_' $
      else filehead = gracedir+sYear+'/ascii/Density_graceB_3deg_'
      champfile = filehead+$
                  tostr((iYear mod 100),2)+'_'+ $
                  chopr('00'+tostr(doy),3)+'.ascii'

      print, 'Reading '+champfile
      read_champ_density, champfile, ChampTime, ChampData, ChampVars, note

      if (cDay eq 0) then begin
         ChampAlt = reform(champdata( 6,*))
         ChampRho = reform(champdata(11,*))
         ChampErr = reform(champdata(15,*))
         ChampLat = reform(champdata(3,*))
         ChampTim = ChampTime
      endif else begin
         ChampAlt = [ChampAlt,reform(champdata( 6,*))]
         ChampRho = [ChampRho,reform(champdata(11,*))]
         ChampErr = [ChampErr,reform(champdata(15,*))]
         ChampLat = [ChampLat,reform(champdata(3,*))]
         ChampTim = [ChampTim,ChampTime]
      endelse

   endif else if (type eq 1) then begin

      itime = [iyear,1,doy,0,0,0]
      c_a_to_r, itime, rtime
      c_r_to_a, itime, rtime
      iMon = itime(1)
      iDay = itime(2)
      champfile = champdir+'/ElectronDensity/'+tostr(iYear)+'/'+ $
                  'CH-ME-2-PLP+'+tostr(iyear)+'-'+tostr(iMon,2)+'-'+$
                  tostr(iDay,2)+'_?.dat'

      print, 'Reading file : ',champfile
      read_champ_electron, champfile, champtime, champdata

      if (cDay eq 0) then begin
         ChampAlt = reform(champdata(0,*))-6378.0
         ChampRho = reform(champdata(3,*))
         ChampErr = reform(champdata(3,*))*0.5
         ChampTim = ChampTime
      endif else begin
         ChampAlt = [ChampAlt,reform(champdata( 0,*))-6378.0]
         ChampRho = [ChampRho,reform(champdata(3,*))]
         ChampErr = [ChampErr,reform(champdata(3,*))*0.5]
         ChampTim = [ChampTim,ChampTime]
      endelse

   endif else begin

      champfile = champdir+'Winds/'+sYear+'/Wind_3deg_'+$
                  tostr((iYear mod 100),2)+'_'+ $
                  chopr('00'+tostr(doy),3)+'.ascii'

      print, 'Reading '+champfile
      read_champ_winds, champfile, ChampTime, ChampData

      if (cDay eq 0) then begin
         ChampAlt = reform(champdata(0,*))
         ChampRho = reform(champdata(1,*))
         ChampEas = reform(champdata(2,*))
         ChampNor = reform(champdata(3,*))
         ChampErr = reform(champdata(4,*))
         ChampTim = ChampTime
      endif else begin
         ChampAlt = [ChampAlt,reform(champdata(0,*))]
         ChampRho = [ChampRho,reform(champdata(1,*))]
         ChampEas = [ChampEas,reform(champdata(2,*))]
         ChampNor = [ChampNor,reform(champdata(3,*))]
         ChampErr = [ChampErr,reform(champdata(4,*))]
         ChampTim = [ChampTim,ChampTime]
      endelse

   endelse

endfor

; Let's adjust the Champ Altitude, since GITM assumes a
; spherical Earth.
if (n_elements(ChampLat) gt 0) then begin
   re = 6378.137
   rp = 6356.75
   diff = re-rp
   rgitm = 6372.0
   r = rp+diff*cos(ChampLat*!dtor)
   correction = r - rgitm
   ChampAlt = ChampAlt + correction
endif

ChampTime = ChampTim

; interpolar on to gitm times

ChampIntAlt = fltarr(nTimes)
ChampIntRho = fltarr(nTimes)
ChampIntEas = fltarr(nTimes)
ChampIntNor = fltarr(nTimes)
ChampIntErr = fltarr(nTimes)
GitmIntRho  = fltarr(nTimes)
GitmIntEas  = fltarr(nTimes)
GitmIntNor  = fltarr(nTimes)

if (type eq 0) then factor = 1.0e12 else factor = 1.0

for itime = 0L, nTimes-1 do begin

   dt = GitmTime(iTime) - ChampTime
   loc = where(dt lt 0, c)
   champdt = 0.0
   if (c eq 0) then begin
           ; haven't reached the start of the file yet
      i = 1
      x = 1.0
   endif else begin
      if (GitmTime(iTime) gt max(ChampTime)) then begin
         i = n_elements(ChampTime)
         x = 0.0
      endif else begin
         i = loc(0)
         if (i eq 0) then i=1
         x = (champtime(i) - GitmTime(itime)) / (champtime(i)-champtime(i-1))
      endelse
   endelse

   ChampIntRho(iTime)   = $
      ((1-x) * ChampRho(i) + x*ChampRho(i-1))*factor
   if (type eq 2) then begin
      ChampIntEas(iTime)   = $
         ((1-x) * ChampEas(i) + x*ChampEas(i-1))*factor
      ChampIntNor(iTime)   = $
         ((1-x) * ChampNor(i) + x*ChampNor(i-1))*factor
   endif
   ChampIntAlt(iTime)  = $
      (1-x) * ChampAlt(i) + x * ChampAlt(i-1)
   ChampIntErr(iTime)      = $
      (1-x) * ChampErr(i) + x * ChampErr(i-1)

   loc = where(GitmAlts gt ChampIntAlt(iTime),c)
   if (c gt 1) then begin

      i = loc(0)
      x = (ChampIntAlt(iTime) - GitmAlts(i-1)) / $
          (GitmAlts(i) - GitmAlts(i-1))

      if (type eq 2) then begin

         GitmIntEas(iTime) = $
            (1.0 - x) * east(iTime,i-1) + $
            (      x) * east(iTIme,i)*factor

         GitmIntNor(iTime) = $
            (1.0 - x) * north(iTime,i-1) + $
            (      x) * north(iTIme,i)*factor

      endif else begin

         GitmIntRho(iTime) = $
            exp((1.0 - x) * alog(GitmRho(iTime,i-1)) + $
                (      x) * alog(GitmRho(iTIme,i)))*factor

      endelse

   endif else begin

      if (type eq 2) then begin

         GitmIntEas(iTime) = east(iTime)*factor
         GitmIntNor(iTime) = north(iTime)*factor

      endif else begin

         GitmIntRho(iTime) = GitmRho(iTime)*factor

      endelse

   endelse

endfor

if (type eq 2) then begin

   ;GitmIntRho = (GitmIntNor * ChampIntNor + GitmIntEas * ChampIntEas)
   GitmIntRho = -GitmIntEas * abs(ChampIntEas)
   GitmIntRho = GitmIntRho -GitmIntNor * abs(ChampIntNor)

endif

d = [ChampIntRho,GitmIntRho]
maxi = max(abs(d))
n = n_elements(d)
l = where(abs(d) gt maxi,c)
while c lt 0.05*n do begin
   maxi = maxi*0.99
   l = where(abs(d) gt maxi,c)
endwhile

maxi = maxi*1.2

if (type eq 2) then yrange = [-maxi,maxi] else yrange = [0,maxi]

ppp = 3
space = 0.04
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 0, pos0, /rect
pos0(0) = pos0(0) + 0.05
pos0(2) = pos0(2) - 0.05

get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

get_position, ppp, space, sizes, 2, pos2, /rect
pos2(0) = pos2(0) + 0.05
pos2(2) = pos2(2) - 0.05

stime = min(GitmTime)
etime = max(GitmTime)
time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn

c_r_to_a, itime, stime
c_a_to_y_m_d, itime, ymd
c_r_to_a, itime, etime
c_a_to_y_m_d, itime, ymd2
psfile_raw  = SatType+'_raw_'
psfile_ave  = SatType+'_ave_'
psfile_mapa = SatType+'_mapa_'
psfile_mapd = SatType+'_mapd_'

if (type eq 0) then begin
   psfile_raw = psfile_raw+'rho'
   psfile_ave = psfile_ave+'rho'
   psfile_mapa = psfile_mapa+'rho'
   psfile_mapd = psfile_mapd+'rho'
endif

if (type eq 1) then begin
   psfile_raw = psfile_raw+'electron'
   psfile_ave = psfile_ave+'electron'
   psfile_mapa = psfile_mapa+'electron'
   psfile_mapd = psfile_mapd+'electron'
endif

if (type eq 2) then begin
   psfile_raw = psfile_raw+'wind'
   psfile_ave = psfile_ave+'wind'
   psfile_mapa = psfile_mapa+'wind'
   psfile_mapd = psfile_mapd+'wind'
endif

psfile_raw = psfile_raw+'_'+ymd+'-'+ymd2+'.ps'
psfile_ave = psfile_ave+'_'+ymd+'-'+ymd2+'.ps'
psfile_mapa = psfile_mapa+'_'+ymd+'-'+ymd2+'.ps'
psfile_mapd = psfile_mapd+'_'+ymd+'-'+ymd2+'.ps'

AltMean = mean(ChampIntAlt)

Re = 6378.0D
r = (Re + AltMean) * 1000.0D
G = 6.67259d-11
Me = 5.9722d24                  ; reference 1
mu = G * Me
v0 = sqrt(mu/r)
period = 2*!pi*r / v0 / 60.0
print, period

GitmAve = GitmIntRho*0.0
ChampAve = ChampIntRho*0.0

for iTime = 0L, nTimes-1 do begin
   loc = where(abs(gitmtime-gitmtime(iTime)) lt period*60.0/2, count)
   if (count gt 0) then begin
      GitmAve(iTime) = mean(GitmIntRho(loc))
      ChampAve(iTime) = mean(ChampIntRho(loc))
   endif
endfor

dLat = 5.0
nX = fix((etime-stime)/(period*60.0))+2
nY = fix(180.0/dlat) + 1

GitmMapa = fltarr(nX, nY)
ChampMapa = fltarr(nX, nY)

GitmMapd = fltarr(nX, nY)
ChampMapd = fltarr(nX, nY)

x = round((GitmTime-sTime)/(period*60.0))
y = round((Lats+90.0)/dLat)
dl = Lats(1:nTimes-1)-Lats(0:nTimes-2)

maptimea = fltarr(nX, nY)
maplat   = fltarr(nX, nY)
maptimed = fltarr(nX, nY)
maplat   = fltarr(nX, nY)

maplat = fltarr(nX,nY)
for i=0,nX-1 do maplat(i,*) = findgen(nY)*dLat-90.0 + dLat/2.0

l = where(dl gt 0 and abs(lats) lt 45.0, c)
if (c gt 0) then lta = mean(localtime(l))

l = where(dl lt 0 and abs(lats) lt 45.0, c)
if (c gt 0) then ltd = mean(localtime(l))

for i=0L,nTimes-2 do begin
   if (dl(i) gt 0) then begin
      GitmMapa(x(i),y(i)) = gitmIntRho(i)
      ChampMapa(x(i),y(i)) = champIntRho(i)
      maptimea(x(i),y(i)) = gitmtime(i)
   endif else begin
      GitmMapd(x(i),y(i)) = gitmIntRho(i)
      ChampMapd(x(i),y(i)) = champIntRho(i)
      maptimed(x(i),y(i)) = gitmtime(i)
   endelse
endfor

for i=0,nY-1 do maptimea(*,i) = dindgen(nX)*period*60.0
for i=0,nY-1 do maptimed(*,i) = dindgen(nX)*period*60.0

; -------------------------------
; Ascending Map
; -------------------------------

setdevice, psfile_mapa, 'p', 5

makect,'mid'

if (type eq 2) then levels = 2*findgen(31)*maxi/30.0-maxi $
else levels = findgen(31)*maxi*1.1/30.0

l = where(champmapa gt levels(29),c)
if (c gt 0) then champmapa(l) = levels(29)
l = where(gitmmapa gt levels(29),c)
if (c gt 0) then gitmmapa(l) = levels(29)

l = where(champmapa lt levels(1),c)
if (c gt 0) then champmapa(l) = levels(1)
l = where(gitmmapa lt levels(1),c)
if (c gt 0) then gitmmapa(l) = levels(1)

cPos = pos0
cPos(2) = cPos(2)-0.05
contour, champmapa, maptimea, maplat, pos = cPos, $
         xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
         xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
         ytitle = 'Latitude (deg)', thick = 3, $
         xrange = [btr,etr], levels = levels, /fill

xyouts, cPos(2), cPos(3)+0.01, 'Local Time : '+string(lta,format='(f5.2)'),$
        /norm, align = 1.0

ctpos = cPos
ctpos(0) = cPos(2)+0.01
ctpos(2) = ctpos(0)+0.03

if (type eq 0) then title = ' Mass Density (10!E-12!N kg/m!E3!N)'
if (type eq 1) then title = ' Electron Density (/m!E3!N)'
if (type eq 2) then title = ' Cross Track Wind (m/s)'
plotct,254,ctpos,mm(levels),Satellite+title,/right

cPos = pos
cPos(2) = cPos(2)-0.05
contour, gitmmapa, maptimea, maplat, pos = cPos, $
         xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
         xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
         ytitle = 'Latitude (deg)', thick = 3, /noerase, $
         xrange = [btr,etr], levels = levels, /fill, $
         xtitle = xtitle

;linelevels = (findgen(9)+1.0)/10.0 * maxi *1.1
;
;contour, champmapa, maptimea, maplat, pos = cPos, $
;         xstyle = 5, thick = 1, /noerase, ystyle = 5, $
;         xrange = [btr,etr], levels = linelevels

ctpos = cPos
ctpos(0) = cPos(2)+0.01
ctpos(2) = ctpos(0)+0.03

plotct,254,ctpos,mm(levels),'GITM'+title,/right

closedevice

; -------------------------------
; Descending Map
; -------------------------------

l = where(champmapd gt levels(29),c)
if (c gt 0) then champmapd(l) = levels(29)
l = where(gitmmapd gt levels(29),c)
if (c gt 0) then gitmmapd(l) = levels(29)

l = where(champmapd lt levels(1),c)
if (c gt 0) then champmapd(l) = levels(1)
l = where(gitmmapd lt levels(1),c)
if (c gt 0) then gitmmapd(l) = levels(1)

setdevice, psfile_mapd, 'p', 5

makect,'mid'

cPos = pos0
cPos(2) = cPos(2)-0.05
contour, champmapd, maptimed, maplat, pos = cPos, $
         xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
         xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
         ytitle = 'Latitude (deg)', thick = 3, $
         xrange = [btr,etr], levels = levels, /fill

xyouts, cPos(2), cPos(3)+0.01, 'Local Time : '+string(ltd,format='(f5.2)'), $
        /norm, align = 1.0

ctpos = cPos
ctpos(0) = cPos(2)+0.01
ctpos(2) = ctpos(0)+0.03

plotct,254,ctpos,mm(levels),Satellite+title,/right

cPos = pos
cPos(2) = cPos(2)-0.05
contour, gitmmapd, maptimed, maplat, pos = cPos, $
         xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
         xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
         ytitle = 'Latitude (deg)', thick = 3, /noerase, $
         xrange = [btr,etr], levels = levels, /fill, $
         xtitle = xtitle

;linelevels = (findgen(9)+1.0)/10.0 * maxi * 1.1
;
;contour, champmapd, maptimed, maplat, pos = cPos, $
;         xstyle = 5, thick = 1, /noerase, ystyle = 5, $
;         xrange = [btr,etr], levels = linelevels

ctpos = cPos
ctpos(0) = cPos(2)+0.01
ctpos(2) = ctpos(0)+0.03

plotct,254,ctpos,mm(levels),'GITM'+title,/right

closedevice

setdevice, psfile_raw, 'p', 5

if (type eq 0) then ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)'
if (type eq 1) then ytitle = 'Electron Density (/m!E3!N)'
if (type eq 2) then ytitle = 'Cross Track Wind (m/s)'

plot, Gitmtime-stime, ChampIntRho, yrange = yrange*1.2, pos = pos, $
      xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = ytitle,   $
      thick = 3, /noerase
oplot, gitmtime-stime, gitmIntRho, linestyle = 2, thick = 3
    
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

d = ChampIntRho-gitmIntRho 
yr = [-maxi,maxi]

plot, Gitmtime-stime, d, $
      yrange = yr, pos = pos2, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = ytitle,   $
      thick = 3, /noerase, ystyle = 1
    
rmse = sqrt(mean((ChampIntRho-gitmIntRho)^2))
rmsd = sqrt(mean((ChampIntRho)^2))
nrms = rmse/rmsd * 100.0
    
pdif = mean((ChampIntRho-gitmIntRho)/ChampIntRho) * 100.0
    
srms = +' (nRMS: '+string(nrms,format = '(f5.1)')+'%, '
srms = srms+string(pdif,format = '(f5.1)')+'% Difference)'
    
xyouts, pos2(2)+0.01, (pos2(1)+pos2(3))/2.0, srms, $
        alignment = 0.5, orient=270, /norm 

xyouts, pos2(2)+0.035, (pos2(1)+pos2(3))/2.0, Satellite+' - GITM', $
        alignment = 0.5, orient=270, /norm 

oplot, [btr,etr], [0.0,0.0], linestyle = 1
    
spawn, 'pwd', dir
dir = dir(0)
xyouts, 0.0, -0.05, dir, /norm, charsize = 0.8
    
closedevice

setdevice, psfile_ave, 'p', 5

plot, Gitmtime-stime, ChampAve, yrange = yrange*1.2, pos = pos, $
      xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = ytitle,   $
      thick = 3, /noerase
oplot, gitmtime-stime, gitmAve, linestyle = 2, thick = 3
    
t1 = (etr-btr)*0.05
t2 = (etr-btr)*0.10
t3 = (etr-btr)*0.11

ChampMedian = median(ChampAve)
GitmMedian = median(GitmAve)
oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3, linestyle = 2
xyouts, t3, max(yrange)*1.05, 'GITM (Median: '+ $
        string(GitmMedian,format='(f5.3)')+' x10!U-12!N kg/m!U3!N)'
    
t1 = (etr-btr)*0.55
t2 = (etr-btr)*0.60
t3 = (etr-btr)*0.61
    
oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3
xyouts, t3, max(yrange)*1.05, Satellite+' ('+ $
        string(ChampMedian,format='(f5.3)')+' x10!U-12!N kg/m!U3!N)'

d = ChampAve-gitmAve 
yr = [-maxi,maxi]

plot, Gitmtime-stime, d, $
      yrange = yr, pos = pos2, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = ytitle,   $
      thick = 3, /noerase, ystyle = 1
    
rmse = sqrt(mean((ChampAve-gitmAve)^2))
rmsd = sqrt(mean((ChampAve)^2))
nrms = rmse/rmsd * 100.0
    
pdif = mean((ChampAve-gitmAve)/ChampAve) * 100.0
    
srms = +' (nRMS: '+string(nrms,format = '(f5.1)')+'%, '
srms = srms+string(pdif,format = '(f5.1)')+'% Difference)'
    
xyouts, pos2(2)+0.01, (pos2(1)+pos2(3))/2.0, srms, $
        alignment = 0.5, orient=270, /norm 

xyouts, pos2(2)+0.035, (pos2(1)+pos2(3))/2.0, Satellite+' - GITM', $
        alignment = 0.5, orient=270, /norm 

oplot, [btr,etr], [0.0,0.0], linestyle = 1
    
spawn, 'pwd', dir
dir = dir(0)
xyouts, 0.0, -0.05, dir, /norm, charsize = 0.8
    
closedevice

end

