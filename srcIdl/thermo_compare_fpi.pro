
pro plot_values, TimesNe, TimesSw, TimesMinMax, $
                 FpiNe, FpiSw, GitmNe, GitmSw, $
                 Label, plotnum, minmax=minmax

  ymin = 0.05+(4-plotnum)*0.20
  ymax = ymin+0.19
  pos = [0.1,ymin,0.95,ymax]

  if (n_elements(minmax) eq 0) then begin
     All = [FpiNe,FpiSw,GitmNe,GitmSw]
     l = where(abs(All) lt 2000.0,c)
     if (c gt 0) then begin
        if (min(All(l)) lt 0) then begin
           maxi = max(abs(All(l)))*1.1
           mini = -maxi
        endif else begin
           maxi = max(All(l))*1.05
           mini = min(All(l))*0.90
        endelse
     endif
  endif else begin
     mini = min(minmax)
     maxi = max(minmax)
  endelse

  stime = min(TimesMinMax)
  etime = max(TimesMinMax)
  time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn
  
  if plotnum lt 4 then begin
     xtn = strarr(60)+' '
     xt  = ' '
  endif else begin
     xtn = xtickname
     xt  = xtitle
  endelse

  plot, [btr,etr], [mini,maxi], pos = pos, /noerase, /nodata, $
        xtickname = xtn,			$
        xtitle = xt,			$
        xtickv = xtickv,			$
        xminor = xminor,			$
        xticks = xtickn,   $
        ytitle = Label, $
        ystyle = 1, xstyle = 1

  oplot, [btr,etr], [0.0,0.0], linestyle = 1, thick = 3
  
  oplot, TimesNe, FpiNe, min_value = mini, linestyle = 0, thick = 5
  oplot, TimesNe, GitmNe, min_value = mini, linestyle = 0, thick = 5, color = 200

  if (TimesSw(0) ne TimesNe(0)) then begin
     oplot, TimesSw, FpiSw, min_value = mini, linestyle = 2, thick = 5
     oplot, TimesSw, GitmSw, min_value = mini, linestyle = 2, thick = 5, color = 200
     s = strpos(label,'/')
  endif else s = strpos(label,' ')

  e = strmid(label,0,s)
  w = strmid(label,s+1,5)

  r = (maxi-mini)*0.05
  y = mini+[r,r]

  r = (etr-btr)*0.05
  x = btr + [r,2*r]
  oplot, x, y, thick = 5
  xyouts, x(1)+0.2*r, y(1), 'FPI '+e, charsize = 0.9

  x = btr + [5*r,6*r]
  oplot, x, y, thick = 5, color = 200
  xyouts, x(1)+0.2*r, y(1), 'GITM '+e, charsize = 0.9, color = 200

  if (TimesSw(0) ne TimesNe(0)) then begin
     x = btr + [10*r,11*r]
     oplot, x, y, thick = 5, linestyle = 2
     xyouts, x(1)+0.2*r, y(1), 'FPI '+w, charsize = 0.9

     x = btr + [14*r,15*r]
     oplot, x, y, thick = 5, color = 200, linestyle = 2
     xyouts, x(1)+0.2*r, y(1), 'GITM '+w, charsize = 0.9, color = 200
  endif
end

pro write_info, unit, label, data, loc

  nPts = n_elements(loc)
  printf,unit,label
  if (loc(0) ge 0) then $
     for i=0,nPts-1 do printf,unit,data(loc(i)) 
  printf,unit,''

end

FpiAlt = 245.0
DistToDeg = 360.0/((6372.0+FpiAlt)*2*!pi)

Stations = ['Arecibo','Cajazeiras','Cariri',$
            'Jicamarca','Kentucky','Millstone',$
            'Movil','Urbana','PeachMountain',$
            'Pari']
centerlats = [ 18.344167, -6.88, -7.38,-11.96, 37.75, 42.619,-14.97, 40.13, 42.27, 35.20]
centerlons = [-66.752778,-38.56,-36.52,-76.86,-84.29,-79.491,-74.89,-88.20,-83.75,-82.85]+360.0

print, 'Where is the FPI located?'
display, Stations
if (n_elements(loc) eq 0) then loc = '0'
loc = ask('location ',loc)
iLoc = fix(loc)

centerlat = centerlats(iLoc)
centerlon = centerlons(iLoc)

if (n_elements(whichyear) eq 0) then whichyear = '13'
whichyear = ask('which year you would like to compare', whichyear)
iYear = fix(whichyear) mod 1000
sYear = tostr(iYear,2)

if (n_elements(whichmonth) eq 0) then whichmonth = '03'
whichmonth = ask('which month you would like to compare', whichmonth)
iMonth = fix(whichmonth)
sMonth = tostr(iMonth,2)

if (n_elements(whichday) eq 0) then whichday = '5'
whichday = ask('which day you would like to compare', whichday)
iDay = fix(whichday)
sDay = tostr(iDay,2)

MadrigalFileList = '/raid3/Gitm/Runs/FabryPerot/'+Stations(iLoc)+'/*'+sYear+sMonth+sDay+'*'
print, 'Looking for file : ',MadrigalFileList
filelist = findfile(MadrigalFileList)
MadrigalFile = filelist(0)
if (strlen(MadrigalFile) eq 0) then begin
   print, "Can't find file! Please download and put in the right place!"
   stop
endif

read_madrigal, MadrigalFile, MadrigalTimes, MadrigalData, MadrigalVars

iEl = -1
iAz = -1
iTemp = -1
iTempE = -1
iLos = -1
iLosE = -1

nMadrigalVars = n_elements(MadrigalVars)
CorrectForHorizontal = 0
for iVar = 0, nMadrigalVars-1 do begin

   if (strpos(MadrigalVars(iVar),'ELM') eq 0)   then iEl    = iVar
   if (strpos(MadrigalVars(iVar),'AZM') eq 0)   then iAz    = iVar

   if (strpos(MadrigalVars(iVar),'TN') eq 0)    then iTemp  = iVar
   if (strpos(MadrigalVars(iVar),'DTN') eq 0)   then iTempE = iVar

   if (strpos(MadrigalVars(iVar),'VNU') eq 0)   then begin
      iLos   = iVar
      CorrectForHorizontal = 1
   endif
   if (strpos(MadrigalVars(iVar),'VNHLU') eq 0) then iLos   = iVar

   if (strpos(MadrigalVars(iVar),'DVNU') eq 0)  then iLosE  = iVar
   if (strpos(MadrigalVars(iVar),'DVNHLU') eq 0)then iLosE  = iVar

endfor

iError = 0
if (iEl eq -1) then begin
   print, "Can't find elevation variable in Madrigal file."
   iError = 1
endif
if (iAz eq -1) then begin
   print, "Can't find Azimuth variable in Madrigal file."
   iError = 1
endif
if (iTemp eq -1) then begin
   print, "Can't find Temperature variable in Madrigal file."
   iError = 1
endif
if (iTempE eq -1) then begin
   print, "Can't find Temperature Error variable in Madrigal file."
   iError = 1
endif
if (iLos eq -1) then begin
   print, "Can't find LOS variable in Madrigal file."
   iError = 1
endif
if (iLosE eq -1) then begin
   print, "Can't find LOS Error variable in Madrigal file."
   iError = 1
endif

if (iError) then begin
   print, "Here is a list of variables:"
   display, MadrigalVars
   stop
endif

FpiAz   = reform(MadrigalData(iAz,*))
FpiEl   = reform(MadrigalData(iEl,*))
FpiLos  = reform(MadrigalData(iLos,*))
FpiTemp = reform(MadrigalData(iTemp,*))

nTimes  = n_elements(MadrigalTimes)

; Find the location of each measurement point:
FpiLats = fltarr(nTimes)+centerlat
FpiLons = fltarr(nTimes)+centerlon
FpiAlts = fltarr(nTimes)+FpiAlt

; If the FPI is not looking up, then the measurement is at a
; different lat/lon: 
l = where(FpiEl lt 89.0,c)
if (c gt 0) then begin
   FpiLats(l) = centerlat+FpiAlts(l)/tan(FpiEl(l)*!dtor)*cos(FpiAz(l)*!dtor)*DistToDeg
   FpiLons(l) = centerlon+FpiAlts(l)/tan(FpiEl(l)*!dtor)*sin(FpiAz(l)*!dtor)*DistToDeg
endif

; Find Cardinal Direction Points:
LocNorth = where(FpiEl lt 85.0 and cos(FpiAz*!dtor) gt 0.99, nNorth)
LocSouth = where(FpiEl lt 85.0 and cos(FpiAz*!dtor+!pi) gt 0.99, nSouth)
LocEast  = where(FpiEl lt 85.0 and cos(FpiAz*!dtor-!pi/2) gt 0.99, nEast)
LocWest  = where(FpiEl lt 85.0 and cos(FpiAz*!dtor+!pi/2) gt 0.99, nWest)
LocVertical  = where(FpiEl gt 85.0, nVertical)

nCardinal = nNorth+nSouth+nEast+nWest+nVertical

if (nCardinal lt nTimes) then begin
   print, "Non-cardinal directions found!"
   ; There are a lot of rounding errors here, so let's try to do this
   ; within about 5 degrees of the different values
   DummyAz = float(fix(FpiAz/5.0+0.5))*5.0
   AllAzimuthsRounded = DummyAz
   if (nNorth gt 0) then DummyAz(LocNorth) = -1.0e32
   if (nSouth gt 0) then DummyAz(LocSouth) = -1.0e32
   if (nEast gt 0) then DummyAz(LocEast) = -1.0e32
   if (nWest gt 0) then DummyAz(LocWest) = -1.0e32
   if (nVertical gt 0) then DummyAz(LocVertical) = -1.0e32
   IsFirstTime = 1
   while (max(DummyAz) gt -1.0e32) do begin
      l = where(DummyAz eq max(DummyAz))
      if (IsFirstTime) then UniqueAz = max(DummyAz) $
      else UniqueAz = [UniqueAz,max(DummyAz)]
      DummyAz(l) = -1.0e32
      IsFirstTime = 0
   endwhile
   nUnique = n_elements(UniqueAz)
   print, "Number of unique azimuths found: ",nUnique
endif else nUnique = 0

; Correct for Horizontal in some stations
if (CorrectForHorizontal) then begin
   print, "Correcting for horizontal projection!"
   l = where(FpiEl lt 85.0 and FpiLos gt -1.0e31,c)
   if (c gt 0) then FpiLos(l) = FpiLos(l)/cos(FpiEl(l)*!dtor)
   ; If we are looking South, then a positive LOS is North!
   if (nSouth gt 0) then FpiLos(LocSouth) = -FpiLos(LocSouth)
   ; If we are looking West, then a positive LOS is East!
   if (nWest gt 0) then FpiLos(LocWest) = -FpiLos(LocWest)
endif else begin
   ; If we are looking South, then a positive LOS is North!
   if (nNorth gt 0) then FpiLos(LocNorth) = -FpiLos(LocNorth)
   ; If we are looking West, then a positive LOS is East!
   if (nEast gt 0) then FpiLos(LocEast) = -FpiLos(LocEast)
;   ; If we are looking South, then a positive LOS is North!
;   if (nSouth gt 0) then FpiLos(LocSouth) = -FpiLos(LocSouth)
;   ; If we are looking West, then a positive LOS is East!
;   if (nWest gt 0) then FpiLos(LocWest) = -FpiLos(LocWest)
endelse

; shit, we may have turned some -1.0e32 into 1.0e32
l = where(FpiLos gt 1.0e31,c)
if (c gt 0) then FpiLos(l) = -FpiLos(l)

l = where(FpiLos gt -1.0e31,c)
if (c gt 0) then begin
   maxiLos = max(abs(FpiLos(l)))*1.2
   miniLos = -maxiLos
endif else begin
   miniLos = -1.0
   maxiLos = 1.0
endelse

if (n_elements(LosMaxi) eq 0) then LosMaxi = maxiLos
if (n_elements(LosMini) eq 0) then LosMini = miniLos

LosMini = float(ask('min value for LOS',string(LosMini)))
LosMaxi = float(ask('max value for LOS',string(LosMaxi)))

l = where(FpiTemp gt -1.0e31,c)
if (c gt 0) then begin
   maxiTemp = max(FpiTemp(l))*1.2
   miniTemp = min(FpiTemp(l))*0.9
endif else begin
   miniTemp = 0.0
   maxiTemp = 1.0
endelse

if (n_elements(TempMaxi) eq 0) then TempMaxi = maxiTemp
if (n_elements(TempMini) eq 0) then TempMini = miniTemp

TempMini = float(ask('min value for Temp',string(TempMini)))
TempMaxi = float(ask('max value for Temp',string(TempMaxi)))

VarsToGet = [15,16,17,18]

get_gitm_points, MadrigalTimes, FpiLons, FpiLats, FpiAlts, VarsToGet, GitmData

Ue       = reform(GitmData(*,1))
Un       = reform(GitmData(*,2))
Uv       = reform(GitmData(*,3))
GitmTemp = reform(GitmData(*,0))

; assume that the vertical component is zero at this point:
GitmLos = Un*cos(FpiAz*!dtor) + Ue*sin(FpiAz*!dtor)
; except in the case of actually being vertical:
if (nVertical gt 0) then GitmLos(LocVertical) = Uv(LocVertical)

; Now, we have to correct this for the cardinal directions, since
; we just computed the LOS and we need to actual direction, which
; is somewhat stupid.  The easiest way to do this is:
if (nNorth gt 0) then GitmLos(LocNorth) = Un(LocNorth)
if (nSouth gt 0) then GitmLos(LocSouth) = Un(LocSouth)
if (nEast gt 0) then GitmLos(LocEast) = Ue(LocEast)
if (nWest gt 0) then GitmLos(LocWest) = Ue(LocWest)
; I know.  Stupid.  Why even compute the LOS, right?

sTime = min(MadrigalTimes)
c_r_to_a, itimeStart, sTime
dTimes = MadrigalTimes-stime

outfile = Stations(iLoc)+'_GITM_'+sYear+sMonth+sDay+'_compare.dat'
psfile  = Stations(iLoc)+'_GITM_'+sYear+sMonth+sDay+'_compare.ps'

close,2
openw,2,outfile

printf,2,itimeStart
printf,2,''

write_info, 2, 'FPI North Times', dTimes, LocNorth
write_info, 2, 'FPI South Times', dTimes, LocSouth
write_info, 2, 'FPI East Times',  dTimes, LocEast
write_info, 2, 'FPI West Times',  dTimes, LocWest
write_info, 2, 'FPI Vertical Times',  dTimes, LocVertical
if (nUnique gt 0) then begin
   for i=0,nUnique-1 do begin
      LocUnique  = where(AllAzimuthsRounded eq UniqueAz(i))
      write_info, 2, 'FPI Azimuth '+tostr(UniqueAz(i))+' Times', dTimes, LocUnique
   endfor
endif

write_info, 2, 'FPI North Winds', FpiLos, LocNorth
write_info, 2, 'FPI South Winds', FpiLos, LocSouth
write_info, 2, 'FPI East Winds',  FpiLos, LocEast
write_info, 2, 'FPI West Winds',  FpiLos, LocWest
write_info, 2, 'FPI Vertical Winds',  FpiLos, LocVertical
if (nUnique gt 0) then begin
   for i=0,nUnique-1 do begin
      LocUnique  = where(AllAzimuthsRounded eq UniqueAz(i))
      write_info, 2, 'FPI Azimuth '+tostr(UniqueAz(i))+' Winds', FpiLos, LocUnique
   endfor
endif

write_info, 2, 'FPI North Temps', FpiTemp, LocNorth
write_info, 2, 'FPI South Temps', FpiTemp, LocSouth
write_info, 2, 'FPI East Temps',  FpiTemp, LocEast
write_info, 2, 'FPI West Temps',  FpiTemp, LocWest
write_info, 2, 'FPI Vertical Temps',  FpiTemp, LocVertical
if (nUnique gt 0) then begin
   for i=0,nUnique-1 do begin
      LocUnique  = where(AllAzimuthsRounded eq UniqueAz(i))
      write_info, 2, 'FPI Azimuth '+tostr(UniqueAz(i))+' Temps', FpiTemp, LocUnique
   endfor
endif

write_info, 2, 'GITM North Winds', GitmLos, LocNorth
write_info, 2, 'GITM South Winds', GitmLos, LocSouth
write_info, 2, 'GITM East Winds',  GitmLos, LocEast
write_info, 2, 'GITM West Winds',  GitmLos, LocWest
write_info, 2, 'GITM Vertical Winds',  GitmLos, LocVertical
if (nUnique gt 0) then begin
   for i=0,nUnique-1 do begin
      LocUnique  = where(AllAzimuthsRounded eq UniqueAz(i))
      write_info, 2, 'GITM Azimuth '+tostr(UniqueAz(i))+' Winds', GitmLos, LocUnique
   endfor
endif

write_info, 2, 'GITM North Temps', GitmTemp, LocNorth
write_info, 2, 'GITM South Temps', GitmTemp, LocSouth
write_info, 2, 'GITM East Temps',  GitmTemp, LocEast
write_info, 2, 'GITM West Temps',  GitmTemp, LocWest
write_info, 2, 'GITM Vertical Temps',  GitmTemp, LocVertical
if (nUnique gt 0) then begin
   for i=0,nUnique-1 do begin
      LocUnique  = where(AllAzimuthsRounded eq UniqueAz(i))
      write_info, 2, 'GITM Azimuth '+tostr(UniqueAz(i))+' Temps', GitmTemp, LocUnique
   endfor
endif

close,2

setdevice, psfile, 'p', 5
plotdumb
makect,'mid'

if (nNorth gt 0) then begin
   TimesNe = dTimes(LocNorth)
   FpiNe = FpiLos(LocNorth)
   GitmNe = GitmLos(LocNorth)
   FpiTempNe = FpiTemp(LocNorth)
   GitmTempNe = GitmTemp(LocNorth)
endif else begin
   TimesNe = [0,0]
   FpiNe = [-1.0,1.0]
   GitmNe = [-1.0,1.0]
   FpiTempNe = [0.0,1.0]
   GitmTempNe = [0.0,1.0]
endelse

if (nSouth gt 0) then begin
   TimesSw = dTimes(LocSouth)
   FpiSw = FpiLos(LocSouth)
   GitmSw = GitmLos(LocSouth)
   FpiTempSw = FpiTemp(LocSouth)
   GitmTempSw = GitmTemp(LocSouth)
endif else begin
   TimesSw = [0,0]
   FpiSw = [-1.0,1.0]
   GitmSw = [-1.0,1.0]
   FpiTempSw = [0.0,1.0]
   GitmTempSw = [0.0,1.0]
endelse

Label = 'North/South LOS (m/s)'
plotnum = 0
TimesMinMax = [min(MadrigalTimes),max(MadrigalTimes)]
plot_values, TimesNe, TimesSw, TimesMinMax, $
             FpiNe, FpiSw, GitmNe, GitmSw, $
             Label, plotnum, minmax = [LosMini, LosMaxi]

Label = 'North/South Temp (K)'
plotnum = 2
TimesMinMax = [min(MadrigalTimes),max(MadrigalTimes)]
plot_values, TimesNe, TimesSw, TimesMinMax, $
             FpiTempNe, FpiTempSw, $
             GitmTempNe, GitmTempSw, $
             Label, plotnum, minmax = [TempMini, TempMaxi]

if (nEast gt 0) then begin
   TimesNe = dTimes(LocEast)
   FpiNe = FpiLos(LocEast)
   GitmNe = GitmLos(LocEast)
   FpiTempNe = FpiTemp(LocEast)
   GitmTempNe = GitmTemp(LocEast)
endif else begin
   TimesNe = [0,0]
   FpiNe = [-1.0,1.0]
   GitmNe = [-1.0,1.0]
   FpiTempNe = [0.0,1.0]
   GitmTempNe = [0.0,1.0]
endelse

if (nWest gt 0) then begin
   TimesSw = dTimes(LocWest)
   FpiSw = FpiLos(LocWest)
   GitmSw = GitmLos(LocWest)
   FpiTempSw = FpiTemp(LocWest)
   GitmTempSw = GitmTemp(LocWest)
endif else begin
   TimesSw = [0,0]
   FpiSw = [-1.0,1.0]
   GitmSw = [-1.0,1.0]
   FpiTempSw = [0.0,1.0]
   GitmTempSw = [0.0,1.0]
endelse

Label = 'East/West LOS (m/s)'
plotnum = 1
TimesMinMax = [min(MadrigalTimes),max(MadrigalTimes)]
plot_values, TimesNe, TimesSw, TimesMinMax, $
             FpiNe, FpiSw, GitmNe, GitmSw, $
             Label, plotnum, minmax = [LosMini, LosMaxi]

Label = 'East/West Temp (K)'
plotnum = 3
TimesMinMax = [min(MadrigalTimes),max(MadrigalTimes)]
plot_values, TimesNe, TimesSw, TimesMinMax, $
             FpiTempNe, FpiTempSw, GitmTempNe, GitmTempSw, $
             Label, plotnum, minmax = [TempMini, TempMaxi]

if (nVertical gt 0) then begin
   TimesNe = dTimes(LocVertical)
   FpiTempNe = FpiTemp(LocVertical)
   GitmTempNe = GitmTemp(LocVertical)
endif else begin
   TimesNe = [0,0]
   FpiTempNe = [800.0,900.0]
   GitmTempNe = [800.0,900.0]
endelse

Label = 'Vertical Temp (K)'
plotnum = 4
TimesSw = TimesNe
FpiTempSw = FpiTempNe
GitmTempSw = GitmTempNe

TimesMinMax = [min(MadrigalTimes),max(MadrigalTimes)]
plot_values, TimesNe, TimesSw, TimesMinMax, $
             FpiTempNe, FpiTempSw, GitmTempNe, GitmTempSw, $
             Label, plotnum, minmax = [TempMini, TempMaxi]

xyouts, 0.5, 1.01, Stations(iLoc), align = 0.5, /norm, charsize = 1.1

closedevice

end


