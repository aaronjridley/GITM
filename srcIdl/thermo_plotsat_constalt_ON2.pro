
GetNewData = 1

if (n_elements(data) gt 0) then begin
    filelist_new = findfile("GUVI*/????_*.dat")
    nfiles_new = n_elements(filelist_new)
    if (nfiles_new eq nfiles) then default = 'n' else default='y'
    GetNewData = mklower(strmid(ask('whether to reread data',default),0,1))
    if (GetNewData eq 'n') then GetNewData = 0 else GetNewData = 1
endif else filelist_new = findfile("GUVI*/????_*.dat")


if (GetNewData) then begin
    thermo_readsat, filelist_new,data, time, nTimes, Vars, nAlts, nSats, nFiles
endif

nPts = nTimes

GitmAlts = reform(data(0,0,2,0:nalts-1))/1000.0
Alts = reform(data(0,0:nPts-1,2,0:nalts-1))/1000.0
Lons = reform(data(0,0:nPts-1,0,0)) * 180.0 / !pi
Lats = reform(data(0,0:nPts-1,1,0)) * 180.0 / !pi
Lat_grad= reform(lats(1:nPts-1)-lats(0:nPts-2))

t  = reform(data(0,0:nPts-1,4,0:nalts-1))

o  = reform(data(0,0:nPts-1,5,0:nalts-1))
o2 = reform(data(0,0:nPts-1,6,0:nalts-1))
n2 = reform(data(0,0:nPts-1,7,0:nalts-1))
n4s = reform(data(0,0:nPts-1,9,0:nalts-1))
n = o + n2 + o2 + n4s
k = 1.3807e-23
mp = 1.6726e-27
rho = o*mp*16 + o2*mp*32 + n2*mp*14
data(0,0:nPts-1,3,0:nalts-1) = rho


    oInt = fltarr(nPts)
    n2Int = fltarr(nPts)
    AltInt = fltarr(nPts)
    ratio =fltarr(nPts)
    MaxValN2 = 1.0e21

    for i=0,nPts-1 do begin

        iAlt = nalts-1
        Done = 0
        while (Done eq 0) do begin
            dAlt = (Alts(i,iAlt)-Alts(i,iAlt-1))*1000.0
            n2Mid = (n2(i,iAlt) + n2(i,iAlt-1))/2.0
            oMid  = ( o(i,iAlt) +  o(i,iAlt-1))/2.0
            if (n2Int(i) + n2Mid*dAlt lt MaxValN2) then begin
                n2Int(i) = n2Int(i) + n2Mid*dAlt
                oInt(i)  =  oInt(i) +  oMid*dAlt
                iAlt = iAlt - 1
            endif else begin
                dAlt = (MaxValN2 - n2Int(i)) / n2Mid
                n2Int(i) = n2Int(i) + n2Mid*dAlt
                oInt(i)  =  oInt(i) +  oMid*dAlt
                AltInt(i) = Alts(i,iAlt) - dAlt/1000.0
                Done = 1
            endelse
        endwhile
        ratio(i)=oInt(i)/n2Int(i)
    endfor

;une = reform(data(0,0:nPts-1,11,0:nalts-1))
;unn = reform(data(0,0:nPts-1,12,0:nalts-1))
;
;

c_r_to_a, itime, time(0)
itime(3:5) = 0
c_a_to_r, itime, basetime
hour = (time/3600.0 mod 24.0) + fix((time-basetime)/(24.0*3600.0))*24.0
localtime = (Lons/15.0 + hour) mod 24.0

c_r_to_a, iTime, time(0)
year = tostr(iTime(0))
juli = chopr('000'+tostr(jday(iTime(0),iTime(1),iTime(2))),3)

loadct,39

setdevice, 'gitm_ON2.ps'
max_value=1.2
min_value=0
loc=where(reform(lat_grad(*)) lt 0.0)

contour,ratio(loc),lons(loc),lats(loc),/irr,levels=findgen(31)/30.*max_value,$
/fill,xstyle=1,ystyle=1,xrange=[0,360],yrange=[-90,90],pos=[0.1,0.1,0.8,0.6],title='O/N!D2!N ratio'
pos=[0.82,0.1,0.86,0.6]
    title=''
    maxmin=[min_value,max_value]
    plotct,254,pos,maxmin,title,/right,color=color

;contour,localtime(loc),lons(loc),lats(loc),/irr,/fill,levels=findgen(31)/30*24.,$
;xstyle=1,ystyle=1,xrange=[0,360],yrange=[-90,90],pos=[0.1,0.1,0.8,0.6]
;pos=[0.82,0.1,0.86,0.6]
;    title=''
;    maxmin=[0,24]
;    plotct,254,pos,maxmin,title,/right,color=color
;
closedevice

;champdir = '/remotehome/ridley/GITM/200310/Data'
;;posfile  = year+'-'+juli+'.dat'
;;massfile = 'champ_acc_'+juli+'.dat'
;posfile  = 'champ_pos_all.dat'
;massfile = 'champ_acc_all.dat'
;
;read_champ, massfile, posfile, ChampTime, ChampPosition, MassDensity, $
;  dir = champdir, iYear = 2003, iDay=jday(2003,10,29)
;
;;;  champdir = '/Volumes/Secondary/ridley/Gitm/Runs/20031029/Data'
;;;  posfile  = year+'-'+juli+'.dat'
;;;  massfile = 'champ_acc_'+juli+'.dat'
;;;  
;;;  read_champ, massfile, posfile, ChampTime, ChampPosition, MassDensity, $
;;;    dir = champdir
;
;ChampPosition(2,*) = 400.0
;
;ChampDensity    = fltarr(nTimes)
;GitmDensity     = fltarr(nTimes)
;GitmUne         = fltarr(nTimes)
;GitmUnn         = fltarr(nTimes)
;GitmDensityHigh = fltarr(nTimes)
;GitmDensityLow  = fltarr(nTimes)
;ChampAltitude   = fltarr(nTimes)
;
;for iTime = 0, nTimes-1 do begin
;
;    dt = abs(time(iTime)-ChampTime)
;    loc = where(dt eq min(dt))
;
;    i = loc(0)
;
;    ChampDensity(iTime)  = MassDensity(i)
;    ChampAltitude(iTime) = ChampPosition(2,i)
;
;    loc = where(GitmAlts gt ChampAltitude(iTime))
;    i = loc(0)
;    x = (ChampAltitude(iTime) - GitmAlts(i-1)) / $
;      (GitmAlts(i) - GitmAlts(i-1))
;    GitmDensity(iTime) = (1.0 - x) * rho(iTime,i-1) + $
;                         (      x) * rho(iTIme,i)
;
;    GitmUne(iTime) = (1.0 - x) * une(iTime,i-1) + $
;                     (      x) * une(iTIme,i)
;    GitmUnn(iTime) = (1.0 - x) * unn(iTime,i-1) + $
;                     (      x) * unn(iTIme,i)
;
;    h = (GitmAlts(i+1) - GitmAlts(i-1))/2.0
;
;    loc = where(GitmAlts gt ChampAltitude(iTime)+h)
;    i = loc(0)
;    x = ((ChampAltitude(iTime)+h) - GitmAlts(i-1)) / $
;      (GitmAlts(i) - GitmAlts(i-1))
;    GitmDensityHigh(iTime) = (1.0 - x) * rho(iTime,i-1) + $
;                             (      x) * rho(iTIme,i)
;
;    loc = where(GitmAlts gt ChampAltitude(iTime)-h)
;    i = loc(0)
;    x = ((ChampAltitude(iTime)-h) - GitmAlts(i-1)) / $
;      (GitmAlts(i) - GitmAlts(i-1))
;    GitmDensityLow(iTime) = (1.0 - x) * rho(iTime,i-1) + $
;                            (      x) * rho(iTIme,i)
;
;endfor
;
;
;nOrbits = 0
;
;day   = where(localtime gt 6.0 and localtime lt 18.0,nPtsDay)
;night = where(localtime lt 6.0 or  localtime gt 18.0,nPtsNight)
;
;for i = 1, nPtsDay-1 do begin
;
;    if (day(i)-day(i-1) gt 1) then begin
;        if (nOrbits eq 0) then begin
;            DayOrbitStart = day(i)
;            DayOrbitEnd   = day(i-1)
;        endif else begin
;            DayOrbitStart = [DayOrbitStart,day(i)]
;            DayOrbitEnd   = [DayOrbitEnd  ,day(i-1)]
;        endelse
;        nOrbits = nOrbits+1
;    endif
;
;endfor
;
;nY = max(DayOrbitStart - DayOrbitEnd)
;
;xDay = fltarr(nOrbits,nY)
;yDay = fltarr(nOrbits,nY)
;vDay = fltarr(nOrbits,nY)
;
;GitmVariable = GitmDensity/1.0e-12
;GitmVariable = GitmDensity/1.0e-12

;iOrbit = 0
;iY = 0
;iFound = 0
;for i = 1, nPtsDay-1 do begin
;
;    if (day(i)-day(i-1) gt 1) then begin
;        iOrbit = iOrbit+1
;        iY = 0
;        iFound = 1
;    endif else iY = iY + 1
;
;;print,iorbit-1,norbits,iy,ny
;
;    if (iFound) then begin
;        xDay(iOrbit-1, iY) = hour(DayOrbitStart(iOrbit-1))
;        yDay(iOrbit-1, iY) = Lats(Day(i))
;        vDay(iOrbit-1, iY) = oInt(Day(i))/n2Int(Day(i))
;    endif
;
;endfor
;
;for iOrbit = 0, nOrbits-2 do begin
;    l = where(xday(iOrbit,*) eq 0,c)
;    if (c gt 0) then begin
;        for j = 0,c-1 do begin
;            xDay(iOrbit,l(j)) = xDay(iOrbit,l(j)-1)
;            yDay(iOrbit,l(j)) = yDay(iOrbit,l(j)-1)
;            vDay(iOrbit,l(j)) = vDay(iOrbit,l(j)-1)
;        endfor
;    endif
;endfor
;
;
;setdevice, 'gitm_ON2_daynight.ps', 'p', 5
;
;ppp = 2
;space = 0.1
;pos_space, ppp, space, sizes, ny = ppp
;
;get_position, ppp, space, sizes, 0, pos1, /rect
;get_position, ppp, space, sizes, 1, pos2, /rect
;
;pos1(0) = pos1(0)+0.07
;pos2(0) = pos2(0)+0.07
;
;pos1(2) = pos1(2)-0.11
;pos2(2) = pos2(2)-0.11
;
;;makect, 'bgr'
;;makect, 'bwr'
;
;
;ytickv = [-90,-60,-30,0,30,60,90]
;
;c_r_to_a, itime, time(0)
;c_a_to_s, itime, stime
;xtitle = 'Hours After '+ strmid(stime,0,9) 
;
;;mini = -1000.0
;;maxi = 1000.0
;mini = 0
;maxi = 1.20
;nLevels = 31
;nLevels2 = 5
;levels = findgen(nLevels)*(maxi-mini)/(nLevels-1) + mini
;levels2 = findgen(nLevels2)*(maxi-mini)/(nLevels2-1) + mini
;
;contour, vday(0:nOrbits-2,0:nY-3), $
;  xday(0:nOrbits-2,0:nY-3),yday(0:nOrbits-2,0:nY-3), $
;  /fill, nlevels = nLevels, levels = levels, $
;  yrange = [-90,90], ystyle = 1, $
;  ytickv = ytickv, yticks = 7, yminor = 6, $
;  xstyle = 1, xtitle = xtitle, ytitle = 'Geo. Latitude', $
;  pos = pos1
;
;contour, vday(0:nOrbits-2,0:nY-3), $
;  xday(0:nOrbits-2,0:nY-3),yday(0:nOrbits-2,0:nY-3), $
;  nlevels = nLevels2, levels = levels2, $
;  yrange = [-90,90], ystyle = 1, $
;  ytickv = ytickv, yticks = 7, yminor = 6, $
;  xstyle = 1, $
;  pos = pos1, /over
;
;for i=1,10 do oplot, [i,i]*24.0, [-1000,1000], linestyle = 1
;
;
;xyouts, 2.0, 95.0, '(A) Dayside'
;
;ctpos = [pos1(2)+0.01,pos1(1),pos1(2)+0.03,pos1(3)]
;plotct,254,ctpos,mm(levels),'Mass Density',/right
;;plotct,254,ctpos,mm(levels),'Un (East)',/right
;
;;closedevice
;
;nNOrbits = 0
;
;for i = 1, nPtsNight-1 do begin
;
;    if (night(i)-night(i-1) gt 1) then begin
;        if (nNorbits eq 0) then begin
;            NightorbitStart = night(i)
;            NightOrbitEnd   = night(i-1)
;        endif else begin
;            NightOrbitStart = [NightOrbitStart,night(i)]
;            NightOrbitEnd   = [NightOrbitEnd  ,night(i-1)]
;        endelse
;        nNorbits = nNorbits+1
;    endif
;
;endfor
;
;nY = max(NightOrbitStart - NightOrbitEnd)
;
;xNight = fltarr(nNorbits,nY)
;yNight = fltarr(nNorbits,nY)
;vNight = fltarr(nNorbits,nY)
;
;iNorbit = 0
;iY = 0
;iFound = 0
;for i = 1, nPtsNight-1 do begin
;
;    if (night(i)-night(i-1) gt 1) then begin
;        iNorbit = iNorbit+1
;        iY = 0
;        iFound = 1
;    endif else iY = iY + 1
;
;    if (iFound) then begin
;        xNight(iNorbit-1, iY) = hour(NightorbitStart(iNorbit-1))
;        yNight(iNorbit-1, iY) = Lats(Night(i))
;        vNight(iNorbit-1, iY) = oInt(Night(i))/n2Int(Night(i))
;    endif
;
;endfor
;
;for iOrbit = 0, nOrbits-2 do begin
;    l = where(xNight(iOrbit,*) eq 0,c)
;    if (c gt 0) then begin
;        for j = 0,c-1 do begin
;            xNight(iOrbit,l(j)) = xNight(iOrbit,l(j)-1)
;            yNight(iOrbit,l(j)) = yNight(iOrbit,l(j)-1)
;            vNight(iOrbit,l(j)) = vNight(iOrbit,l(j)-1)
;        endfor
;    endif
;endfor
;
;ytickv = [-90,-60,-30,0,30,60,90]
;
;c_r_to_a, itime, time(0)
;c_a_to_s, itime, stime
;xtitle = 'Hours After '+ strmid(stime,0,9) 
;
;contour, vnight(0:nNorbits-2,0:nY-3), $
;  xnight(0:nNorbits-2,0:nY-3),ynight(0:nNorbits-2,0:nY-3), $
;  /fill, nlevels = nLevels, levels = levels, $
;  yrange = [-90,90], ystyle = 1, $
;  ytickv = ytickv, yticks = 7, yminor = 6, $
;  xstyle = 1, xtitle = xtitle, ytitle = 'Geo. Latitude', $
;  pos = pos2, /noerase
;
;contour, vnight(0:nOrbits-2,0:nY-3), $
;  xnight(0:nOrbits-2,0:nY-3),ynight(0:nOrbits-2,0:nY-3), $
;  nlevels = nLevels2, levels = levels2, $
;  yrange = [-90,90], ystyle = 1, $
;  ytickv = ytickv, yticks = 7, yminor = 6, $
;  xstyle = 1, $
;  pos = pos2, /over
;
;for i=1,10 do oplot, [i,i]*24.0, [-1000,1000], linestyle = 1
;
;xyouts, 2.0, 95.0, '(B) Night Side'
;
;ctpos = [pos2(2)+0.01,pos2(1),pos2(2)+0.03,pos2(3)]
;plotct,254,ctpos,mm(levels),'Mass Density',/right
;
;
;closedevice
;
end
