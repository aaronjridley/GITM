
GetNewData = 1

if (n_elements(data) gt 0) then begin
    filelist_new = findfile("*/cham_*.bin")
    nfiles_new = n_elements(filelist_new)
    if (nfiles_new eq nfiles) then default = 'n' else default='y'
    GetNewData = mklower(strmid(ask('whether to reread data',default),0,1))
    if (GetNewData eq 'n') then GetNewData = 0 else GetNewData = 1
endif else filelist_new = findfile("*/cham_*.bin")

if (GetNewData) then begin
    thermo_readsat, filelist_new,data, time, nTimes, Vars, nAlts, nSats, nFiles
    nfiles = n_elements(filelist_new)
endif

nPts = nTimes

GitmAlts = reform(data(0,0,2,0:nalts-1))/1000.0
Lons = reform(data(0,0:nPts-1,0,0)) * 180.0 / !pi
Lats = reform(data(0,0:nPts-1,1,0)) * 180.0 / !pi

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

c_r_to_a, iTime, time(0)
year = tostr(iTime(0))
juli = chopr('000'+tostr(jday(iTime(0),iTime(1),iTime(2))),3)

champdir = '/remotehome/ridley/GITM/200310/Data'
;posfile  = year+'-'+juli+'.dat'
;massfile = 'champ_acc_'+juli+'.dat'
posfile  = 'champ_pos_all.dat'
massfile = 'champ_acc_all.dat'

read_champ, massfile, posfile, ChampTime, ChampPosition, MassDensity, $
  dir = champdir, iYear = 2003, iDay=jday(2003,10,29)

ChampPosition(2,*) = 400.0

ChampDensity    = fltarr(nTimes)
GitmDensity     = fltarr(nTimes)
GitmDensityHigh = fltarr(nTimes)
GitmDensityLow  = fltarr(nTimes)
ChampAltitude   = fltarr(nTimes)

for iTime = 0, nTimes-1 do begin

    dt = abs(time(iTime)-ChampTime)
    loc = where(dt eq min(dt))

    i = loc(0)

    ChampDensity(iTime)  = MassDensity(i)
    ChampAltitude(iTime) = ChampPosition(2,i)

    loc = where(GitmAlts gt ChampAltitude(iTime))
    i = loc(0)
    x = (ChampAltitude(iTime) - GitmAlts(i-1)) / $
      (GitmAlts(i) - GitmAlts(i-1))
    GitmDensity(iTime) = (1.0 - x) * rho(iTime,i-1) + $
                         (      x) * rho(iTIme,i)

    h = (GitmAlts(i+1) - GitmAlts(i-1))/2.0

    loc = where(GitmAlts gt ChampAltitude(iTime)+h)
    i = loc(0)
    x = ((ChampAltitude(iTime)+h) - GitmAlts(i-1)) / $
      (GitmAlts(i) - GitmAlts(i-1))
    GitmDensityHigh(iTime) = (1.0 - x) * rho(iTime,i-1) + $
                             (      x) * rho(iTIme,i)

    loc = where(GitmAlts gt ChampAltitude(iTime)-h)
    i = loc(0)
    x = ((ChampAltitude(iTime)-h) - GitmAlts(i-1)) / $
      (GitmAlts(i) - GitmAlts(i-1))
    GitmDensityLow(iTime) = (1.0 - x) * rho(iTime,i-1) + $
                            (      x) * rho(iTIme,i)

endfor

print, h

stime = time(0)
etime = max(time)
time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn

c_r_to_a, iTime, mean(Time)
c_a_to_ymd, iTime, ymd

setdevice, 'gitm_champ400_'+ymd+'.ps', 'p', 5

makect, 'mid'

plot, Time-sTime, champdensity/1.0e-12, $
  ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)',   $
;  xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
  xtickname = strarr(10)+' ', xtickv = xtickv, $
  xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
  thick = 3, pos = [0.12, 0.49, 0.84, 0.95], $
  yrange = [0.0,30.0]

oplot, Time-sTime, gitmdensity/1.0e-12, linestyle = 2, color = 250, thick = 3
;oplot, Time-sTime, gitmdensityHigh, linestyle = 1, color = 50, thick = 1
;oplot, Time-sTime, gitmdensityLow, linestyle = 1, color = 50, thick = 1

t1 = (etr-btr)*0.05
t2 = (etr-btr)*0.10
t3 = (etr-btr)*0.11

oplot, [t1,t2], [1.0,1.0], thick = 3, linestyle = 2, color = 250
xyouts, t3, 1.0, 'GITM', color = 250

oplot, [t1,t2], [2.0,2.0], thick = 3
xyouts, t3, 2.0, 'CHAMP'

plot, Time-sTime, 100*(champdensity-gitmdensity)/champdensity, $
  ytitle = 'Percent Difference',   $
  xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
  xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
  thick = 3, pos = [0.12, 0.05, 0.84, 0.48], $
  yrange = [-100.0,100.0], /noerase

oplot, mm(time-stime), [0,0], linestyle = 2

for i=1,10 do oplot, [i,i]*24.0*3600, [-1000,1000], linestyle = 1

rms = 100.0*sqrt(mean((gitmdensity-champdensity)^2))/sqrt(mean(champdensity^2))

rmss = "Normalized RMS : "+string(rms,format="(f4.1)")+"%"
xyouts, 0.86, (0.05+0.48)/2, rmss, /norm, orient=270, align = 0.5


closedevice

setdevice, 'gitm_champ400s_'+ymd+'.ps', 'p', 5

makect, 'mid'

dt = 3600.0*1.5

smooth_1d, time, champdensity, dt, outdensity

t = (Time-sTime)/3600.0
champ = outdensity/1.0e-12

plot, Time-sTime, outdensity/1.0e-12, $
  ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)',   $
  xtickname = strarr(10)+' ', xtickv = xtickv, $
  xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
  thick = 3, pos = [0.12, 0.49, 0.84, 0.95], $
  yrange = [0.0,30.0]

for i=1,10 do oplot, [i,i]*24.0*3600, [-1000,1000], linestyle = 1

smooth_1d, time, gitmdensity, dt, outdensity

gitm = outdensity/1.0e-12

oplot, Time-sTime, outdensity/1.0e-12, linestyle = 2, color = 250, thick = 3

smooth_1d, time, gitmdensityHigh, dt, outdensity
oplot, Time-sTime, outdensity/1.0e-12, linestyle = 1, color = 50, thick = 3

smooth_1d, time, gitmdensityLow, dt, outdensity
oplot, Time-sTime, outdensity/1.0e-12, linestyle = 1, color = 50, thick = 3

t1 = (etr-btr)*0.05
t2 = (etr-btr)*0.10
t3 = (etr-btr)*0.11

smooth = ' ('+tostr(dt/60.0)+' min. smooth)'

oplot, [t1,t2], [4.0,4.0], thick = 3
xyouts, t3, 4.0, 'CHAMP'+smooth

oplot, [t1,t2], [2.5,2.5], thick = 3, linestyle = 2, color = 250
xyouts, t3, 2.5, 'GITM'+smooth, color = 250

oplot, [t1,t2], [1.0,1.0], thick = 1, linestyle = 1, color = 50
xyouts, t3, 1.0, 'GITM (+/- 1 grid height)'+smooth, color = 50

t = (Time-sTime)/3600.0

smooth_1d, time, champdensity, dt, outdensity
champ = outdensity/1.0e-12

smooth_1d, time, gitmdensity, dt, outdensity
gitm = outdensity/1.0e-12

plot, Time-sTime, 100*(champ-gitm)/champ, $
  ytitle = 'Percent Difference',   $
  xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
  xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
  thick = 3, pos = [0.12, 0.05, 0.84, 0.48], $
  yrange = [-100.0,100.0], /noerase

oplot, mm(time-stime), [0,0], linestyle = 2

for i=1,10 do oplot, [i,i]*24.0*3600, [-1000,1000], linestyle = 1

rms = 100.0*sqrt(mean((gitm-champ)^2))/sqrt(mean(champ^2))

rmss = "Normalized RMS : "+string(rms,format="(f4.1)")+"%"
xyouts, 0.86, (0.05+0.48)/2, rmss, /norm, orient=270, align = 0.5

closedevice

end
