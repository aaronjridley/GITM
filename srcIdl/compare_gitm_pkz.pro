
;file = '/Users/ridley/Dropbox/Alaska2010Data/Alaska2010Results Analyzed May2011/09-10Jan2010/ResultsPKZ10012010.txt'
;itime = [2010,01,09,0,0,0]

file = '/Users/ridley/Dropbox/Alaska2010Data/Alaska2010Results Analyzed May2011/10-11Jan2010/ResultsPKZ11012010.txt'
itime = [2010,01,10,0,0,0]

;file = '/Users/ridley/Dropbox/Alaska2010Data/Alaska2010Results Analyzed May2011/10-11Feb2010/ResultsPKZ11022010.txt'
;itime = [2010,02,10,0,0,0]

;file = '/Users/ridley/Dropbox/Alaska2010Data/Alaska2010Results Analyzed May2011/11-12Feb2010/ResultsPKZ12022010.txt'
;itime = [2010,02,11,0,0,0]

;file = '/Users/ridley/Dropbox/Alaska2010Data/Alaska2010Results Analyzed May2011/15-16Feb2010/ResultsPKZ16022010.txt'
;itime = [2010,02,15,0,0,0]


c_a_to_r, itime, basetime

psfile = 'pkz_gitm_'+tostr(itime(0),4)+tostr(itime(1),2)+tostr(itime(2),2)+'.ps'

read_meriwether_fpi, file, pkztime, vz, vze, tn, tne, in

nTimes = n_elements(pkztime)

pkztime = pkztime + basetime

print, mm(pkztime)

mintime = min(pkztime)
maxtime = max(pkztime)

c_r_to_a, itimemin, mintime
c_a_to_ymd, itimemin, ymd

ymd = strmid(ymd,2,6)

filelist = findfile('pkr*'+ymd+'*.bin')

print, filelist

stop

gitm_read_bin, filelist, data, gitmtime, nVars, Vars, version

alts = reform(data(0,2,0,0,*)/1000.0)

d = abs(alts-240.0)
iAlt240 = where(d eq min(d))
iAlt240 = iAlt240(0)

d = abs(alts-350.0)
iAlt350 = where(d eq min(d))
iAlt350 = iAlt350(0)

d = abs(alts-110.0)
iAlt110 = where(d eq min(d))
iAlt110 = iAlt110(0)

gitmVz240 = fltarr(nTimes)
gitmVz350 = fltarr(nTimes)

gitmTn240 = fltarr(nTimes)
gitmTn350 = fltarr(nTimes)

gitme = fltarr(nTimes)

iVz = 19
iTn = 15
iE = 33

for i=0,nTimes-1 do begin

   l = where(gitmtime ge pkztime(i),c)
   if (c gt 0) then l = l(0) else l = n_elements(gitmtime)-1

   gitmvz240(i) = data(l,iVz,0,0,iAlt240)
   gitmvz350(i) = data(l,iVz,0,0,iAlt350)

   gitmtn240(i) = data(l,iTn,0,0,iAlt240)
   gitmtn350(i) = data(l,iTn,0,0,iAlt350)

   gitme(i) = data(l,iE,0,0,iAlt110)

endfor


stime = min(pkztime)
etime = max(pkztime)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

setdevice, psfile, 'p', 5

space = 0.04
ppp = 3
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0)+0.05
pos(2) = pos(2)-0.05

maxi = max(abs(vz))
plot, pkztime-stime, vz-mean(vz), linestyle = 1, thick = 4, $
      xtickname = xtickname,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      pos = pos, ystyle = 4, yrange = [-maxi,maxi]

axis, yaxis=0, ytitle = 'PKZ Vz (m/s)'

maxi = max(abs(gitmvz240))
plot, pkztime-stime, gitmvz240-mean(gitmvz240), thick = 4, $
      xtickname = xtickname,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      pos = pos, yrange = [-maxi,maxi], ystyle = 4, /noerase

axis, yaxis=1, ytitle = 'GITM Vz (m/s)'


get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0)+0.05
pos(2) = pos(2)-0.05
plot, pkztime-stime, tn, linestyle = 1, thick = 4, $
      xtickname = xtickname,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      pos = pos, ytitle = 'Temperature (K)', /noerase, $
      yrange = [500,1000]

oplot, pkztime-stime, gitmtn240, thick = 4

get_position, ppp, space, sizes, 2, pos, /rect
pos(0) = pos(0)+0.05
pos(2) = pos(2)-0.05
plot, pkztime-stime, in, linestyle = 1, thick = 4, $
      xtickname = xtickname,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      xtitle = xtitle,   $
      pos = pos, /noerase, ystyle = 4

axis, yaxis=0, ytitle = 'Intensity (R?)'

plot, pkztime-stime, gitme/1e11, linestyle = 0, thick = 4, ystyle = 4, $
      pos = pos, /noerase, $
      xtickname = xtickname,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn

axis, yaxis=1, ytitle = 'Electron Density (x10!U11!N m!U-3!N)'

closedevice

end


