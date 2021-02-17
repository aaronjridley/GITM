
filelist = findfile("t98050[2-5]*.save")

nfiles = n_elements(filelist)

tn_  = 14
alt_ = 2
pot_ = 28

temp = fltarr(nfiles, 22)
cpcp = fltarr(nfiles)

for ifile = 0, nfiles-1 do begin

  file = filelist(ifile)

  print, "Reading File : ", file

  read_thermosphere_file, file, nvars, nalts, nlats, nlons,vars,data

  nl = nlons/2

  for i=0,21 do temp(ifile, i) = max(data(tn_,*,nl:nlons-1,i))

  if (n_elements(data(*,0,0,0)) gt 28) then begin
      cpcp(ifile) = max(data(pot_,*,nl:nlons-1,*)) - $
        min(data(pot_,*,nl:nlons-1,*))
  endif else cpcp(ifile) = 0.0

endfor

cpcp = cpcp/1000.0

alts = reform(data(alt_,0,0,*))

time = findgen(nfiles)

setdevice, "temp.ps", "p", 4, 0.95

plotdumb

ppp = 3
space = 0.01
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0) + 0.075

plot, time, temp(*,7), xstyle = 1, ystyle = 1, $
  ytitle = "Temp. @ 125 km (K)", $
  pos = pos, /noerase, xtickname = strarr(20)+' ', $
  yrange = [300,550], charsize = 1.2, thick = 3, xrange = [0,96]

for t=1,10 do begin
    oplot, 24.0*[t,t], [-1.0e6, 1.0e6], linestyle = 1
endfor

get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0) + 0.075

plot, time, temp(*,15), xstyle = 1, ystyle = 1, $
  ytitle = "Temp. @ 290 km (K)", $
  pos = pos, /noerase, xtickname = strarr(20)+' ', $
  yrange = [1200,3200], charsize = 1.2, thick = 3, xrange = [0,96]

for t=1,10 do begin
    oplot, 24.0*[t,t], [-1.0e6, 1.0e6], linestyle = 1
endfor

get_position, ppp, space, sizes, 2, pos, /rect
pos(0) = pos(0) + 0.075

plot, time, cpcp(*), xstyle = 1, ystyle = 1, $
  ytitle = "CPCP (kV)", xtitle = "Hours since May 2, 1998 00 UT", $
  pos = pos, /noerase, $
  yrange = [0,375], charsize = 1.2, thick = 3, xrange = [0,96]

for t=1,10 do begin
    oplot, 24.0*[t,t], [-1.0e6, 1.0e6], linestyle = 1
endfor

closedevice

end
