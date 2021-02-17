
file = findfile("*.save")

tdfile = $
  '/Users/ridley/Data/GITM/200204/Outputs/17/3D/t020417_003000.3DALL.save'
tdfile = ask('3d file',tdfile)
restore, tdfile

tddata = data

p = strpos(tdfile,"DALL")-15
itime = [ $
          fix(strmid(tdfile,p,2)), $
          fix(strmid(tdfile,p+2,2)), $
          fix(strmid(tdfile,p+4,2)), $
          fix(strmid(tdfile,p+7,2)), $
          fix(strmid(tdfile,p+9,2)), $
          fix(strmid(tdfile,p+11,2))]

print, itime
c_a_to_r, itime, tdtime
c_a_to_s, itime, tdstime

restore, file(0)

datafile = '/Users/ridley/Data/GITM/200204/Inputs/GUVI_im_disk_v008r01_2002107_REV01929.L1B'
datafile = ask('data file to plot',datafile)

if (strpos(datafile,'GUVI') gt -1) then begin
    read_guvi, datafile, sattime, satdata, satlatitude, satlongitude
    filters = ['lyman alpha (121.6nm)', $
               'oxygen 130.4 nm', $
               'oxygen 135.6 nm', $
               'N2 LBH band about 140 to 150 nm', $
               'N2 LBH band about 165 to 180 nm']

    nfilters = n_elements(filters)
    display, filters
    ifilter = fix(ask('filter to examine','2'))

    loc = where(sattime ge min(time) and sattime le max(time), count)

    if (count le 0) then begin
        print, "No satellite data found during this period!!"
        stop
    endif

    satdata = reform(satdata(ifilter,0,*,loc))
    satlat  = reform(satlatitude(0,*,loc))
    satlon  = reform(satlongitude(0,*,loc))
    
    ntsat = count
    npsat = n_elements(satdata(*,0))

    satvalue  = fltarr(nTsat,nPsat)
    sattimes  = dblarr(nTsat,nPsat)
    satpos    = fltarr(nTsat,nPsat)
    for i=0,nPsat-1 do satvalue(*,i) = satdata(i,*)

    for i=0,nTsat-1 do sattimes(i,*) = sattime(loc(i)) - time(0)
    for i=0,nTsat-1 do satpos(i,*)   = findgen(nPsat)

endif

if (nSwaths eq 1) then begin

    nPts = nTimes

    Alts = reform(data(0,0:nPts-1,2,2:nAlts_t-3))/1000.0
    Lons = reform(data(0,0:nPts-1,0,0)) * 180.0 / !pi
    Lats = reform(data(0,0:nPts-1,1,0)) * 180.0 / !pi

    time2d = dblarr(nPts,nAlts_t-4)
    for i=0,nPts-1 do time2d(i,*) = time(i)- time(0)

    display, vars
    iVar = 3
    iVar = fix(ask('variable to plot',tostr(iVar)))

    value = reform(data(0,0:nPts-1,iVar,2:nAlts_t-3))

    if (min(value) gt 0) then begin
        an = ask('whether you would like variable to be alog10','y')
        if (strpos(mklower(an),'y') eq 0) then begin
            value = alog10(value)
            title = 'alog10('+vars(ivar)+')'
        endif else title = vars(ivar)
    endif else title = vars(ivar)

    setdevice, 'test.ps', 'p', 5, 0.95

    makect, 'mid'

    ppp = 8
    space = 0.01
    pos_space, ppp, space, sizes, ny = ppp
    
    get_position, ppp, space, sizes, 4, pos1, /rect
    get_position, ppp, space, sizes, 7, pos2, /rect
    pos = [pos1(0)+0.05,pos2(1), pos1(2)-0.07,pos1(3)]

    mini = min(value)
    maxi = max(value)
    range = (maxi-mini)
    mini = mini - 0.1*range
    maxi = maxi + 0.1*range
    levels = findgen(31) * (maxi-mini) / 30 + mini

    stime = time(0)
    etime = max(time)
    time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn
    xtitle = strmid(xtitle,0,12)

    contour, value, time2d, Alts, /follow, /fill, $
      nlevels = 30, pos = pos, levels = levels, $
      yrange = [0,max(alts)], ystyle = 1, ytitle = 'Altitude (km)', $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1

    ctpos = pos
    ctpos(0) = pos(2)+0.01
    ctpos(2) = ctpos(0)+0.03
    maxmin = [mini,maxi]
    plotct, 255, ctpos, maxmin, title, /right

; Put the max and min on the plot

    mini_tmp = min(value)
    maxi_tmp = max(value)

    r = (maxi_tmp - mini_tmp)/50.0

    plots, [0.0,1.0], [mini_tmp, mini_tmp], thick = 5
    plots, [1.0,0.6], [mini_tmp, mini_tmp+r], thick = 2
    plots, [1.0,0.6], [mini_tmp, mini_tmp-r], thick = 2

    plots, [0.0,1.0], [maxi_tmp, maxi_tmp], thick = 5
    plots, [1.0,0.6], [maxi_tmp, maxi_tmp+r], thick = 2
    plots, [1.0,0.6], [maxi_tmp, maxi_tmp-r], thick = 2

    get_position, ppp, space, sizes, 3, pos1, /rect
    pos = [pos1(0)+0.05,pos1(1), pos1(2)-0.07,pos1(3)]

    plot, time-stime, Lons, ytitle = 'Longitude', /noerase, $
      xtickname = strarr(10)+' ', xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, pos = pos, $
      yrange = [0.0,360.0], ystyle = 1, thick = 3
  
    get_position, ppp, space, sizes, 2, pos1, /rect
    pos = [pos1(0)+0.05,pos1(1), pos1(2)-0.07,pos1(3)]

    plot, time-stime, Lats, ytitle = 'Latitude', /noerase, $
      xtickname = strarr(10)+' ', xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, pos = pos, $
      yrange = [-90.0,90.0], ystyle = 1, thick = 3
  
    get_position, ppp, space, sizes, 0, pos1, /rect
    get_position, ppp, space, sizes, 1, pos2, /rect
    pos = [pos1(0)+0.05,pos2(1), pos1(2)-0.07,pos1(3)]

    !p.position = pos

    map_set, /noerase
    map_continents, color = 0

    oplot, lons, lats, psym = 1

    xyouts, lons(0), lats(0), 'Start', charsize = 1.2
    xyouts, lons(nPts-1), lats(nPts-1), 'End', charsize = 1.2

    closedevice

endif else begin

    lons = reform(data(*,*, 0, 0))*180.0/!pi
    lats = reform(data(*,*, 1, 0))*180.0/!pi

    display, vars
    iVar = 3
    iVar = fix(ask('variable to plot',tostr(iVar)))

    iLevel = fix(ask('ilevel to plot','25'))

    alt = data(0,0,2,iLevel)/1000.0

    value_tmp = reform(data(*,*,iVar,iLevel))

    loc = where(lons gt 180,count)
    if count gt 0 then lons(loc) = lons(loc)-360.0

    for i = 0, nPtsSw-1 do begin
        for j = 2, nTimes-1 do begin
            if (abs(lons(i,j)-lons(i,j-1)) gt 10.0 and $
                abs(lons(i,j)+360.0-lons(i,j-1)) gt 10.0) then begin
                lons(i,j) = lons(i,j-1) + (lons(i,j-1)-lons(i,j-2))
                value_tmp(i,j) = value_tmp(i,j-1)
            endif
        endfor
    endfor

    value  = fltarr(nTimes,nPtsSw)
    for i=0,nPtsSw-1 do value(*,i) = value_tmp(i,*)

    time2d = dblarr(nTimes,nPtsSw)
    dummy  = fltarr(nTimes,nPtsSw)
    for i=0,nTimes-1 do time2d(i,*) = time(i)- time(0)
    for i=0,nTimes-1 do dummy(i,*)  = findgen(nPtsSw)

    IsLog = 0
    if (min(value) gt 0) then begin
        an = ask('whether you would like variable to be alog10','y')
        if (strpos(mklower(an),'y') eq 0) then begin
            IsLog = 1
            value_tmp = alog10(value_tmp)
            value = alog10(value)
            title = 'alog10('+vars(ivar)+')'
        endif else title = vars(ivar)
    endif else title = vars(ivar)

    setdevice, 'test.ps', 'p', 5, 0.95

    makect, 'mid'

    ppp = 8
    space = 0.01
    pos_space, ppp, space, sizes, ny = ppp
    
    get_position, ppp, space, sizes, 4, pos1, /rect
    get_position, ppp, space, sizes, 5, pos2, /rect
    pos = [pos1(0)+0.05,pos2(1), pos1(2)-0.07,pos1(3)]

    mini = min(value)
    maxi = max(value)
    range = (maxi-mini)
    mini = mini - 0.1*range
    maxi = maxi + 0.1*range
    levels = findgen(31) * (maxi-mini) / 30 + mini

    stime = time(0)
    etime = max(time)
    time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn
    xtitle = strmid(xtitle,0,12)

    contour, value, time2d, dummy, /follow, /fill, $
      nlevels = 30, pos = pos, levels = levels, $
      ystyle = 1, ytickname=strarr(10)+' ', ytitle = 'Model Results', $
      xtickname = strarr(10)+' ', xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1

    oplot, [tdtime,tdtime]-stime, mm(dummy), linestyle = 1

    ctpos = pos
    ctpos(0) = pos(2)+0.01
    ctpos(2) = ctpos(0)+0.03
    maxmin = [mini,maxi]
    plotct, 255, ctpos, maxmin, title, /right

; Put the max and min on the plot

    mini_tmp = min(value)
    maxi_tmp = max(value)

    r = (maxi_tmp - mini_tmp)/50.0

    plots, [0.0,1.0], [mini_tmp, mini_tmp], thick = 5
    plots, [1.0,0.6], [mini_tmp, mini_tmp+r], thick = 2
    plots, [1.0,0.6], [mini_tmp, mini_tmp-r], thick = 2

    plots, [0.0,1.0], [maxi_tmp, maxi_tmp], thick = 5
    plots, [1.0,0.6], [maxi_tmp, maxi_tmp+r], thick = 2
    plots, [1.0,0.6], [maxi_tmp, maxi_tmp-r], thick = 2

    ;-------------------------------------------------
    ; sat data
    ;-------------------------------------------------

    get_position, ppp, space, sizes, 6, pos1, /rect
    get_position, ppp, space, sizes, 7, pos2, /rect
    pos = [pos1(0)+0.05,pos2(1), pos1(2)-0.07,pos1(3)]

    if (IsLog) then begin
        loc = where(satvalue le 0.0,count)
        if count gt 0 then begin
            loc2 = where(satvalue gt 0.0, count2)
            if (count2 gt 0) then begin
                satvalue(loc) = min(satvalue(loc2))
            endif else satvalue(loc) = 1.0
        endif
        satvalue = alog10(satvalue)
        filters(ifilter) = 'alog10('+filters(ifilter)+')'
    endif

    mini = min(satvalue)
    maxi = max(satvalue)
    range = (maxi-mini)
    mini = mini - 0.1*range
    maxi = maxi + 0.1*range
    levels = findgen(31) * (maxi-mini) / 30 + mini

;    stime = time(0)
;    etime = max(time)
;    time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn
;    xtitle = strmid(xtitle,0,12)

    image = (satvalue-mini) * 255.0 / (maxi-mini)

    tvscl, image, pos(0), pos(1), /norm, $
      xsize = pos(2)-pos(0), ysize=pos(3)-pos(1)

    plot, mm(sattimes), mm(satpos), /nodata, /noerase, $
      ystyle = 1,xrange = mm(time2d),$
      ytickname = strarr(10)+' ', ytitle = 'GUVI Data', $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, pos = pos

;    contour, satvalue, sattimes, satpos, /follow, /fill, $
;      nlevels = 30, pos = pos, levels = levels, $
;      ystyle = 1, /noerase, xrange = mm(time2d), $
;      ytickname = strarr(10)+' ', ytitle = 'Position in Swath', $
;      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
;      xminor = xminor, xticks = xtickn, xstyle = 1

    ctpos = pos
    ctpos(0) = pos(2)+0.01
    ctpos(2) = ctpos(0)+0.03
    maxmin = [mini,maxi]
    plotct, 255, ctpos, maxmin, filters(ifilter), /right

; Put the max and min on the plot

    mini_tmp = min(satvalue)
    maxi_tmp = max(satvalue)

    r = (maxi_tmp - mini_tmp)/50.0

    plots, [0.0,1.0], [mini_tmp, mini_tmp], thick = 5
    plots, [1.0,0.6], [mini_tmp, mini_tmp+r], thick = 2
    plots, [1.0,0.6], [mini_tmp, mini_tmp-r], thick = 2

    plots, [0.0,1.0], [maxi_tmp, maxi_tmp], thick = 5
    plots, [1.0,0.6], [maxi_tmp, maxi_tmp+r], thick = 2
    plots, [1.0,0.6], [maxi_tmp, maxi_tmp-r], thick = 2

    ;-------------------------------------------------
    ; model on Earth
    ;-------------------------------------------------

    mini = min(value_tmp)
    maxi = max(value_tmp)
    range = (maxi-mini)
    mini = mini - 0.1*range
    maxi = maxi + 0.1*range
    levels = findgen(31) * (maxi-mini) / 30 + mini

    get_position, ppp, space, sizes, 2, pos1, /rect
    get_position, ppp, space, sizes, 3, pos2, /rect
    pos = [pos1(0)+0.05,pos2(1), pos1(2)-0.07,pos1(3)]

    !p.position = pos

    map_set, /noerase

    iEnd = 0
    Stop = 0
    while not Stop do begin

        if (max(abs(lons(*,iEnd))) gt 150.0) then begin
            if (max(lons(*,iEnd)) gt 150.0 and $
                min(lons(*,iEnd)) lt 0) then begin
                Stop= 1
                print, "lon stop"
            endif
        endif

        if (max(abs(lats(*,iEnd))) gt 80.0) then begin
            for i = 0, nPtsSw-1 do begin
                if (abs(lats(i,iEnd)) - abs(lats(i,iEnd+1)) lt 0) then begin
                    if (lons(i,iEnd+1) - lons(i,iEnd) gt 20) then begin
                        stop = 1
                        print, "Lat stop", iEnd
                    endif
                endif
            endfor
        endif

        if (not stop) then iEnd = iEnd +1

    endwhile
    iEnd = iEnd-1

    contour, value_tmp(1:nPtsSw-2,0:iEnd-1), $
      lons(1:nPtsSw-2,0:iEnd-1), lats(1:nPtsSw-2,0:iEnd-1), $
      /over, /cell_fill, nlevels = 30, levels = levels

    contour, value_tmp(1:nPtsSw-2,iEnd+2:nTimes-1), $
      lons(1:nPtsSw-2,iEnd+2:nTimes-1), $
      lats(1:nPtsSw-2,iEnd+2:nTimes-1), $
      /over, /cell_fill, nlevels = 30, levels = levels

    map_continents, color = 0

    xyouts, lons(0,0), lats(0,0), 'Start', charsize = 1.2
    xyouts, lons(nPtsSw-1,nTimes-1), lats(nPtsSw-1,nTimes-1), 'End', $
      charsize = 1.2

    ;-------------------------------------------------
    ; 3d model on Earth
    ;-------------------------------------------------

    newqua = reform(tddata(iVar, 1:nLons-2, 1:nLats-2, iLevel))
    newlon = reform(tddata(0, 1:nLons-2, 1:nLats-2, iLevel)) * 180/!pi
    newlat = reform(tddata(1, 1:nLons-2, 1:nLats-2, iLevel)) * 180/!pi

    if (IsLog) then begin
        newqua = alog10(newqua)
    endif

    nLons  = nLons-2
    nLats  = nLats-2

    newqua(0,*)       = (newqua(1,*)+newqua(nLons-2,*))/2.0
    newqua(nLons-1,*) = (newqua(1,*)+newqua(nLons-2,*))/2.0
    newqua(*,0)       = mean(newqua(*,1))
    newqua(*,nLats-1) = mean(newqua(*,nLats-2))

    newlon(0,*)       = 0.0
    newlon(nLons-1,*) = 360.0
    newlat(*,0) = -90.0
    newlat(*,nLats-1) =  90.0

    mini = min(value_tmp)
    maxi = max(value_tmp)
    range = (maxi-mini)
    mini = mini - 0.1*range
    maxi = maxi + 0.1*range
    levels = findgen(31) * (maxi-mini) / 30 + mini

    get_position, ppp, space, sizes, 0, pos1, /rect
    get_position, ppp, space, sizes, 1, pos2, /rect
    pos = [pos1(0)+0.05,pos2(1), pos1(2)-0.07,pos1(3)]

    !p.position = pos

    map_set, /noerase, title = Vars(iVar)+' at '+tostr(fix(alt))+' km '+ $
      strmid(tdstime,0,15)

    contour, newqua, newlon, newlat, $
      /over, /cell_fill, nlevels = 30, levels = levels

    map_continents, color = 0
    
    oplot, lons(1,*), lats(1,*), psym = 3
    oplot, lons(nPtsSw/2,*), lats(nPtsSw/2,*), psym = 3
    oplot, lons(nPtsSw-2,*), lats(nPtsSw-2,*), psym = 3

    closedevice

endelse

!p.position = -1

end
