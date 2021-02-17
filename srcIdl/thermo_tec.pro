
filelist = findfile('3D*.bin')

north = 1

nfiles = n_elements(filelist)

ff = 0
lf = nfiles-1

for ifile = ff, lf do begin

    file = filelist(iFile)

    print, 'Reading file : ',file

    read_thermosphere_file, file, nvars, nalts, nlats, nlons,vars,data, $
      nBLKlat, nBLKlon, nBLK, iTime, Version

    p = strpos(file, '3DALL')
    if (p eq 0) then p = strlen(file)
    tecfile = strmid(file,0,p)+'tec'

    setdevice, tecfile+'.ps','p',5

    e   = reform(data(32,*,*,*))
    Alt = reform(data(2,*,*,*))
    lats = reform(data(1,*,*,0)) / !dtor
    lons = reform(data(0,*,*,0)) / !dtor
    dalt = alt(*,*,1:nAlts-1) - alt(*,*,0:nAlts-2)
    dalte = dalt * e(*,*,0:nAlts-2)

    tec = fltarr(nLons, nLats)

    for i=2,nAlts-3 do tec = tec + dalte(*,*,i)/ 1.0e16

    print, mm(tec)

    nl = 31
    levels = findgen(nl) * 75.0 / (nl-1)

    ppp = 2
    space = 0.05
    pos_space, ppp, space, sizes
    get_position, ppp, space, sizes, 0, pos

    ctpos = pos
    ctpos(0) = pos(2)+0.01
    ctpos(2) = ctpos(0)+0.03

    !p.position = pos

    utime = itime(3)*3600.0 + $
      itime(4)*60.0 + $
      itime(5)
    utime = utime(0)
    
    p0lon = utime/3600.0 * 360.0 / 24.0
    makect,'all'

    if (north) then map_set, 30.0, 180.0-p0lon, /orthographic, /cont $
    else map_set, 0.0, 180.0-p0lon, /orthographic, /cont

    xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.01, file, /norm, align=0.5
    !p.position = -1

    loc = where(tec gt levels(nl-2),count)
    if (count gt 0) then tec(loc) = levels(nl-2)

    newrat = tec(1:nLons-2,1:nLats-2)

    newlat = lats(1:nLons-2,1:nLats-2)
    newlon = lons(1:nLons-2,1:nLats-2)
    nLons  = nLons-2
    nLats  = nLats-2

    newrat(0,*)       = (newrat(1,*)+newrat(nLons-2,*))/2.0
    newrat(nLons-1,*) = (newrat(1,*)+newrat(nLons-2,*))/2.0
    newrat(*,0)       = mean(newrat(*,1))
    newrat(*,nLats-1) = mean(newrat(*,nLats-2))

    newlon(0,*)       = 0.0
    newlon(nLons-1,*) = 360.0
    newlat(*,0) = -90.0
    newlat(*,nLats-1) =  90.0

    contour, newrat, newlon, newlat, $
      /follow, nlevels = nl, /cell_fill, /over, $
      levels = levels, title = file

;;    contour, tec, lons, lats, $
;;      /follow, nlevels = nl, /cell_fill, /over, $
;;      levels = levels, title = file

        map_continents
        map_grid, lats = findgen(19)*10-90, glinethick=3

        plotct, 255, ctpos, mm(levels), 'TECU', /right

    get_position, ppp, space, sizes, 1, pos

    ctpos = pos
    ctpos(0) = pos(2)+0.01
    ctpos(2) = ctpos(0)+0.03

    !p.position = pos

    p0lon = utime/3600.0 * 360.0 / 24.0
    makect,'all'

    if (north) then map_set, 30.0, -p0lon, /orthographic, /cont, /noerase $
    else map_set, 0.0, 180.0-p0lon, /orthographic, /cont, /noerase

    contour, newrat, newlon, newlat, $
      /follow, nlevels = nl, /cell_fill, /over, $
      levels = levels, title = file

;;    contour, tec, lons, lats, $
;;      /follow, nlevels = nl, /cell_fill, /over, $
;;      levels = levels, title = file

        map_continents
        map_grid, lats = findgen(19)*10-90, glinethick=3

        plotct, 255, ctpos, mm(levels), 'TECU', /right

endfor


end
