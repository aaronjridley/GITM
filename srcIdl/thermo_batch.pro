dir = "/grid/swmf/Oct26/SixCompFixed3/UA"

filelist = findfile("-t "+dir+"/*.save")
if (strlen(filelist(0)) eq 0) then filelist = findfile("-t *.bin")

filelist = ask('filename to plot',filelist(0))

filelist = findfile(filelist)

nfiles = n_elements(filelist)

for iFile = 0, nFiles-1 do begin

    filename = filelist(iFile)

    print, 'Reading file ',filename

    read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
      vars, data, rb, cb, bl_cnt, iTime, Version

    alt = reform(data(2,*,*,*)) / 1000.0
    lat = reform(data(1,*,*,*)) / !dtor
    lon = reform(data(0,*,*,*)) / !dtor

    if (iFile eq 0) then begin

        for i=0,nvars-1 do print, tostr(i)+'. '+vars(i)
        sel = fix(ask('which var to plot','9'))

        plotlog = ask('whether you want log or not (y/n)','n')
        if (strpos(plotlog,'y') eq 0) then plotlog = 1 else plotlog = 0

        psfile = filename+'.ps'
        psfile = ask('ps file name',psfile)

        print, '1. Constant Altitude Plot'
        print, '2. Constant Longitude Plot'
        print, '3. Constant Latitude Plot'
        slice = fix(ask('type of plot to make','1'))

        cnt1 = 0
        cnt2 = 0
        cnt3 = 0

;cnt1 is a lat/lon plot
        if (slice eq 1) then cnt1 = 1

;cnt1 is a lat/alt plot
        if (slice eq 2) then cnt3 = 1

;cnt1 is a lon/alt plot
        if (slice eq 3) then cnt2 = 1

        if (slice eq 1) then begin
            for i=0,nalts-1 do print, tostr(i)+'. '+string(alt(2,2,i))
            selset = fix(ask('which altitude to plot','0'))
        endif

        if (slice eq 2) then begin
            for i=0,nlons-1 do print, tostr(i)+'. '+string(lon(i,2,2))
            selset = fix(ask('which longitude to plot','0'))
        endif

        if (slice eq 3) then begin
            for i=0,nlats-1 do print, tostr(i)+'. '+string(lat(2,i,2))
            selset = fix(ask('which latitude to plot','0'))
        endif

        smini = ask('minimum (0.0 for automatic)','0.0')
        smaxi = ask('maximum (0.0 for automatic)','0.0')

        plotvector = ask('whether you want vectors or not (y/n)','y')
        if strpos(plotvector,'y') eq 0 then plotvector=1 else plotvector = 0

        if (plotvector) then begin
            print,'-1  : automatic selection'
            factors = [1.0, 5.0, 10.0, 20.0, 25.0, $
                       50.0, 75.0, 100.0, 150.0, 200.0,300.0]
            nfacs = n_elements(factors)
            for i=0,nfacs-1 do print, tostr(i)+'. '+string(factors(i)*10.0)
            vector_factor = fix(ask('velocity factor','-1'))
        endif else vector_factor = 0

; cursor position variables, which don't matter at this point
        cursor_x = 0.0
        cursor_y = 0.0
        strx = '0.0'
        stry = '0.0'

; yes is whether ghostcells are plotted or not:
        yes = 0
        no  = 1

; yeslog is whether variable should be logged or not:
        if (plotlog) then begin 
            yeslog = 1
            nolog  = 0
        endif else begin
            yeslog = 0
            nolog = 1
        endelse

; yeswrite_cnt is whether we have to output to a ps file or not.
        yeswrite_cnt = 1

; polar is variable to say whether we have polar plots or not
        polar = 1

; npolar is whether we are doing the northern or southern hemisphere
        npolar = 1

; MinLat is for polar plots:
        MinLat = 50.0

; showgridyes says whether to plot the grid or not.
        showgridyes = 0

;plotvectoryes says whether to plot vectors or not
        plotvectoryes = plotvector

; number of points to skip when plotting vectors:
        step = 1

; vi_cnt is whether to plot vectors of Vi
        vi_cnt = 1

; vn_cnt is whether to plot vectors of Vn
        vn_cnt = 1-vi_cnt

        cursor_cnt = 0

        xrange = [0.0,0.0]

        yrange = [0.0,0.0]

    endif

    if (nFiles gt 1) then begin
        p = strpos(psfile,'.ps')
        if (p gt -1) then psfile = strmid(psfile,0,p)
        psfile_final = psfile+'_'+chopr('000'+tostr(iFile),4)+'.ps'
    endif else begin
        psfile_final = psfile
    endelse

    smini_final = smini
    smaxi_final = smaxi

    thermo_plot,cursor_x,cursor_y,strx,stry,step,nvars,sel,nfiles,	$
		     cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
		     1-yeslog,nalts,nlats,nlons,yeswrite_cnt,$
		     polar,npolar,MinLat,showgridyes,	  $
		     plotvectoryes,vi_cnt,vn_cnt,vector_factor,	  $
		     cursor_cnt,data,alt,lat,lon,	  $
		     xrange,yrange,selset,smini_final,smaxi_final,	  $
		     filename,vars, psfile_final, 0, 'mid', itime

endfor


end
