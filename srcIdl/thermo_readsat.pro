
pro thermo_readsat, filelist, $
                    data, time, nTimes, vars_t, nAlts_t, nSats, nFiles

    nfiles = n_elements(filelist)

    if (strlen(filelist(0)) gt 0) then begin
        f = filelist(0)
        iSatPosStart = 5
        if (strpos(filelist(0),'ALL') gt 0) then iSatPosStart = 1
        while (strpos(f,'/') gt -1) do begin
            sp = strpos(f,'/')+1
            iSatPosStart = iSatPosStart + sp
            f = strmid(f,sp,strlen(f))
        endwhile
    endif

    iSatPosLength = 3
    if (strpos(filelist(0),'ALL') gt 0) then iSatPosLength = 4
    if (strpos(filelist(0),'ALL') gt 0 and $
        strpos(filelist(0),'_b0') gt 0) then iSatPosLength = 4+6

    iCurSat = 0
    iFile = 0

    nSats = 0
    sats    = strarr(100)
    nSwaths = intarr(100)
    nPts    = intarr(100,100)
    nPts(0,0) = 1
    
    sats(0) = strmid(filelist(0),0,4)
    
    n = 1
    for i=1,nFiles-1 do begin
        cFile = filelist(i)
        if (strpos(cFile, sats(nSats)) ne 0) then begin
            n = n + 1
            nSats = nSats + 1
            sats(nSats) = strmid(filelist(i),0,4)
        endif else begin
;        n = fix(strmid(cFile,5,3))
            if (n gt nSwaths(nSats)) then nSwaths(nSats) = n
            nPts(nSats,n-1) = nPts(nSats,n-1) + 1
        endelse
    endfor
    nSats = nSats + 1
    
    if (nSats gt 1) then begin
        for i=0,nSats-1 do print, tostr(i)+'. '+sats(i)
        if (n_elements(isat) eq 0) then isat = 0
        isat = fix(ask('satellite to plot',tostr(isat)))
    endif else isat = 0

    if (nSats gt 1) then begin
        nPtsSw = nSwaths(iSat)
    endif else nPtsSw = 1

    nTimes = max(npts(iSat,*))

    it = 0

    for i=0,nFiles-1 do begin

        cFile = filelist(i)

        if (strpos(cFile, sats(isat)) eq 0) then begin

            read_thermosphere_file, cFile,nvars_t, nalts_t, nlats_t, nlons_t, $
              vars_t, data_t, nBLKlat_t, nBLKlon_t, nBLK_t, itime, version

            if (it eq 0) then begin
                data = fltarr(nPtsSw, nTimes, nvars_t, nalts_t)
                time = dblarr(nTimes)
            endif
            
;            if (strpos(cFile,'1DALL') lt 0) then begin
;                if (n ne fix(strmid(cFile,iSatPosStart,iSatPosLength))-1 and $
;                    nSats gt 1) then it = 0
;                if (nSats gt 1) then $
;                  n = fix(strmid(cFile,iSatPosStart,iSatPosLength))-1 else n = 0
;            endif else begin
                it = i
                n = 0
;            endelse

            data(n, it, *, *) = data_t(*,0,0,*)

;            if (n eq 0) then begin
;                iSP = iSatPosStart + iSatPosLength + 2
;                itime = [ $
;                          fix(strmid(cfile,iSP   ,2)), $
;                          fix(strmid(cfile,iSP+ 2,2)), $
;                          fix(strmid(cfile,iSP+ 4,2)), $
;                          fix(strmid(cfile,iSP+ 7,2)), $
;                          fix(strmid(cfile,iSP+ 9,2)), $
;                          fix(strmid(cfile,iSP+11,2))]
                c_a_to_r, itime, rtime
                time(it) = rtime
;            endif

            it = it + 1

        endif

    endfor

end

