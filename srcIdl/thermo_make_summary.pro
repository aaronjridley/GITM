
iRho  = 3
iTemp = 15
iVV   = 18
iAlt_  = 2

if (n_elements(dir) eq 0) then dir = 'data.01'
dir = ask('directory',dir)

filelist = findfile(dir+'/3DALL*.bin')

nFiles = n_elements(filelist)

for iFile = 0, nFiles-1 do begin

    file = filelist(iFile)
    print, "Reading File : ",file
    gitm_read_bin, file, data, time, nVars, Vars, version

    if (iFile eq 0) then begin

        nLons = n_elements(data(0,*,0,0))
        nLats = n_elements(data(0,0,*,0))
        nAlts = n_elements(data(0,0,0,*))

        global_rho  = fltarr(4, nFiles, nAlts-4)
        global_temp = fltarr(4, nFiles, nAlts-4)
        global_vv   = fltarr(4, nFiles, nAlts-4)
        global_alt  = fltarr(4, nFiles, nAlts-4)

        hem_temp    = fltarr(4, 2, nFiles, nAlts-4)
        hem_rho     = fltarr(4, 2, nFiles, nAlts-4)
        hem_vv      = fltarr(4, 2, nFiles, nAlts-4)
        hem_alt     = fltarr(4, 2, nFiles, nAlts-4)

        time_all = dblarr(nFiles)

    endif

    time_all(iFile) = time

    for iAlt = 2, nAlts-3 do begin

        global_rho(0, iFile, iAlt-2) = $
          mean(data(iRho, 2:nLons-3, 2:nLats-3, iAlt))
        global_temp(0, iFile, iAlt-2) = $
          mean(data(iTemp, 2:nLons-3, 2:nLats-3, iAlt))
        global_vv(0, iFile, iAlt-2) = $
          mean(data(iVV, 2:nLons-3, 2:nLats-3, iAlt))
        global_alt(0, iFile, iAlt-2) = $
          mean(data(iAlt_, 2:nLons-3, 2:nLats-3, iAlt))

        global_rho(1, iFile, iAlt-2) = $
          min(data(iRho, 2:nLons-3, 2:nLats-3, iAlt))
        global_temp(1, iFile, iAlt-2) = $
          min(data(iTemp, 2:nLons-3, 2:nLats-3, iAlt))
        global_vv(1, iFile, iAlt-2) = $
          min(data(iVV, 2:nLons-3, 2:nLats-3, iAlt))
        global_alt(1, iFile, iAlt-2) = $
          min(data(iAlt_, 2:nLons-3, 2:nLats-3, iAlt))

        global_rho(2, iFile, iAlt-2) = $
          max(data(iRho, 2:nLons-3, 2:nLats-3, iAlt))
        global_temp(2, iFile, iAlt-2) = $
          max(data(iTemp, 2:nLons-3, 2:nLats-3, iAlt))
        global_vv(2, iFile, iAlt-2) = $
          max(data(iVV, 2:nLons-3, 2:nLats-3, iAlt))
        global_alt(2, iFile, iAlt-2) = $
          max(data(iAlt_, 2:nLons-3, 2:nLats-3, iAlt))

        global_rho(3, iFile, iAlt-2) = $
          stddev(data(iRho, 2:nLons-3, 2:nLats-3, iAlt))
        global_temp(3, iFile, iAlt-2) = $
          stddev(data(iTemp, 2:nLons-3, 2:nLats-3, iAlt))
        global_vv(3, iFile, iAlt-2) = $
          stddev(data(iVV, 2:nLons-3, 2:nLats-3, iAlt))
        global_alt(3, iFile, iAlt-2) = $
          stddev(data(iAlt_, 2:nLons-3, 2:nLats-3, iAlt))

    endfor

    iLs = 2
    iLe = (nLons-4)/4+2-1

    ; South -----------------
    for iAlt = 2, nAlts-3 do begin
        hem_rho(0,0,iFile, iAlt-2) = $
          mean(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(0,0,iFile, iAlt-2) = $
          mean(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(0,0,iFile, iAlt-2) = $
          mean(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(0,0,iFile, iAlt-2) = $
          mean(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(1,0,iFile, iAlt-2) = $
          min(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(1,0,iFile, iAlt-2) = $
          min(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(1,0,iFile, iAlt-2) = $
          min(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(1,0,iFile, iAlt-2) = $
          min(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(2,0,iFile, iAlt-2) = $
          max(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(2,0,iFile, iAlt-2) = $
          max(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(2,0,iFile, iAlt-2) = $
          max(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(2,0,iFile, iAlt-2) = $
          max(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(3,0,iFile, iAlt-2) = $
          stddev(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(3,0,iFile, iAlt-2) = $
          stddev(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(3,0,iFile, iAlt-2) = $
          stddev(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(3,0,iFile, iAlt-2) = $
          stddev(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))
    endfor

    iLs = 3*(nLons-4)/4+2
    iLe = nLons-3

    ; North -----------------
    for iAlt = 2, nAlts-3 do begin

        hem_rho(0,1,iFile, iAlt-2) = $
          mean(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(0,1,iFile, iAlt-2) = $
          mean(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(0,1,iFile, iAlt-2) = $
          mean(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(0,1,iFile, iAlt-2) = $
          mean(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(1,1,iFile, iAlt-2) = $
          min(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(1,1,iFile, iAlt-2) = $
          min(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(1,1,iFile, iAlt-2) = $
          min(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(1,1,iFile, iAlt-2) = $
          min(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(2,1,iFile, iAlt-2) = $
          max(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(2,1,iFile, iAlt-2) = $
          max(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(2,1,iFile, iAlt-2) = $
          max(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(2,1,iFile, iAlt-2) = $
          max(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(3,1,iFile, iAlt-2) = $
          stddev(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(3,1,iFile, iAlt-2) = $
          stddev(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(3,1,iFile, iAlt-2) = $
          stddev(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(3,1,iFile, iAlt-2) = $
          stddev(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

    endfor

endfor

save, file = dir+"/summary.save", hem_rho, hem_temp, hem_vv, hem_alt, $
  global_rho, global_temp, global_vv, global_alt, time_all

end

