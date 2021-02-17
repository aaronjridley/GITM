
filelist = findfile("b0001_*.*ALL")
if (strlen(filelist(0)) eq 0) then filelist = findfile("b0001_*")

nfiles = n_elements(filelist)

if (nFiles gt 1 or strlen(filelist(0)) gt 1) then begin

    for ifile = 0, nfiles - 1 do begin

        nBLKlat = 0
        nBLKlon = 0
        nBLK = 0

        file = filelist(ifile)
        read_thermosphere_file, file, nvars, nalts, nlats, nlons,vars,data, $
          nBLKlat, nBLKlon, nBLK, iTime, Version

        ofile = strmid(file,6,strlen(file)-6)+".save"
        print, ofile
        save, file=ofile, nvars, nalts, nlats, nlons, vars, data, $
          iTime, Version

    endfor

endif


filelist = findfile("*.bin")

nfiles = n_elements(filelist)

if (nFiles gt 1 or strlen(filelist(0)) gt 1) then begin

    for ifile = 0, nfiles - 1 do begin

        nBLKlat = 0
        nBLKlon = 0
        nBLK = 0

        file = filelist(ifile)
        read_thermosphere_file, file, nvars, nalts, nlats, nlons,vars,data, $
          nBLKlat, nBLKlon, nBLK, iTime, Version

        ofile = strmid(file,0,strlen(file)-4)+".save"
        print, ofile
        save, file=ofile, nvars, nalts, nlats, nlons, vars, data, $
          iTime, Version

    endfor

endif


filelist_long = findfile('*0000.dat')

if (n_elements(filelist_long) gt 1) then begin

    done = 0
    iFile = 1
    iFileCurrent = 0
    nFiles = n_elements(filelist_long)

    while (not done) do begin

        Old = strmid(filelist_long(iFileCurrent),0,4)
        New = strmid(filelist_long(iFile),0,4)

        if ((strpos(Old,New) eq -1) or (iFile eq 1)) then begin

            print, Old, " -> ", New
            filelist_new = findfile(New+"*.dat")
            thermo_readsat, filelist_new, data, time, nTimes, Vars, $
              nAlts, nSats, Files
            save, file=New+".save",data, time, nTimes, Vars, nAlts, nSats
            iFileCurrent = iFile

        endif

        iFile = iFile + 1
        if (iFile gt nFiles-1) then Done = 1

    endwhile

endif

end
