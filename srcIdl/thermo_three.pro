
if (n_elements(filelist) eq 0) then filelist = findfile("-t *.bin")

filelist = ask('filename to plot',filelist(0))
filelist = findfile(filelist)

nfiles = n_elements(filelist)

if (nFiles eq 1 and strlen(filelist(0)) eq 0) then begin
   print, 'can not find file!'
   stop
endif

if (n_elements(psfile) eq 0) then psfile = filelist(0)+'.ps'
psfile = ask('ps file name',psfile)

diffdir = ask('directory to use to make differences (return for no diffs)','')

if (n_elements(percentdiff) eq 0) then percentdiff = 0
if (strlen(diffdir) gt 0) then begin
   percentdiff = fix(ask('whether to plot percent diff (1) or not (0)', $
                         tostr(percentdiff)))
endif

file = filelist(0)
VarsToGet=[0,1,2]
gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                        VarsToGet = VarsToGet

alts = reform(data(2,0,0,*)) / 1000.0
lats = reform(data(1,*,*,0))
lons = reform(data(0,*,*,0))
nAlts = n_elements(alts)
nLons = n_elements(lons(*,0))
nLats = n_elements(lats(0,*))

CalcTEC = 0
display, vars
print,tostr(nVars),'. TEC'
if (n_elements(iVar) eq 0) then iVar = '3' else iVar = tostr(iVar)
iVar = fix(ask('which var to plot',iVar))
if (iVar eq nVars) then begin
   CalcTEC = 1
   iVar = 34
endif

for i=0,nalts-1 do print, tostr(i)+'. '+string(alts(i))
if (n_elements(iAlt) eq 0) then iAlt='0' else iAlt=tostr(iAlt)
iAlt = fix(ask('which altitude to plot',iAlt))

if (n_elements(plotVectors) eq 0) then plotVectors=0
plotVectors = fix(ask('whether you want vectors or not (0/1)',tostr(plotVectors)))

AllValues = fltarr(nFiles, nLons, nLats)
AllTimes = dblarr(nFiles)

VarsToGet=[iVar]

if (plotVectors) then begin
   VarsToGet = [VarsToGet,16,17]
   ve = fltarr(nFiles, nLons, nLats)
   vn = fltarr(nFiles, nLons, nLats)
endif

for iFile = 0, nFiles-1 do begin

   filename = filelist(iFile)

   print, 'Reading file ',filename

   file = filelist(iFile)
   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet
   AllTimes(iFile) = time

   if (CalcTEC eq 1) then begin
      tec = fltarr(nLons,nLats)

      gitm_read_bin_1var, file, He, timeHe, nVars, Vars, version, $
                          VarsToGet = [31,32]

      for iAlt = 2,nAlts-3 do begin
         tec = tec + (data(0,*,*,iAlt)-He(0,*,*,iAlt)-He(1,*,*,iAlt)) * $
               (alts(iAlt)-alts(iAlt-1))*1000.0
      endfor
      AllValues(iFile,*,*) = tec/1e16
   endif else AllValues(iFile,*,*) = data(0,*,*,iAlt)

   if (plotVectors) then begin
      ve(iFile,*,*) = data(1,*,*,iAlt)
      vn(iFile,*,*) = data(2,*,*,iAlt)
   endif

   if (strlen(diffdir) gt 0) then begin
      l1 = strpos(file,'_t')
      l2 = strpos(file,'.bin')
      length = l2-l1 +5 -2 ; +5 is for 3DALL, -2 is seconds
      filename = diffdir+'/'+strmid(file,l1-5,length)+'??.bin'
      gitm_read_bin_1var, filename, data_diff, time, nVars, Vars, version, $
                          VarsToGet = VarsToGet
      if (CalcTEC eq 1) then begin
         tec = fltarr(nLons,nLats)
         for iAlt = 2,nAlts-3 do begin
            tec = tec + data_diff(0,*,*,iAlt) * $
                  (alts(iAlt)-alts(iAlt-1))*1000.0
         endfor
         AllValues(iFile,*,*) = AllValues(iFile,*,*) - tec/1e16
         if (percentdiff) then $
            AllValues(iFile,*,*) = 100*AllValues(iFile,*,*)/(tec/1e16)

      endif else begin
         AllValues(iFile,*,*) = AllValues(iFile,*,*) - data_diff(0,*,*,iAlt)
         if (percentdiff) then $
            AllValues(iFile,*,*) = 100*AllValues(iFile,*,*)/data_diff(0,*,*,iAlt)
      endelse

      if (plotVectors) then begin
         ve(iFile,*,*) = ve(iFile,*,*) - data_diff(1,*,*,iAlt)
         vn(iFile,*,*) = vn(iFile,*,*) - data_diff(2,*,*,iAlt)
      endif

   endif

endfor

if (plotVectors) then begin
   speed = sqrt(ve^2 + vn^2)
   determine_min_max, speed, minSpeed, maxSpeed, percent=99.999
endif

determine_min_max, AllValues, mini, maxi, percent=99.999

if (CalcTEC) then colortitle = 'TEC (TECU)' else colortitle = Vars(iVar)
maxrange = 40.0

for iFile = 0,nFiles-1 do begin

   if (nFiles gt 1) then begin
      p = strpos(psfile,'.ps')
      if (p gt -1) then psfile = strmid(psfile,0,p)
      psfile_final = psfile+'_'+tostr(iFile,4)+'.ps'
   endif else begin
      psfile_final = psfile
   endelse

   value = reform(AllValues(iFile,*,*))
   sAlt = string(alts(iAlt),format='(f7.1)')
   title = colortitle + ' at '+ sAlt+' km Alt.'

   print, 'Writing file ',psfile_final
   if (plotVectors) then begin
      thermo_threeplot, psfile_final, $
                        value, AllTimes(iFile), lons, lats, mini, maxi, $
                        title, colortitle, maxrange, $
                        ve = reform(ve(iFile,*,*)), vn = reform(vn(iFile,*,*)), $
                        factor = maxSpeed/10.0
   endif else begin
      thermo_threeplot, psfile_final, $
                        value, AllTimes(iFile), lons, lats, mini, maxi, $
                        title, colortitle, maxrange
   endelse

endfor

end

