
filelist = findfile('3DALL_t020321_120[678]*.bin')
nFiles = n_elements(filelist)

integral = fltarr(nFiles)

for i=0,nFiles-1 do begin

   file = filelist(i)
   print, file

   gitm_read_bin, file, data, time, nVars, Vars, version

   rho = reform(data(3,*,*,*))
   alt = reform(data(2,0,0,*))
   lat = reform(data(1,0,*,0))*6372.0
   lon = reform(data(0,*,0,0))*6372.0

   nLons = n_elements(lon)
   nLats = n_elements(lat)
   nAlts = n_elements(alt)

   integral(i) = total(rho(nLons/2,2:nLats-3,48))

   if (i eq 0) then plot, rho(nLons/2, *, 48), yrange=[5e-12,10e-12] $
   else oplot, rho(nLons/2, *, 48)

endfor

end
