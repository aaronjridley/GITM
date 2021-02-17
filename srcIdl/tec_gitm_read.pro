
pro tec_gitm_read, file, lats, lons, localtime, tec, gitm, time

  openr,1,file
  itime = intarr(6)
  readf,1,itime
  c_a_to_r, itime, time

  ut = (time mod 86400.0) / 3600.0

  readf,1,nPts
  tmp = fltarr(4,nPts)
  readf,1,tmp
  lats = reform(tmp(0,*))
  lons = reform(tmp(1,*))
  tec  = reform(tmp(2,*))
  gitm = reform(tmp(3,*))
  
  localtime = (lons/15.0 + ut) mod 24.0

  close,1

end
