;file = '2DGEL_t020321_000100.bin'
;  read_thermosphere_file, file, nvars, nalts, nlats, nlons, $
;                          vars, data, rb, cb, bl_cnt, iTime, Version

;stop
  
file = '2DMEL_t020321_000000.bin'
;file = '2DMEL_t130308_180001.bin'
  read_thermosphere_file, file, nvars, nalts, nlats, nlons, $
                          vars, data, rb, cb, bl_cnt, iTime, Version

  mlat = reform(data(1,*,*)) / !dtor
  mlon = reform(data(0,*,*)) / !dtor
  mlt = reform(data(3,*,*)) / !dtor

  pot = reform(data(24,*,*))

  iLat = 90

  mlt1d = reform(mlt(*,iLat))
  pot1d = reform(pot(*,iLat))

;  sort_a, mlt1d, pot1d
  
  plot, mlt1d, pot1d

;stop
  
  n = n_elements(pot1d)
  e = fltarr(5,n)
  

  for i=88,92 do begin
     pot1d = reform(pot(*,i))
     
     e(i-88,1:n-2) = -(pot1d(2:n-1) - pot1d(0:n-3))/2
     e(i-88,0) = -(pot1d(1)-pot1d(0))
     e(i-88,n-1) = -(pot1d(n-1)-pot1d(n-2))
  endfor
  
  plot, mlt1d, e(3,*)
  

end

  
