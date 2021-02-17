
filename = '3DALL_t021221_000500.bin'

read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
      vars, data, rb, cb, bl_cnt, iTime, Version

o  = reform(data(4,*,*,*))
o2 = reform(data(5,*,*,*))
n2 = reform(data(6,*,*,*))
n  = reform(data(7,*,*,*))
no = reform(data(8,*,*,*))

op  = reform(data(24,*,*,*))
o2p = reform(data(25,*,*,*))
n2p = reform(data(26,*,*,*))
np  = reform(data(27,*,*,*))
nop = reform(data(28,*,*,*))

nden = o+o2+n2+n+no
eden = reform(data(33,*,*,*))

mmm  = (16*o + 32*o2 + 28*n2 + 14*n + 30*no)/nden
rhoi = (16*op + 32*o2p + 28*n2p + 14*np + 30*nop) * cmp_
rhon = reform(data(3,*,*,*))

vin = 2.6e-15 * (nden + eden)/sqrt(mmm)

fac = vin * rhoi ;/rhon

iondrag = fac * 0.0

for i=0,2 do begin
   iondrag = iondrag + fac * reform(data(36+i,*,*,*)-data(16+i,*,*,*))
endfor

r = 6372000.0 + reform(data(2,*,*,*))
dlat = data(1,0,1,0) - data(1,0,0,0)
dlon = data(0,1,0,0) - data(0,0,0,0)
lat = reform(data(1,*,*,*))


totaldrag = 0.0
for i=2,nAlts-3 do begin

   r1 = reform(r(*,*,i))
   l1 = reform(lat(*,*,i))
   id1 = reform(iondrag(*,*,i))

   dr = (r(0,0,i+1)-r(0,0,i-1))/2.0

   area = r1*r1*dr * dlat * dlon * cos(l1)

   totaldrag = totaldrag + total(id1*area)

   print, total(id1*area), totaldrag

   contour, id1, nlevels = 20

;   wait, 0.1

endfor

swforce = !pi*(cRe_*10.0)^2 * 1.0e-9

print, totaldrag/swforce

end
