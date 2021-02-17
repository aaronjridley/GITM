
filelist = "b0016_t980321_030001.3DALL"

read_thermosphere_file, filelist, nvars, nalts, nlats, nlons,vars,data

op_  = 24
n2p_ = 26
o2p_ = 25
nop_ = 28
np_  = 27
e_   = 33

numden = fltarr(nlons, nlats, nalts)
B0     = fltarr(nlons, nlats, nalts)

lat_save = reform(data(1,0,*,0,0))

for i=4,8 do numden = numden + reform(data(i,*,*,*))
mmm = 0.0
mass = [16.0,32.0,28.0,14.0,28.0]
for i=4,8 do mmm = mmm + mass(i-3) * reform(data(i,*,*,*)) / numden

altitude = reform(data( 2,*,*,*))
mmd      = reform(data( 3,*,*,*))
tn       = reform(data(15,*,*,*))
electron = reform(data(33,*,*,*))
te       = reform(data(34,*,*,*))

ec = 1.602e-19
e2 = ec ^ 2
mi = mmm * 1.6726e-27             ; pretend that mass ions = mass neutrals
me = 9.1094e-31

Vi = 2.6e-15 * (numden + electron)*(mmm^(-0.5))
Ve = 5.4e-16 * (numden)*(TE^0.5)

MeVe = me * ve
MiVi = mi * vi

B0_1d = 31000.0e-9 * (1.0 + 3.0*sin(lat_save)^2)^0.5

for i=0,nlons-1 do for k=0,nalts-1 do B0(i,*,k) = B0_1d

GyroFrequency_Ion = ec*B0/Mi
GyroFrequency_Electron = ec*B0/me

VeOe = Ve^2 + GyroFrequency_Electron^2
ViOi = Vi^2 + GyroFrequency_Ion^2


Sigma_Pedersen = ((1.0/MeVe) * (Ve*Ve/VeOe) + $
                  (1.0/MiVi) * (Vi*Vi/ViOi)) * electron * e2

Sigma_Hall = ((1.0/MeVe) * (Ve*GyroFrequency_Electron/VeOe) - $
              (1.0/MiVi) * (Vi*GyroFrequency_Ion/ViOi)) * electron * e2

Cond_Pedersen = fltarr(nlons,nlats)
Cond_Hall     = fltarr(nlons,nlats)

for k=2,nalts-3 do begin
  da = reform(altitude(*,*,k+1) - altitude(*,*,k-1))/2
  cond_pedersen(*,*) = cond_pedersen(*,*) + da * Sigma_Pedersen(*,*,k)
  cond_hall(*,*)     = cond_hall(*,*)     + da * Sigma_Hall(*,*,k)
endfor

end
