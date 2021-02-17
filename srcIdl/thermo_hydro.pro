
function limiter, dUp, dDown

  beta = 1.6
  if (dUp gt 0.0) then begin
     if (dDown gt 0.0) then begin
        limiter = min([beta*dUp, beta*dDown, (dUp+dDown)*0.5])
     endif else limiter = 0.0
  endif else begin
     if (dDown lt 0.0) then begin
        limiter = min([beta*dUp, beta*dDown, (dUp+dDown)*0.5])
     endif else limiter = 0.0
  endelse

end

pro calc_grad, v, dA, gradv

  InvDAlt = 1.0 / dA

  nV = n_elements(v)

  ; in fortran v is defined from -1 to nV+2
  ; in IDL                        0 to nV-1
  Factor1 = 0.6250000
  Factor2 = 0.0416667

  dVarLimited = v * 0.0
  vLeft = v*0.0
  vRight = v*0.0

  for i=2,nV-3 do begin  ; 1,nV in fortran
     
     h = invDAlt(i+1)*2
     dVarUp = h * (Factor1 * (v(i+1)-v(i)) - Factor2*(v(i+2)-v(i-1)))
     h = invDAlt(i)*2
     dVarDown = h * (Factor1 * (v(i)-v(i-1)) - Factor2*(v(i+1)-v(i-2)))
     dVarLimited(i) = limiter(dVarUp, dVarDown)

  endfor

  i = 1
  dVarUp            = (v(i+1) - v(i))   * InvDAlt(i+1)
  dVarDown          = (v(i)   - v(i-1)) * InvDAlt(i)
  dVarLimited(i) = Limiter(dVarUp, dVarDown)

  i = nV-2
  dVarUp            = (v(i+1) - v(i))   * InvDAlt(i+1)
  dVarDown          = (v(i)   - v(i-1)) * InvDAlt(i)
  dVarLimited(i) = Limiter(dVarUp, dVarDown)

  for i=2, nV-2 do begin
     vLeft(i)  = v(i-1) + 0.5 * dVarLimited(i-1) * dA(i)
     vRight(i) = v(i  ) + 0.5 * dVarLimited(i  ) * dA(i)
  endfor

  gradv = 0.5 * (vLeft(1:nV-1)+vRight(1:nV-1)-vLeft(0:nV-2)-vRight(0:nV-2))/dA(0:nV-2)

  gradv(0:1) = 0.0
  gradv(nV-2) = 0.0

end


file = '1DALL_t020321_*.bin'

gitm_read_bin, file, data, time, nVars, Vars, version

nTimes = n_elements(time)

rho = reform(data(*,3,0,0,*))
alt = reform(data(0,2,0,0,*))
nAlts = n_elements(alt)

da = alt
da(1:nAlts-2) = (alt(2:nAlts-1) - alt(0:nAlts-3))/2
da(0) = da(1)
da(nAlts-1) = da(nAlts-2)

o  = reform(data(*,4,0,0,*))
o2 = reform(data(*,5,0,0,*))
n2 = reform(data(*,6,0,0,*))
n1 = reform(data(*,7,0,0,*))
no = reform(data(*,8,0,0,*)) 
n = o + o2 + n2 + n1 + no
mm = (o * 16.0 + o2 * 32.0 + n2 * 28 + n1 * 14 + no * 30)/n * cmp_

t = reform(data(*,15,0,0,*))

rBody = 6372.0*1000.0
gc = 9.8
r = rBody + alt
g = -gc * (rBody/r)^2

boltz = 1.38e-23

aAll = fltarr(nTimes,nAlts-4)

for iT = 0, nTimes-1 do begin

   ts = reform(t(iT,*))
   ns = alog(reform(n(iT,*)))
   ms = reform(mm(iT,*))
   calc_grad, ns, da, gradn
   calc_grad, ts, da, gradt
   a = - (ts * gradn * boltz/ms + gradt * boltz / ms - g)

   aAll(iT,*) = a(2:nAlts-3)

endfor



end

