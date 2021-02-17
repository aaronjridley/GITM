;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf

pro makect, color

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

  if not keyword_set(color) then begin

    print,'red   - white to red'
    print,'wyr   - white to red'
    print,'blue  - white to blue'
    print,'rwb   - red white blue'
    print,'bwr   - blue white red'
    print,'bw    - black white'
    print,'wb    - white black'
    print,'mid   - blue green white yellow red'

    color = ''
    read,'Enter color table from list above : ', color

  endif

  color = mklower(color)

  ; Set read, green, blue to values normalized to the 0.0 -- 1.0 range.

  case color of
    'red' : begin
              r(*) = 1.
              g(*) = 1. - findgen(n)/(n-1)
              b(*) = 1. - findgen(n)/(n-1)
            end

    'wb' : begin
              r(*) = 1. - findgen(n)/(n-1)
              g(*) = 1. - findgen(n)/(n-1)
              b(*) = 1. - findgen(n)/(n-1)
            end

    'bw' : begin
              r(*) = 0. + findgen(n)/(n-1)
              g(*) = 0. + findgen(n)/(n-1)
              b(*) = 0. + findgen(n)/(n-1)
            end

    'blue' : begin
               r(*) = 1. - findgen(n)/(n-1)
               b(*) = 1.
               g(*) = 1. - findgen(n)/(n-1)
             end

    'rwb' : begin
              half=n/2
              r(0:half-1) = 1.
              g(0:half-1) = findgen(half)/(half-1)
              b(0:half-1) = findgen(half)/(half-1)

              r(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              b(half:n-1) = 1.
            end

    'wyr' : begin
              half=n/2
              r(0:half-1) = 1.
              g(0:half-1) = 1.
              b(0:half-1) = 1. - findgen(half)/(half-1)

              r(half:n-1) = 1.
              g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              b(half:n-1) = 0.
            end

    'bwr' : begin
              half=n/2
              b(0:half-1) = 1.
              g(0:half-1) = findgen(half)/(half-1)
              r(0:half-1) = findgen(half)/(half-1)

              b(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              r(half:n-1) = 1.
            end

    'mid' : begin
              r(0:n/3-1)     = 0.0
              r(n/3:n/2-1)   = findgen(n/2-n/3)/(n/2-n/3-1)
              r(n/2:n-1)     = 1.0

              b(0:n/2-1)      = 1.
              b(n/2:2*n/3-1)  = 1. - findgen(2*n/3-n/2)/(2*n/3-n/2-1)
              b(2*n/3-1:n-1)  = 0.

              g(0:n/3-1)      = findgen(n/3)/(n/3-1)
              g(n/3:2*n/3-1)  = 1.
              g(2*n/3:n-1)    = 1. - findgen(n-2*n/3)/(n-2*n/3-1)

	      r(n/2) = 1.0
	      g(n/2) = 1.0
	      b(n/2) = 1.0

            end

    'bgr' : begin

        if (n gt 255) then n = 255

        r = fltarr(n)
        g = fltarr(n)
        b = fltarr(n)

        r(*) = 0.0
        b(*) = 0.0
        g(*) = 0.0
        ff = findgen(n/3)/(n/3)

        n13 = 1*n/3
        n23 = 2*n/3

        r(  0:n13-1)   = 0.0
        r(n13:n23-1)   = ff
        r(n23:n  -1)   = 1.0

        g(  0:n13-1)   = ff
        g(n13:n23-1)   = 1.0
        g(n23:n  -1)   = 1.0-ff

        b(  0:n13-1)   = 1.0
        b(n13:n23-1)   = 1.0-ff
        b(n23:n  -1)   = 0.0

    end

    'all' : begin

        if (n gt 255) then n = 255

        r = fltarr(n)
        g = fltarr(n)
        b = fltarr(n)

        r(*) = 0.0
        b(*) = 0.0
        g(*) = 0.0
        ff = findgen(n/5)/(n/5)

        n25 = 2*n/5
        n35 = 3*n/5
        n45 = 4*n/5

        ; bottom 1/5 Blue to Green
        b(0:n/5-1) = 1.0 - ff
        g(0:n/5-1) = ff

        ; Next 1/5 Green to Cyan
        b(n/5:n25-1) = ff
        g(n/5:n25-1) = 1.0

        ; Next 1/5 Cyan to white
        b(n25:n35-1) = 1.0
        g(n25:n35-1) = 1.0
        r(n25:n35-1) = ff

        ; Next 1/5 White to Yellow
        b(n35:n45-1) = 1.0 - ff
        g(n35:n45-1) = 1.0
        r(n35:n45-1) = 1.0

        ; Next 1/5 Yellow to Red
        g(n45:n-1) = 1.0 - ff
        r(n45:n-1) = 1.0

    end

    else : begin
             print, "Unknown value for color=",color
             r(*) = findgen(n)
             g(*) = findgen(n)
             b(*) = findgen(n)
           end

  endcase

  r(0) = 0.0
  g(0) = 0.0
  b(0) = 0.0

  r(n-1) = 1.0
  g(n-1) = 1.0
  b(n-1) = 1.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end
