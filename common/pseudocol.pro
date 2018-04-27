pro pseudocol, on=on, off=off, quiet=quiet

; switches the use of pseudocolors on or off
; the color table is:
;  0: black
;  1: red
;  2: green
;  3: blue
;  4; yellow
;  5; cyan
;  6; mangenta
;  7: dark gray  (50%)
;  8: light gray (75%)
;  9: white
; 10:253: a white-yellow-orange-red-magenta-blue color table
; 254,255: white


if (keyword_set(off) && (n_elements(on) EQ 0)) then begin	; 'off' is set and 'on' is not set
  on=0
endif else begin						; 'off' is not set or 'on' is set as well ('on' rules over 'off') 
  if (n_elements(on) EQ 0) then on=1				; if 'on' is not set: turn it on
endelse


if (on EQ 1) then begin
  if ~keyword_set(quiet) then print, 'pseudocolors turned on'

  ; switch pseudocolors on
  device,true=24,decomp=0
  window,xs=5,ys=5,col=getenv('idlcols'),1,/pix
  wdelete,1

  ; initial the color table
  red256=intarr(256)
  green256=intarr(256)
  blue256=intarr(256)

  ; define the first part of the color table [0:9]
  ;         black, red, green, blue, yellow, cyan, magenta, dark gray, light gray, white
  redarr  =[    0, 255,     0,    0,    255,    0,     255,       127,        191,   255]
  greenarr=[    0,   0,   255,    0,    255,  255,       0,       127,        191,   255]
  bluearr =[    0,   0,     0,  255,      0,  255,     255,       127,        191,   255]

  red256[0:9]  =redarr
  green256[0:9]=greenarr
  blue256[0:9] =bluearr

  ; define the second part of the color table (the white-yellow-orange-red-magenta-blue color table)
  up    = findgen(61)*255.0/60.0
  black = fltarr(61)
  down  = 255 - findgen(61)*255.0/60.0
  white = fltarr(61) + 255
  redarr  = round([white, white, white, down ])
  greenarr= round([white, down , black, black])
  bluearr = round([down , black, up   , white])

  red256[10:253]  =redarr
  green256[10:253]=greenarr
  blue256[10:253] =bluearr

  ; define last part of the color table (white)
  red256[254:255]   = 255
  green256[254:255] = 255
  blue256[254:255]  = 255

  ;now load to display:
  tvlct,red256,green256,blue256
endif else begin
  if ~keyword_set(quiet) then print, 'pseudocolors turned off'

  ; switch pseudocolors off
  device,true=24,decomp=1
  window,xs=5,ys=5,1,/pix
  wdelete,1
endelse

end