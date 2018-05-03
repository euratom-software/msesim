;; Copyright 2018 Euratom/CCFE/TuE

;; Permission is hereby granted, free of charge, to any person obtaining a copy of this
;; software and associated documentation files (the "Software"), to deal in the Software
;; without restriction, including without limitation the rights to use, copy, modify,
;; merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
;; permit persons to whom the Software is furnished to do so, subject to the following
;; conditions:

;; The above copyright notice and this permission notice shall be included in all copies
;; or substantial portions of the Software.

;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
;; INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
;; PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;; HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
;; CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
;; OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function ps_restore,system,noreset=noreset,close=close,_extra=ex
  
  IF N_TAGS(system) EQ 6 THEN BEGIN
    res = tag_names(system)
    found_tag = 0
    FOR i=0,5 DO BEGIN
      IF res[i] EQ 'P' THEN found_tag=found_tag+1
      IF res[i] EQ 'Y' THEN found_tag=found_tag+1
      IF res[i] EQ 'X' THEN found_tag=found_tag+1
      IF res[i] EQ 'Z' THEN found_tag=found_tag+1
      IF res[i] EQ 'C' THEN found_tag=found_tag+1
      IF res[i] EQ 'CT' THEN found_tag=found_tag+1
    ENDFOR
    IF found_tag EQ 6 THEN BEGIN
      pres = !p
      xres=!x
      yres=!y
      zres=!z
      cres=!c
  
      if not keyword_set(noreset) then begin
          device,/close
          set_plot,'X'
          print,'Closing PS device'

          !p=system.p
          !y=system.y
          !x=system.x
          !z=system.z
          !c=system.c
          IF N_TAGS(system.ct) EQ 3 THEN BEGIN
              res = tag_names(system.ct)
              found_tag=0
              FOR i=0,2 DO BEGIN
                  IF res[i] EQ 'R' THEN found_tag=found_tag+1
                  IF res[i] EQ 'G' THEN found_tag=found_tag+1
                  IF res[i] EQ 'B' THEN found_tag=found_tag+1
              ENDFOR
              IF found_tag EQ 3 THEN BEGIN
                  print,'Restoring color table'
                  tvlct,system.ct.r,system.ct.g,system.ct.b
              ENDIF ELSE print,'Can not restore color table-> wrong structure'
          ENDIF ELSE print,$
            'Can not restore color table-> needs to be a structure'
      endif else begin
          if keyword_set(close) then begin
              device,/close
              set_plot,'X'
              print,'Closing PS device'
          endif
      endelse
      return, {p:pres,y:yres,x:xres,z:zres,c:cres,ct:system.ct}
    ENDIF ELSE BEGIN
      print,'Wrong structure'
    ENDELSE
  ENDIF ELSE BEGIN
    print,'Argument needs to be a structure returned by ps_init'
  ENDELSE
end
