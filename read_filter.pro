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

function read_filter, fltfile, mapfile, chanID
; this function reads a filter functions from a filter file.
; Which filter to use for which channel is determined by a mapping file.
; 
; - input:  fltfile : an ASCII spreadsheat file with (number of filters + 1) columns
;                      column 1   : [number of filters,  0,    0,           0,             wavelengths in Angstrom (float)]
;                      column i>1 : [partno. filter,   CWL, FWHM, peak trans., transmission (float in between 0.0 and 1.0)]
;           mapfile : an ASCII spreadsheat file with 3 rows and (number of fibre bundles) columns
;                      row 1 : fibre bundle number (for MAST i.e. 142 to 181 => fibreID without the MD)
;                      row 2 : MSE channel number (number of the filterscope: not used, just for completeness)
;                      row 3 : Part numbers of the installed filters
;
; -returns: fltstruct: a structure with fields
;                 flt   : filter transmission for each channel/fibre bundle [nwl x nchan]
;                 lambda: wavelength vector of the filter measurement [nwl x 1]
;                 cwl   : central wavelength of the filter [nchan x 1]
;                 fwhm  : FWHM of the filter [nchan x 1]
;                 trans : peak transmission of the filter [nchan x 1]
;
;  v1.0 mdebock, 12/07/2007
;
;  v2.0 mdebock, 06/08/2008
;    Uses different filter functions for different channels (just like in the real MSE diagnostic).
;    As a result 2 input files are needed: one with the measurement of every filter
;                                          one that maps the channel/fibre bundle to the correct filter
;    And the output is a 2D array [wavelength x channels]


;-------------------------------------------------------------------------------
; read in the filter file
;-------------------------------------------------------------------------------
filter = read_txt(fltfile)
nfilt  = filter[0]
ntot   = n_elements(filter)
filter = reform(filter,nfilt+1,ntot/(nfilt+1)) ; filter file has 'nfilt+1' columns

; the filter wavelength vector
fltlambda  = filter[0,4:*]
nlambda    = n_elements(fltlambda)
;-------------------------------------------------------------------------------
; read in the mapping file
;-------------------------------------------------------------------------------
map = read_txt(mapfile)
map = reform(map, n_elements(map)/3,3) ; mapping file has 3 rows


;-------------------------------------------------------------------------------
; find the relevant filter functions from the mapping and the chanID array
;-------------------------------------------------------------------------------
nchan  = n_elements(chanID)
flt    = fltarr(n_elements(filter[0,*])-4,nchan)
cwl    = fltarr(nchan)
fwhm   = fltarr(nchan)
trans  = fltarr(nchan)
; turns the chanID string array into a numeric list (it is assumed that
; the chanID list starts with a string-prefix followed by a number. For MAST e.g. MD156)
chanNo = uint_extract(chanID)
for k=0,nchan-1 do begin
  ; find the filter partno. in the mapping array
  idx = where(map[*,0] eq chanNo[k])
  filterID = map[idx,2]
  filterID = filterID[0]  ; silly IDL pitfall: map[idx,2] is a [1x1]-array => needs converting to scalar

  ; find the filter function, based on the filter partno., in the filter array
  idx = where(filter[1:*,0] eq filterID) + 1  ; ignoring the first column (so the partno. can be equal to the
  flt[*,k] = reform(filter[idx,4:*])          ; number of filters without causing a but)

  ; we need the peak wavelength (centre of mass wavelength) of the filter.
  cumflt = total(flt[*,k],/cumulative)
  cntidx = value_locate(cumflt, 0.5*max(cumflt))
  if cntidx lt 0 then cntidx=0
  if cntidx gt nlambda-1 then cntidx=nlambda-1
  cwl[k] = fltlambda[cntidx]
  ; we need the FWHM
  mxflt = max(flt[*,k],mxidx)
  if mxidx ne 0 then begin
    tmp   = min(abs(flt[0:(mxidx-1),k] - 0.5*mxflt), halfidx1)
  endif else halfidx1=0
  tmp       = min(abs(flt[mxidx:(nlambda-1),k] - 0.5*mxflt), halfidx2)
  halfidx2 += mxidx
  fwhm[k]   = fltlambda[halfidx2]-fltlambda[halfidx1]
  ; peak transmission 
  trans[k]  = mxflt
end

; create the output structure
fltstruct = {flt: flt, lambda: fltlambda, cwl: cwl, fwhm: fwhm, trans: trans}

; and return the result
return, fltstruct

end