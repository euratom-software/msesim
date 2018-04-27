function calc_filter, fltshape, fltwidth, cwl, lambda

; this function calculates a fitler function
;
; - input:  fltshape : a string determining the filter bandshape (currently this has to be 'gauss’, ‘lorentz' or the number of cavities)
;           fltwidth : a float that gives the filter bandwidth in Angstrom (FWHM)
;           cwl      : central wavelength of the filter in Angstrom
;           lambda   : the wavelength vector
;
; -returns: flt      :    an [nlambdax1]-array with the Gaussian filter function
;                      or an [nlambdax1]-array with the Lorenzian filter function
;
;  v1.0 mdebock, 12/07/2007
;
;  v2.0 mdebock, 05/08/2008
;       because filtering goes a lot quicker with the Stokes vector approach,
;       the size of the calculated filter is now made as big as the size of the wavelength vector.
;
;  v3.0 mdebock, 17/09/2008
;       Now also a CWL is given as input parameter.
;  v3.1 agglange, 24/03/2014
;	ncav turned out to be an array [1] after the first channel-loop. now taken the first element
; convert filtershape to lowercase

fltshape = strlowcase(fltshape) 

; make the filter of that is wanted
case fltshape of
  'gauss' :  begin
               ewidth = fltwidth/(2*sqrt(alog(2)))      ; conversion from FWHM to half 1/e width
               flt    = exp(-(cwl-lambda)^2/ewidth^2)
             end
  'lorentz': begin
               flt    = (0.5*fltwidth)^2/((cwl-lambda)^2 + (0.5*fltwidth)^2)
             end
  else:      begin
               ncav   = uint_extract(fltshape)	  ; extract the number of cavities
               if ncav eq -1 then begin
                 print, ''
                 print, '   WARNING: Filter shape not recognized! 1-cavity interference filter used!'
                 print, ''
                 ncav = 1
               endif
	           ncav   = float(ncav)
	           ncav   = ncav[0]
	           flt    = (0.5*fltwidth)^(2*ncav)/((cwl-lambda)^(2*ncav) + (0.5*fltwidth)^(2*ncav))
             end
endcase

; and return the result
return, flt
end
