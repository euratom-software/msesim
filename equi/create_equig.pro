PRO create_equig, gfile, fluxold=fluxold

;Creates an MSESIM equilibrium file from a g.eqdsk file

;optionally, it uses the old flux-coords, because spectrum_calc wants ...
;currently (2013-10-12) it stops at spectrum_calc lines 958, with a fluxgrid that does not look nice.
;
; UPDATE: Now the poloidal flux is normalized and it seems to work. If things keep being allright, in the next update the fluxold-option will be removed.
;
; 2013-10-12, AGGL


;directory separators
;----------------------------------------------------------------------------------
if strcmp( !version.os,'Win32',/fold_case) then sep='\' else sep='/'

; Get the main directory (where the folder "equi" is assumed to be located)
proname = 'create_equig'
help, /source, name=proname,output=output
outputstr=''
for i=0,n_elements(output)-1 do outputstr+=output[i]
first = strpos(outputstr,  strupcase(proname) )+strlen(proname)+1
last  = strpos(outputstr,'Compiled Functions:',/reverse_search)
length= last-first
sourcefile = strtrim(strmid(outputstr,first,length),2)
sourcedir  = strmid(sourcefile,0,strpos(sourcefile, sep, /reverse_search))
maindir = strmid(sourcedir,0,strpos(sourcedir, sep, /reverse_search))


print, 'sourcedir = ' + sourcedir
print, 'sourcefile = ' + sourcefile
print, 'maindir = ' + maindir


; Read the gfile

g=readg(gfile) ; returns the g.eqdsk file as a structure

psirz = g.psirz ; the fluxcoords
arrsize = n_elements(psirz[*,0])


if (keyword_set(fluxold)) then begin

	if arrsize eq 65 then begin

		restore, maindir + sep + 'equi' + sep + 'equi_KSTAR_Te=2500eV_ne=5e19m3_Ip=0.5MA_B=2T.sav'

		fluxcoord=fluxcoord ; we use this one

	endif else begin

		fluxcoord=congrid(fluxcoord, arrsize, arrsize, /interp) ; rebin (non-integer)

	endelse

	fluxoldflag = '_oldflux'

endif

if ~(keyword_set(fluxold)) then begin

		fluxcoord = (psirz-g.ssimag)/(g.ssibry-g.ssimag)

		fluxoldflag = ''
endif




; Calculate the B-field

calculate_bfield, bp,br,bt,bz,g

bfld=fltarr(3,arrsize,arrsize)

bfld[0,*,*]=br
bfld[1,*,*]=bz
bfld[2,*,*]=bt

Rm=g.rzero
R=g.R
Z=g.Z

gfile=strmid(gfile,strpos(gfile,sep,/reverse_search)+1,15)
shot = strmid(gfile,3,4)
time = strmid(gfile,10,4)

; save this into the given file
save, R,Z,Bfld,fluxcoord,Rm, filename=maindir+sep+'equi'+sep+'equi_K'+shot+'_'+time+fluxoldflag+'.sav'


end
