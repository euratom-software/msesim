;Version 1-jun-95  B. Rice
;file readg_jet.pro
;
;Modification of Ferron/Forest readg.pro 
;Anonymous structure is used for g so dimension can be changed
;Uses passed filename so extensions can be added g0 file
;
; 20-May-99 N Hawkes
; Modified to automatically define structure name (for easy concatenation
; of similar structures into an array)
;
;Calling format:
; g = readg(filename,[nldata=nldata])
; note: filename should include directory if g0 file is not in default dir

;Returned values:
; g : A structure called geqdsk containing the entire content of 
;     the G eqdsk file, including r,z grid points
;     The structure is defined below.
;
; g.error : indicates whether there was an error in reading the file.
;           true (or odd) indicates that there was an error.
;           false or 0 indicates no error.
;
; nldata  : anonymous structure of namelist data at the end of a formatted
;           g-file (at JET, at least).
;
;Here is some documentation from rdgeqk.for:
;
;c----------------------------------------------------------------------
;c--  CASE   : descriptive character strings                          --
;c--  XDIM   : radial dimension of rectangular grid in m              --
;c--  ZDIM   : vertical dimension in m                                --
;c--  RGRID(1):minimum major radius in m of grid                      --
;c--  ZMID   : vertical position of box center in m                   --
;c--  BCENTR : toroidal field in T at RZERO                           --
;c--  RMAXIS : major radius of magnetic axis in m                     --
;c--  ZMAXIS : vertical position of magnetic axis in m                --
;c--  SSIMAG : flux at magnetic axis in v-sec/rad                     --
;c--  SSIBRY : boundary flux                                          --
;c--  CPASMA : plasma current in Amp                                  --
;c--  FPOL   : array of poloidal current function RBt                 --
;c--  PRES   : pressure array in nt/m2 at flux grids                  --
;c--  FFPRIM : radial derivative of FPOL                              --
;c--  PPRIME : radial derivative of PRES                              --
;c--  PSIRZ  : poloidal fluxes at the (R,Z) grids                     --
;c--  QPSI   : q array at flux grids                                  --
;c--  NBDRY  : number of boundary points                              --
;c--  LIMITR : number of limiter points                               --
;c--  RBDRY, ZBDRY : (R,Z) of boundary in m                           --
;c--  XLIM, YLIM : (R,Z) of limiter                                   --
;c--  R,Z    :grid points  (calculated from g file data)
;c----------------------------------------------------------------------
;
pro geqdsk__define
  common geqdsksz, mw, mh, limmax, bdrymax

  g = { geqdsk,     $
        shot:long(0),time:long(0),error:long(0),$
        ecase:strarr(6),mw:long(0),mh:long(0),xdim:float(0),$
        zdim:float(0),rzero:float(0),rgrid1:float(0),zmid:float(0),$
        rmaxis:float(0),zmaxis:float(0),ssimag:float(0),ssibry:float(0),$
        bcentr:float(0),cpasma:float(0),$
        fpol:fltarr(mw),pres:fltarr(mw),ffprim:fltarr(mw),$
        pprime:fltarr(mw),psirz:fltarr(mw,mh),qpsi:fltarr(mw),$
        nbdry:long(0),limitr:long(0),bdry:fltarr(2,bdrymax),$
        lim:fltarr(2,limmax),R:fltarr(mw),Z:fltarr(mh)}
end

pro geqdskbig__define
  common geqdsksz, mw, mh, limmax, bdrymax

  g = { geqdskbig,     $
        shot:long(0),time:long(0),error:long(0),$
        ecase:strarr(6),mw:long(0),mh:long(0),xdim:float(0),$
        zdim:float(0),rzero:float(0),rgrid1:float(0),zmid:float(0),$
        rmaxis:float(0),zmaxis:float(0),ssimag:float(0),ssibry:float(0),$
        bcentr:float(0),cpasma:float(0),$
        fpol:fltarr(mw),pres:fltarr(mw),ffprim:fltarr(mw),$
        pprime:fltarr(mw),psirz:fltarr(mw,mh),qpsi:fltarr(mw),$
        nbdry:long(0),limitr:long(0),bdry:fltarr(2,bdrymax),$
        lim:fltarr(2,limmax),R:fltarr(mw),Z:fltarr(mh)}
end


pro geqdsk129__define
  common geqdsksz, mw, mh, limmax, bdrymax

  g = { geqdsk129,     $
        shot:long(0),time:long(0),error:long(0),$
        ecase:strarr(6),mw:long(0),mh:long(0),xdim:float(0),$
        zdim:float(0),rzero:float(0),rgrid1:float(0),zmid:float(0),$
        rmaxis:float(0),zmaxis:float(0),ssimag:float(0),ssibry:float(0),$
        bcentr:float(0),cpasma:float(0),$
        fpol:fltarr(mw),pres:fltarr(mw),ffprim:fltarr(mw),$
        pprime:fltarr(mw),psirz:fltarr(mw,mh),qpsi:fltarr(mw),$
        nbdry:long(0),limitr:long(0),bdry:fltarr(2,bdrymax),$
        lim:fltarr(2,limmax),R:fltarr(mw),Z:fltarr(mh)}
end


function readg,filename,nldata=nldata
common geqdsksz, mw, mh, limmax, bdrymax
;
;Initialize the error flag to say that there was an error.  This will
;be changed later if this routine executes successfully.
;
;----------------------------------------------------------------------
;Open the file for formatted reads.
;
on_ioerror,nofile
openr,lun,filename,/get_lun
fileformatted = 1
;
;Create the variables for the first read.
;
ecase = strarr(6)
idum = long(0)
mw = long(0)
mh = long(0)
;
;Read the first group of variables assuming that the file is formatted.
;If there is an error jump ahead to try to read the file assuming that
;it is unformatted.
;
on_ioerror,itsunformatted
;
readf,lun,ecase,idum,mw,mh,format='(6a8,3i4)'
goto,keepgo2
;
itsunformatted:
;
;Close the file then reopen it for unformatted reads.
;
on_ioerror,miscioerror
close,lun
on_ioerror,nofile
openr,lun,filename,/f77_unformatted,/segmented
;
;Create the variables for the first read.
;
ecase = strarr(6)
ecase(*) = '          '
bdum = bytarr(10,6)
idum = long(0)
mw = long(0)
mh = long(0)
;
on_ioerror,bigerror
;
readu,lun,bdum,idum,mw,mh
fileformatted = 0
ecase = string(bdum)
;
keepgo2:
;
;At this point we should have the values of mw and mh.  A typical error
;seems to be that these are 0 for some as yet undetermined reason.
;
if( (mw eq 0) or (mh eq 0) ) then return,g  ;g.error will be true (=1).
;
;----------------------------------------------------------------------
;Define the structure to be returned using mh and mw.  
;Initially the structure is filled with zeros.
;
; at least 403 and 141 needed for JET 129x129 grid
; allow for three different grid sizes by selecting one of two different
; structures (still allows structures to be combined without the clumsy
; struct_assign method).

bdrymax = 450
limmax = 200

if mw gt 65 then begin
    g = {geqdsk129}
endif else if mw gt 33 then begin
    g = {geqdskbig}
endif else begin
    g = {geqdsk}
endelse

g.error = 1
;
;
;----------------------------------------------------------------------
;create the second group of variables.
;
xdim=double(0)
zdim=double(0)
rzero=double(0)
rgrid1=double(0)
zmid = double(0)
rmaxis = double(0)
zmaxis = double(0)
ssimag = double(0)
ssibry = double(0)
bcentr = double(0)   
cpasma = double(0)
xdum = double(0)
xdum1 = double(0)
fpol = dblarr(mw)
pres = dblarr(mw)
ffprim = dblarr(mw)
pprime = dblarr(mw)
psirz = dblarr(mw,mh)
qpsi = dblarr(mw)
;
;Read the second group of variables.
;
if(fileformatted) then begin
   print,'Reading formatted g0file'
   readf,lun,xdim,zdim,rzero,rgrid1,zmid,format='(5e16.9)'
   readf,lun,rmaxis,zmaxis,ssimag,ssibry,bcentr,format='(5e16.9)'
   readf,lun,cpasma,ssimag,xdum,rmaxis,xdum,format='(5e16.9)'
   readf,lun,zmaxis,xdum,ssibry,xdum,xdum,format='(5e16.9)'
   readf,lun,fpol,format='(5e16.9)'
   readf,lun,pres,format='(5e16.9)'
   readf,lun,ffprim,format='(5e16.9)'
   readf,lun,pprime,format='(5e16.9)'
   readf,lun,psirz,format='(5e16.9)'
   readf,lun,qpsi,format='(5e16.9)'
endif else begin
   print,'Reading unformatted g0file'
   readu,lun,xdim,zdim,rzero,rgrid1,zmid
   readu,lun,rmaxis,zmaxis,ssimag,ssibry,bcentr
   readu,lun,cpasma,ssimag,xdum,rmaxis,xdum
   readu,lun,zmaxis,xdum,ssibry,xdum,xdum
   readu,lun,fpol
   readu,lun,pres
   readu,lun,ffprim
   readu,lun,pprime
   readu,lun,psirz
   readu,lun,qpsi
endelse
;
;----------------------------------------------------------------------
;Find the efit version number that was used to create the file.
;If it is large enough, then boundary and limiter data is present
;in the file so it can be read.
;
;Create the third group of variables.  These are needed even if the
;data cannot be read from the file.
;
nbdry = long(0)
limitr = long(0)
;
nvernum = long(strmid(ecase(3-1),2-1,2)+strmid(ecase(2-1),4-1,2)+$
               strmid(ecase(2-1),7-1,2))
; Y2K bug fix
if fix(strmid(ecase(3-1),2-1,2)) lt 70 then nvernum += 20000000L
;
if(nvernum ge 870520) then begin
;
;Read the third group of variables.
;
   if(fileformatted) then begin
      readf,lun,nbdry,limitr,format='(2i5)'
   endif else begin
      readu,lun,nbdry,limitr
   endelse
;
;Make the fourth group of variables.
;
   if nbdry gt 0 then bdry = dblarr(2,nbdry)
   lim = dblarr(2,limitr)
;
;Read the fourth group of variables.
;
   if(fileformatted) then begin
      if nbdry gt 0 then begin
 	 readf,lun,bdry,format='(5e16.9)'
      endif else begin
	 a='   '
         readf,lun,a		;read a null line when there is no boundary
      endelse
      readf,lun,lim,format='(5e16.9)'
   endif else begin
      if nbdry gt 0 then readu,lun,bdry
      readu,lun,lim
   endelse
;
endif
;
;----------------------------------------------------------------------
;NOTE: at this point the header information is not implemented on 
;non-VAX computers.
;
;header = string(replicate(32b,42))
;print,'header ="'+header+'"'
;;
;if(fileformatted) then begin
;   readf,lun,header,format='(a43)'
;endif else begin
;   readu,lun,header
;endelse
;;
;print,'header ="'+header+'"'
;
;----------------------------------------------------------------------
; skip the remaining variable and read the appended namelist data
;if arg_present(nldata) then begin
;  print, 'reading g-file namelist data'
;  readnl,lun,nldata
;  print, 'done'
;endif
;----------------------------------------------------------------------
;Copy the variables into the structure.
;
g.ecase = ecase
g.mw = mw
g.mh = mh
g.xdim=xdim
g.zdim=zdim
g.rzero=rzero
g.rgrid1=rgrid1
g.zmid=zmid
g.rmaxis=rmaxis
g.zmaxis=zmaxis
g.ssimag=ssimag 
g.ssibry=ssibry 
g.bcentr=bcentr 
g.cpasma=cpasma 
g.fpol=fpol
g.pres=pres 
g.ffprim=ffprim 
g.pprime=pprime 
for i=0,mh-1 do g.psirz(0:mw-1,i) = psirz(0:mw-1,i)
g.qpsi=qpsi 
g.nbdry=nbdry 
g.limitr=limitr 
if(nbdry gt 0) then begin
   if (nbdry le bdrymax) then for i=0,nbdry-1 do $
                                    g.bdry(0:1,i)=bdry(0:1,i)
endif
if(limitr gt 0) then begin
   if(limitr le limmax) then for i=0,limitr-1 do $
                                    g.lim(0:1,i)=lim(0:1,i)
endif
;

;************** added g.r,g.z calc-B. Rice *********
dR = xdim/(mw-1)
dz = zdim/(mh-1)

for i=0,mw-1 do begin
   g.r(i) = rgrid1+i*dR
endfor

for i=0,mh-1 do begin
    g.z(i) = zmid-0.5*zdim+i*dz
endfor

;*********** change sign of flux for reverse current **********
; this is currently commented out for consistency with readg.pro
;  ip_sign=g.cpasma/abs(g.cpasma)
; g.ssimag=-g.ssimag*ip_sign ;change signs here so Bz comes out with right sign
;  g.ssibry=-g.ssibry*ip_sign
;  g.psirz=-g.psirz*ip_sign 

;************************
;
;----------------------------------------------------------------------
;Get the actual shot number and time from the G file data.
;

nine = strmid(ecase(4-1),3-1,1)
if(nine ne ' ') then begin
   nshot6=strmid(ecase(4-1),3-1,6)
   g.shot = long(nshot6)
endif else begin
   nshot5=strmid(ecase(4-1),4-1,5)
   g.shot = long(nshot5)
endelse
;
;nntime = strmid(ecase(5-1),3-1,4)
;g.time = long(nntime)
g.time=long(strmid(ecase(4),1,5))
;
;----------------------------------------------------------------------
;The routine executed successfully if the boundary or limiter arrays
;did not overflow.
;
if( (nbdry le bdrymax) and (limitr le limmax) ) then g.error = 0
;
;----------------------------------------------------------------------
;close the file
;
alldone:
   on_ioerror,miscioerror
   close,lun
   on_ioerror,miscioerror
   free_lun,lun
;
finishup:
return,g
;
;----------------------------------------------------------------------
;errors
;
nofile:
   print,'In readg, the file "'+filename+'" could not be located.'
   goto,finishup
bigerror:
   print,'Error reading the file "'+filename+'" in readg.'
   goto,alldone
miscioerror:
   print,'I/O error in readg'
   goto,finishup
;
end
