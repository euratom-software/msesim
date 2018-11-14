PRO calculate_bfield,bp,br,bt,bz,g

    compile_opt defint32,strictarr,strictarrsubs
    
    mw=g.mw & mh=g.mh
    bp=fltarr(mw,mh) & bt=fltarr(mw,mh) & br=fltarr(mw,mh) & bz=fltarr(mw,mh)
    dpsidx = fltarr(mw,mh)
    dpsidy = fltarr(mw,mh)
    
    ; calculate vertical derivative of psi
    for i = 0,mw-1 do begin
     dpsidy[i,*] = Deriv(g.z[0:mh-1],g.psirz[i,0:mh-1])
    endfor
    
    ; calculate horizontal derivative of psi
    for j = 0,mh-1 do begin
      dpsidx[*,j] = Deriv(g.r[0:mw-1],g.psirz[0:mw-1,j])
    endfor
    
    ; calculate array of Br, Bz, and Bp
    for j = 0,mh-1 do begin
       br[*,j] = dpsidy[0:mw-1,j]/g.r[0:mw-1]
       bz[*,j] = -dpsidx[0:mw-1,j]/g.r[0:mw-1]
    endfor
    bp = sqrt(br*br+bz*bz)
    
    ; WWH get right sign
    if g.cpasma lt 0. then begin
      br=-br & bz=-bz
    end
    
    ; Calculate toroidal field
    ; Original coding was from gfield.for by Peter Politzer,
    ;   translated to IDL by Gary Porter (see BFIELD.PRO).
    ; The code below has be optimized for IDL by Michael D. Brown, 2/16/94
    
    dpsi = (g.ssibry-g.ssimag)/float(mw-1)
    ; first order Bt value.
    for j=0,mh-1 do bt[0:mw-1,j]=g.bcentr*g.rzero/g.r[0:mw-1]  
    k = long((g.psirz - g.ssimag)/dpsi)
    iw=where(k ge 0 and k lt mw-1,n)  ; 1-d indexes where k is a valid index.
    if n gt 0 then begin
      iwr = iw mod mw  ; map matrix 1-d selected indexes to an refit row index.
      bt[iw] = ( g.fpol[k[iw]]+(g.fpol[k[iw]+1]-g.fpol[k[iw]])* $
                 (g.psirz[iw]-(k[iw]*dpsi+g.ssimag))/dpsi ) / g.r[iwr]
    endif
    
    return
end
