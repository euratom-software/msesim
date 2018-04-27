pro create_coll_param, filename=filename
; reads in the IDL xdr-file create by Clive's LEON code, asks for some reader input and 
; outputs it into an xml-file that can be read by the STARK simulation code.

set_plot,'x'
pseudocol

; first the IDL xdr-file is restored
print, ''
print, 'reading file: ', filename
restore, filename
print, ''

; print lens position
p0 = p0/1e3
print, ' xyz: ', p0, ' (machine coordinates [m] of the collection lens)'
print, ''

; print effective focal length
print, ' efl: ', efl, ' mm'
print, ' (effective focal length)'
print, ''

; the hatz-vector is transformed into the k-vector,
; the l- and m-vectors are created and compared to the xhat and yhat vectors
k = -zhat
l = crossp([0.,0.,1.],k)/norm(crossp([0.,0.,1.],k))
m = crossp(k,l)

print, ' k-vector: ', k, ' (optical axis, from emission - towards lens)'
print, ' l-vector: ', l, ' (horizontal & perpendicular to optical axis: l= zxk/|zxk|)'
print, ' m-vector: ', m, ' (perpendicular to optical axis & horizontal: m = kxl)'
print, ''

if abs(total(l-xhat)) gt 0.01 then begin
  print, 'l-vector does not correspond with xhat-vector!'
  print, ' l   : ', l
  print, ' xhat: ', xhat
  ans=''
  read, prompt='Do you wish to continue (Y/N)? (l-vector will be used) ', ans
  if strcmp(ans,'n',1,/fold_case) then return
endif
if abs(total(m-yhat)) gt 0.01 then begin
  print, 'm-vector does not correspond with yhat-vector!'
  print, ' m   : ', m
  print, ' yhat: ', yhat
  ans=''
  read, prompt='Do you wish to continue (Y/N)? (m-vector will be used) ', ans
  if strcmp(ans,'n',1,/fold_case) then return
endif

; print fibre bundle layout
nbund = n_elements(x0[*,0])   ; number of fibre bundles
nfib  = n_elements(x0[0,*])   ; number of fibres/bundle
print, ' # fibre bundles: ', nbund
print, ' # fibres/bundle: ', nfib
print, ''

; convert x0 and y0 to focal plane
x0 = -x0
y0 = -y0

; check whether the first bundle is most inboard or most outboard.
; in the latter case: reverse bundle order 
; (i.e. increasing bundleno. <=> increasing tangency radius)
Rt0 = tangency(x0[0,*],p0,k,l,efl)
Rtn = tangency(x0[nbund-1,*],p0,k,l,efl)
print, ''
print, ' tangency radius first bundle: ', Rt0, ' m'
print, ' tangency radius last bundle : ', Rtn, ' m'
if Rt0 lt Rtn then print, ' => OK' $
else begin	; reversing first dimension of x0 and y0
  x0 = reverse(x0,1)
  y0 = reverse(y0,1)
  print, ' => reversing fibre bundle order (lowest tangency radius first)'
endelse


; ask for the fibre diameter
print, ''
ans=0.0
read, prompt='Please enter the total fibre diameter in micron      : ', ans
df = float(ans)/1e3
read, prompt='Please enter the diameter of the fibre core in micron: ', ans
cf = float(ans)/1e3

; ask for the lens diameter
print, ''
ans=0.0
read, prompt='Please enter the diameter of the collection lens [mm]: ', ans
dl = float(ans)

; plot the current settings
plot_coll, p0, efl, k,l, x0, y0, df


; allow tweaking of the effective focal length
o=1 
while o do begin 
  print, ''
  prompt = string(format='("Are you happy with the effective focal length (",F6.2," mm) ? (Y/N) ")',efl)
  ans=''
  read, prompt=prompt, ans
  if strcmp(ans,'n',1,/fold_case) then begin
    read, ' => Enter a new effective focal length [mm]: ', efl
    plot_coll, p0, efl, k,l, x0, y0, df
  endif else o=0
endwhile 

; allow tweaking of the lens position
o=1 
while o do begin 
  print, ''
  prompt = string(format='("Are you happy with the lens position (",F5.2,",",F5.2,",",F5.2,") [m] ? (Y/N) ")',$
                  p0[0],p0[1],p0[2])
  ans=''
  read, prompt=prompt, ans
  if strcmp(ans,'n',1,/fold_case) then begin
    read, ' => Enter new x,y,z lens coordinates [m]: ', p0
    plot_coll, p0, efl, k,l, x0, y0, df
  endif else o=0
endwhile 


; allow tweaking of the optical axis
o=1 
while o do begin 
  print, ''
  prompt = string(format='("Are you happy with the k-vector/optical axis  (",F5.2,",",F5.2,",",F5.2,") [m] ? (Y/N) ")',$
                  k[0],k[1],k[2])
  ans=''
  read, prompt=prompt, ans
  if strcmp(ans,'n',1,/fold_case) then begin
    read, ' => Enter new x,y,z k-vector coordinates (does not need to be normalised): ', k
    k = k/norm(k)
    l = crossp([0.,0.,1.],k)/norm(crossp([0.,0.,1.],k))
    m = crossp(k,l)
    plot_coll, p0, efl, k,l, x0, y0, df
  endif else o=0
endwhile 

; Get the fibre bundle ID
print, ''
print, format='("There are ",I2," fibre bundles")', nbund
pre   = ''
read, prompt=' => Please enter a prefix for the fibre bundle ID (e.g. "MD"): ', pre
start = 1
read, prompt=' => Please enter the start counter for the fibre bundle ID (e.g. "142"): ', start

; that's enough tweaking, lets save the data in a xml file
print, ''
fname= ''
read, prompt='Please enter an xml-filename to save this data to: ', fname

; open the file for writing
openw, fileID, fname, /GET_LUN 
; and start writing:
printf, fileID, '<stark_settings description="Input parameters for the IDL STARK code">'
printf, fileID, '  <coll description="Data relevant to the Collection optics">'
printf, fileID, '    <xyz'
printf, fileID, '      description="The xyz-coordinates - in the machine coordinate system - of the collection lens."'
printf, fileID, '      unit="m"'
printf, fileID, '      type="float"'
printf, fileID, format='(A,F6.3,", ",F6.3,", ",F6.3,A)',$
                '      value="',p0[0],p0[1],p0[2],'"/>'
printf, fileID, '    <k'
printf, fileID, '      description="The xyz-coordinates - in the machine coordinate system - of the k-vector of the collection lens.' 
printf, fileID, '      This is the vector describing the optical axis and in the direction from the emission volume towards the'
printf, fileID, '      collection lens. REMARK: the collection lens coordinate system is given by: (k,l,m), with k the unit vector'
printf, fileID, '      given here, l=zxk/|zxk| - z being the z of the machine coordinate system - and m=kxl. This means that'
printf, fileID, '      k-coordinates are along the optical axis, l coordinates are horizontal and m-coordinates are - usually -'
printf, fileID, '      more or less vertical (if the optical axis is more or less horizontal)."'
printf, fileID, '      type="float"'
printf, fileID, format='(A,F6.3,", ",F6.3,", ",F6.3,A)',$
                '      value="',k[0],k[1],k[2],'"/>'
printf, fileID, '    <d'
printf, fileID, '      description="The diameter of the collection lens."'
printf, fileID, '      unit="mm"'
printf, fileID, '      type="float"'
printf, fileID, format='(A,F5.2,A)',$
                '      value="',dl,'"/>'
printf, fileID, '    <efl'
printf, fileID, '      description="Effective focal length of the collection lens (This is linked with the alpha - the opening angle'
printf, fileID, '      that corresponds with 1 mm on the fibre plate - through: efl = 1/tan(alpha) )."'
printf, fileID, '      unit="mm"'
printf, fileID, '      type="float"'
printf, fileID, format='(A,F5.2,A)',$
                '      value="',efl,'"/>'
printf, fileID, '    <fibre'
printf, fileID, '      description="Diameter of the optical core of the fibre and the diameter of the total fibre."'
printf, fileID, '      unit="micrometer"'
printf, fileID, '      type="float"'
printf, fileID, format='(A,F6.1,", ",F6.1,A)',$
                '      value="',cf*1e3,df*1e3,'"/>'
printf, fileID, '    <bundleID'
printf, fileID, '      description="The IDs of each fibre bundle. E.g. for MAST the MSE fibre bundles are named: MD142, MD142, ...,'
printf, fileID, '      MD181. REMARK: The xml-nodes further down in this file describing l- and m-coordinates of the fibres in each'
printf, fileID, '      fibrebundle should be named with this ID and the number of nodes should correspond with (twice, one for l-,'
printf, fileID, '      one for m-coordinates) the number of elements in this ID-list."'
printf, fileID, '      type="string"'
printf, fileID, format='($,A)',$
                '      value="'
for i=0,nbund-1 do begin
  printf, fileID, format='($,A,I03)',$
                  pre,start+i
  if i ne nbund-1 then printf, fileID, format='($,", ")' $
  else printf, fileID, '"/>'
endfor

for i=0,nbund-1 do begin
  printf, fileID, format='(A,A,I03,"l")',$
                '    <',pre, start+i
  printf, fileID, '      description="l-coordinates in mm of the fibres in this fibrebundle"'
  printf, fileID, '      unit="mm"'
  printf, fileID, '      type="float"'
  printf, fileID, format='($,A)',$
                  '      value="'
  for j=0,nfib-1 do begin
    printf, fileID, format='($,F7.3)', x0[i,j]
    if j ne nfib-1 then printf, fileID, format='($,", ")' $
    else printf, fileID, '"/>'
  endfor

  printf, fileID, format='(A,A,I03,"m")',$
                '    <',pre, start+i
  printf, fileID, '      description="m-coordinates in mm of the fibres in this fibrebundle"'
  printf, fileID, '      unit="mm"'
  printf, fileID, '      type="float"'
  printf, fileID, format='($,A)',$
                  '      value="'
  for j=0,nfib-1 do begin
    printf, fileID, format='($,F7.3)', y0[i,j]
    if j ne nfib-1 then printf, fileID, format='($,", ")' $
    else printf, fileID, '"/>'
  endfor
endfor
printf, fileID, '  </coll>'
printf, fileID, '</stark_settings>'
; close the file
close, fileID
free_lun, fileID



end



;----------------------------------
; plot routine
;----------------------------------
pro plot_coll, p0, efl, k,l, x0, y0, df
; plots the lines-of-sight, as a double check

window, xsize=1000, ysize=750
!p.multi=2

; plot the vessel contours
Rvessel = [0.2,2.0]
phi     = findgen(120)/119.*2.*!pi
plot, Rvessel[0]*cos(phi),Rvessel[0]*sin(phi),color=9, /iso,$
      xs=1,xr=[-2.1,2.1],ys=1,yr=[-2.3,1.0],xtitle='X [m]',ytitle='Y [m]',$
      position=[0.05,0.3,0.95,0.95]
oplot, Rvessel[1]*cos(phi),Rvessel[1]*sin(phi),color=9

; plot optical axis and lens
d = 0.038	; lens diameter
e = 3.0		; length of the l.o.s.
oplot, [p0[0],p0[0]-e*k[0]], [p0[1],p0[1]-e*k[1]],$
       color=1,thick=3
oplot, [p0[0]-d/2.*l[0], p0[0]+d/2.*l[0]], [p0[1]-d/2.*l[1], p0[1]+d/2.*l[1]],$
       color=5,thick=3

; loop through all bundles and plot topview
nbund = n_elements(x0[*,0])   ; number of fibre bundles
nfib  = n_elements(x0[0,*])   ; number of fibres/bundle
for i=0,nbund-1 do begin
  ; emission-vectors on each line of sight:
  z0  = 0.*x0[i,*] + efl
  klm = [z0,x0[i,*],y0[i,*]]
  klm = klm/1e3                                                   ; convert from mm into m
  vklm = klm/rebin(sqrt(klm[0,*]^2+klm[1,*]^2+klm[2,*]^2),3,nfib) ; unit vectors

  ; convert into klm into xyz-coordinates and vklm into xyz-vectors
  O1   = coordtrans([0.,0.,0.],[[p0],[k],[l]])
  X1   = coordtrans([1.,0.,0.],[[0,0,0],[k],[l]])
  Y1   = coordtrans([0.,1.,0.],[[0,0,0],[k],[l]])
  xyz  = coordtrans(klm, [[O1],[X1],[Y1]])
  vxyz = coordtrans(vklm, [[0.,0.,0.],[X1],[Y1]])

  ; plot all lines-of-sight
  color = (i mod 2) ? 4 : 9
  for j=0,nfib-1 do begin
    oplot, [xyz[0,j],xyz[0,j]-e*vxyz[0,j]], [xyz[1,j],xyz[1,j]-e*vxyz[1,j]], color=color
  endfor
endfor

; loop through all bundles and plot focal plane
phi     = findgen(30)/29.*2.*!pi
plot, [min(x0)-0.1*(max(x0)-min(x0)),max(x0)+0.1*(max(x0)-min(x0))],$
      [min(y0)-0.1*(max(y0)-min(y0)),max(y0)+0.1*(max(y0)-min(y0))],$
      /nodata,/iso,xs=1,ys=1, title='Focal plane', xtitle='l [mm]',ytitle='m [mm]',$
      position=[0.05,0.05,0.95,0.3]
corner = [min(x0)-1.,max(y0)-0.5]
oplot, corner[0]+0.5*cos(phi),corner[1]+0.5*sin(phi), color=9
polyfill, corner[0]+0.1*cos(phi),corner[1]+0.1*sin(phi), color=9
xyouts, min(x0)-1.0,max(y0)-1.5, 'k'
for i=0,nbund-1 do begin
  ; plot all lines-of-sight
  color = (i mod 2) ? 4 : 9
  for j=0,nfib-1 do begin
    polyfill, x0[i,j]+df/2.*cos(phi), y0[i,j]+df/2.*sin(phi), color=color
  endfor
endfor

end


;----------------------------------
; tangency radius function
;----------------------------------
function tangency, x0,p0,k,l,efl
; calculates tangency radius of the fibre bundle

; 'central' fibre emission vector in klm-coordinates
klm = [efl, mean(x0),0.0]
klm = klm/norm(klm)            ; unit vectors
; convert into klm into xyz-coordinates
X1  = coordtrans([1.,0.,0.],[[0,0,0],[k],[l]])
Y1  = coordtrans([0.,1.,0.],[[0,0,0],[k],[l]])
xyz = coordtrans(klm, [[0.,0.,0.],[X1],[Y1]])

; distance to p0:
D   = norm(p0)
; vector to p0:
v   = p0/D
; angle between xyz and v:
alpha = acos(dotp(v,xyz))
; tangency radius:
Rt  = D*sin(alpha)

return, Rt




end
