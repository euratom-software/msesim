pro spectrum_calc, settingfile, outputfile

; routine calculates the raw spectrum for a given inputfile (xml-format)
; and returns the raw spectrum data in the given outputfile (xdr-format)
;
; the raw spectrum data consists of the following
;   - the 'central points' data: channel ID (e.g. ['MD150','MD151',...] for MAST)                  - chanID
;                                coordinates of the central points                                 - xyz0
;                                unmodulated beam velocity                                         - B_v0
;                                beam duct coordinates                                             - B_xyz
;                                half beam sampling width                                          - B_w
;                                (unity) beam vector (beam axis)                                   - B_vec
;                                collection lens coordinates                                       - C_xyz
;                                (unity) collection lens vector (optical axis)                     - C_k
;                                (unity) emission vector of the central points                     - vec0
;                                B-field vector at the central points                              - Bfld0
;                                E-field vector at the central points                              - Efld0
;                                Doppler shift at the central points                               - Dshift0
;                                Max. Stark shift at the central points                            - Sshift0
;                                Polarisation angles at the central points                         - alpha0
;                                normalised psi at the central points                              - psi0
;
;   - the 'grid points' data   : coordinates of the grid points                                    - gp_xyz
;                                (unity) beam velocity vector at the grid points                   - gp_vel
;                                (unity) emission vector at the grid points                        - gp_vec
;                                B-field vector at the grid points                                 - gp_Bfld
;                                E-field vector at the grid points                                 - gp_Efld
;                                Polarisation angles at the grid points                            - gp_alpha
;                                normalised psi at the grid points                                 - gp_psi
;                                emission intensity at the grid points                             - gp_emis
;
;   - the 'resolution' data    : R, Z and psi vectors                                              - R, Z, psi
;                                psi in the RZ-plane and the magnetic axis                         - RZpsi, Rm
;                                emission intensity in the RZ-plane                                - RZ_emis
;                                emission intensity as a function of midplane R                    - R_emis
;                                emission intensity as a function of psi                           - psi_emis
;                                resolution vector as a function of psi                            - psi_res
;                                resolution vector as a function of R                              - R_res
;
;   - the 'spectral' data      : wavelength vector [nwlx1]                                         - lambda
;                                pi Stokes vector [nwlx4xnchan]                                    - pstokes0
;                                sigma Stokes vector [nwlx4xnchan]                                 - sstokes0
;                                total Stokes vector [nwlx4xnchan]                                 - stokes   ( = pstokes0 + sstokes0 )
;                                Stokes vector at the central wavelength [5xnchan]                 - cwlstokes
;                                Stokes vector at the optimal sigma wavelength [5xnchan]           - sigstokes
;                                Stokes vector at the optimal blue shifted pi wavelength [5xnchan] - pibstokes
;                                Stokes vector at the optimal red shifted pi wavelength [5xnchan]  - pirstokes
;                                (5 element contains the central wavelength)
;   - a 'filterflag'           : put to 0 indicating that this data is                             - filterflag
;                                'unfiltered'. The 'spectrum_filter.pro' routine
;                                will set this flag to 1 (after filtering the data)
;
; v1.0, mdebock 17/07/2007
; v2.0, mdebock 21/04/2008 : * The calculation of the Stark spectrum now uses Stokes vectors. This increases
;                              the speed of the code and its accuracy. However, it does mean changes were
;                              made in the structure of both setting and output files. v2.0 is therefore
;                              no longer downwards compatible with v1.0
; v2.1, daussems 29/09/2011: * Stokes vectors a CWL and position optimal pi-red, sigma and pi-blue now determined (for profile plots)
;       mdebock              * Windows/Unix compatibility ensured
;
; v2.2, pgeelen 09/02/2012:  * Inclusion of the Paschen-back effect. 
;       mdebock                Stokes vector (4D) now incorporates circular polarisation term.
;
; v2.3, pgeelen 21/03/2012:  * Inclusion of the correction for the non statistical population of the energy levels. 
;       mdebock                
;
; v2.4, mdebock 10/01/2014:  * Different channels can now have a different number of fibres. 
;                

;----------------------------------------------------------------------------------
; Trap every math error (for debugging)
;----------------------------------------------------------------------------------
!EXCEPT=2

;----------------------------------------------------------------------------------
; Distinguish between Windows and UNIX directory separators
;----------------------------------------------------------------------------------
if strcmp(!version.os,'Win32',/fold_case) then sep='\' else sep='/'


;**********************************************************************************
;* READ INPUT PARAMETERS                                                          *
;**********************************************************************************

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen (initialisation begins)
;----------------------------------------------------------------------------------
runinfo = strarr(9)
runinfo[0:3] = strsplit(settingfile,sep,/extract)
runinfo[4:8] = strsplit(outputfile,sep,/extract)
print,'MSE spectrum calculation:'
print,'-------------------------'
print,'   run         : ' + runinfo[1]
print,'   setting file: ' + runinfo[3]
print,'   output file : ' + runinfo[8]
print,''
print,FORMAT='($,A)',"* Initialising ....................."
print,FORMAT='($,A)',"................................... "

;----------------------------------------------------------------------------------
; Get the physical constants from the 'physic_constants.xml'-file in the
; 'physics' directory
;----------------------------------------------------------------------------------
common physics, c, kb, em, e, pm, nm, RydInf

;fname = ['/home/sgibson/PycharmProjects/msesim/equi/physics/','physic_constants.xml']
cwd = getenv('PWD')
fname = [cwd+'/equi/physics/','physic_constants.xml']

physics  = readxml(fname)
c        = physics.c		; light speed (m/s)
kb       = physics.kb		; Boltzmann constant (J/K)
em       = physics.em		; electron mass (kg)
e        = physics.e		; electron charge (C)
pm       = physics.pm		; proton mass (kg)
nm       = physics.nm		; neutron mass (kg)
RydInf   = physics.RydInf	; Infinity Rydberg constant (1/m)
nAir	 = physics.nAir		; Refractive index of air
eps0     = physics.eps0		; vacuum permittivity [1e-12 F/m]
h        = physics.h/physics.e  ; Planck constant     [eV.s]


;----------------------------------------------------------------------------------
; Read the input settings for calculating the Stark spectrum
input = read_setting(settingfile)

;----------------------------------------------------------------------------------
; Put equilibrium data from the input.equi-structure into variables
;----------------------------------------------------------------------------------
R0       = input.equi.R0		; Major radius of the centre of the plasma (m)
a        = input.equi.a			; Minor radius of the plasma in the midplane (m)

equifile = input.equi.file		; Name of the file that contains the results of the

if ~strcmp('none',equifile,/fold_case) then begin
  equifile = cwd+'/equi/' + STRING(equifile); equilibrium code. If the filename is set to 'none'
endif					; than the build-in equilibrium model is used.

; Following parameters are needed for the build-in equilibrium model 
; (only used if the equifile is set to 'none')
Bphi     = input.equi.Bphi		; Toroidal magnetic field in vacuum at R0 (T)
q0       = input.equi.q0		; q at the magnetic axis
qa       = input.equi.qa		; q at the plasma edge
qidx     = input.equi.qidx		; q-index: q(psi) = q0 + (qa-q0)*sqrt(psi)^qidx
Bp0      = input.equi.Bpar0		; paramagnetic field at the magnetic axis
Bpa      = input.equi.Bpara		; paramagnetic field at the plasma edge
Bpidx    = input.equi.Bparidx		; Bp(psi) = Bp0 + (Bpa-Bp0)*sqrt(psi)^Bpidx
shafr	 = input.equi.shafr		; Shafranov shift (m)
elong    = input.equi.elong		; Elongation


;----------------------------------------------------------------------------------
; Put neutral beam data from the input.beam-structure into variables
;----------------------------------------------------------------------------------
fixset   = input.beam.fixedsettings	; the xml-file that contains beam settings
fixset   ='beam'+ sep + fixset		; that are unlikely to change.
fixed    = read_setting(fixset)
B_xyz    = fixed.beam.xyz		; xyz-coordinates of the beam duct (m)
B_xi     = fixed.beam.xidelta[0]	; Angle between x-axis and beam axis (degrees)
B_delta  = fixed.beam.xidelta[1]	; Angle between z-axis and beam axis (degrees)
B_dS	 = fixed.beam.dS		; Distance beam source to the beam duct (m)
B_hf	 = fixed.beam.focus[0]		; Horizontal focus of the beam, from beam source (m)
B_vf     = fixed.beam.focus[1]		; Vertical focus of the beam, from beam source (m)
pinifile = fixed.beam.pinifile		; file containing the pini-beamlet coordinates
pinifile ='beam'+sep + pinifile
B_w      = fixed.beam.w			; The half sampling width of the beam, also full
					; 1/e-width if internal beam model is used (m)
B_div    = fixed.beam.div		; Half 1/e width of the beam divergence (degrees)
B_M      = fixed.beam.mass		; Mass number of beam particles (1=H, 2=D, 3=T)
B_tEf    = fixed.beam.Vmod_type		; Type of the beam energy fluctuation
B_aEf    = fixed.beam.Vmod_amp		; Peak-to-Peak amplitude of the fluctuation (keV)
B_nEf    = fixed.beam.Vmod_n		; Number of energy samples to average
B_fEf    = fixed.beam.Vmod_file		; file containing a histrogram of the beam ripple (only if Vmod_type is 4)
if ~strcmp('none',B_fEf,/fold_case) then begin
  B_fEf  ='beam'+sep + B_fEf
endif
Qion     = fixed.beam.Qion		; Ionisation rate of beam particles per m (m^2)
Qemit    = fixed.beam.Qemit		; Emission rate of beam particles per s (photons * m^3/s)

; following beam settings are more likely to change so they come from the main setting file:
B_E      = input.beam.V			; Beam energy (keV)
B_E     *= 1e3				; convert from keV to eV
beamfile = input.beam.file		; Name of the file that contains the results
if ~strcmp('none',beamfile,/fold_case) then begin
  beamfile ='beam' + sep + beamfile	; of M. Turnyanskiy's beam code. If the filename is
endif					; set to 'none' than the build-in beam model is used

; Following parameters are needed for the build-in beam model
; (only used if the beamfile is set to 'none')
Bdens0   = input.beam.Bdens0		; Beam density at the plasma edge (m^(-3))
edens0   = input.beam.edens0		; Electron density on axis (m^(-3))

;----------------------------------------------------------------------------------
; Put settings for the spectra from the input.spectrum-structure into variables
;----------------------------------------------------------------------------------
trans    = input.spectrum.trans		; Type of transition (1=H-alpha, 2=H-beta, 3=H-gamma)
atomdata = input.spectrum.atomdata	; Name of the file containing the atomic data.
atomdata = '..'+sep+'physics'+sep+'stark'+sep+atomdata	; This file is located in the 'physics' directory

resfrac  = input.spectrum.resfrac	; The fraction of the emission used for the spatial resolution

Erfile   = input.spectrum.Erfile	; Name of the file that contains the radial electric field in the plasma
if ~strcmp('none',Erfile,/fold_case) then begin
  Erfile ='Er'+sep+ Erfile		; If the filename is set to 'none' than the build-in beam model is used
endif
Ermax    = input.spectrum.Ermax		; The maximum value of the radial electric field. It is assumed this maximum 
					; is located at half radius (psi=0.5), and that the Er-field is zero at the 
					; magnetic axis (psi=0) and at the plasma edge (psi=1)
Eridx    = input.spectrum.Eridx		; Defines the shape of the the radial electric field profile: 
					; Er(psi) = Ermax*(1 - (2*abs(0.5-sqrt(psi)))^Eridx)

;----------------------------------------------------------------------------------
; Put collection lens data from the input.coll-structure into variables
;----------------------------------------------------------------------------------
collfile = input.coll.file		       ; the xml-file that contains collection optics settings
collfile ='coll_optics'+sep + collfile ; (which are unlikely to change).
collset  = read_setting(collfile)
C_xyz    = collset.coll.xyz		       ; xyz-coordinates of the collection lens (m)
C_k      = collset.coll.k		       ; vector 'k' of the optical axis (from the emission region towards
C_k      = C_k/norm(C_k)		       ; the lens. This should be a unit vector)
C_l      = crossp([0.,0.,1.],C_k)	   ; and the other 2 vectors of the collection lens coordinate system:
C_l      = C_l/norm(C_l)		       ; the 'horizontal' vector: l=zxk/|zxk|
C_m      = crossp(C_k,C_l)		       ; the 'vertical' vector  : m= kxl
C_m      = C_m/norm(C_m)
C_w      = collset.coll.d		       ; diameter of the collection lens (mm)
C_efl    = collset.coll.efl		       ; effective focal length of the lens (mm)
C_fibre  = collset.coll.fibre		   ; C_fibre[0] = diameter of the optical core of the fibre (micrometer)
					                   ; C_fibre[1] = diameter of the total fibre (core+clad)   (micrometer)
C_bundID = collset.coll.bundleID       ; string array with the names of all fibre bundles. The number of
nchan    = n_elements(C_bundID)		   ; elements corresponds with the total number of bundles/channels

;----------------------------------------------------------------------------------
; Put settings for the numerical integration from the input.integration-structure into variables
;----------------------------------------------------------------------------------
bundles  = input.integration.bundles	; list of the fibre bundles for which to calculate a Stark spectrum
; convert the into a numeric channel list (first bundle = 0, et cetera)
tmpstr=''
if strcmp(bundles[0],'-1',/fold_case) then begin
  channels=indgen(nchan)                ; if it is set to -1, then calculate all channels/bundles
endif else begin
  j=0
  for i=0,n_elements(bundles)-1 do begin
    idx = where(strcmp(C_bundID,bundles[i],/fold_case),cnt)
    if cnt ne 0 then begin
      if j eq 0 then channels=idx else channels=[channels,idx]
      j+=1
    endif else begin
      print, ''
      print, format='("    WARNING: ",A," is not recognised as a fibre bundle! It will be ignored!")',$
             '"'+bundles[i]+'"'
      tmpstr='                                                                        '
    endelse
  endfor
endelse
print, format='($,A)',tmpstr

nslice   = input.integration.nslice	; number of 'slices' that the fibre image is divided in

an_C     = input.integration.an_C	; number of collection lens annuli 
seg_C    = input.integration.seg_C	; increment in number of collection lens segments/annulus
h_div    = input.integration.nhor_pini	; number of horizontal beamlets used to calculate the divergence
v_div    = input.integration.nvert_pini	; number of vertical beamlets used to calculate the divergence

dlambda  = input.integration.dlambda	; Size of one wavelength bin (Angstrom)
lambdatol= input.integration.lambdatol	; Wavelength tolerance (Angstrom)

npsi     = input.integration.npsi	; the number of normalised flux coordinates

; Setting for just a calculation of the intensity, not the whole spectrum
int_only = input.integration.only_intensity



;**********************************************************************************
;* IDL-related initialisation                                                     *
;**********************************************************************************

;----------------------------------------------------------------------------------
; Convert all angles from degrees to radians
;----------------------------------------------------------------------------------
B_xi     = !dtor*B_xi
B_delta  = !dtor*B_delta
B_div    = !dtor*B_div

;----------------------------------------------------------------------------------
; Convert all (non-wavelength) lengths to meter
;----------------------------------------------------------------------------------
C_fibre     = C_fibre*1e-6		; micrometer converted to meter
C_efl       = C_efl*1e-3		; mm converted to meter
C_w         = C_w*0.5e-3		; diameter in mm converted to radius in meter

;**********************************************************************************
;* Initialisation of the geometry                                                 *
;**********************************************************************************

;----------------------------------------------------------------------------------
; Calculate Beam vector
;----------------------------------------------------------------------------------
B_vec = double([cos(B_xi)*sin(B_delta), sin(B_xi)*sin(B_delta), cos(B_delta)])


;----------------------------------------------------------------------------------
; We also need to know the coordinates of the beam source
;----------------------------------------------------------------------------------
B_src = B_xyz - B_dS*B_vec

;----------------------------------------------------------------------------------
; Get the fibre coordinates for each bundle
;----------------------------------------------------------------------------------
; the names of all fields in the collset.coll structure
names = tag_names(collset.coll)

; the first bundle
; find the field with the l-coordinates
idx = where(strmatch(names,C_bundID[0]+'l', /fold_case))
; some error handling if the field doesn't exist
if idx eq -1 then begin
print, ''
errmsg = ' ERROR! Could not find the field "'+C_bundID[0]+'l" describing the l-coordinates of bundle "'$
		 + C_bundID[0] + '" in the collection optics setting file "'+collfile+'" !'
print, errmsg
print, ''
return
endif
; store the l-coordinates in the structure and convert them from mm to m
C_bund_l = create_struct(C_bundID[0],collset.coll.(idx)*1e-3)
; number of 'grid points' (i.e. fibres times slices) in this 1st channel
gp_max   = n_elements(C_bund_l.(0)) * nslice


; find the field with the m-coordinates
idx = where(strmatch(names,C_bundID[0]+'m', /fold_case))
; some error handling if the field doesn't exist
if idx eq -1 then begin
print, ''
errmsg = ' ERROR! Could not find the field "'+C_bundID[0]+'m" describing the m-coordinates of bundle "'$
		 + C_bundID[0] + '" in the collection optics setting file "'+collfile+'" !'
print, errmsg
print, ''
return
endif
; store the m-coordinates in the structure and convert them from mm to m
C_bund_m = create_struct(C_bundID[0],collset.coll.(idx)*1e-3)

; Loop through the rest of the bundles
if nchan gt 1 then begin
  for i = 1, nchan-1 do begin
    ; find the field with the l-coordinates
    idx = where(strmatch(names,C_bundID[i]+'l', /fold_case))
    ; some error handling if the field doesn't exist
    if idx eq -1 then begin
      print, ''
      errmsg = ' ERROR! Could not find the field "'+C_bundID[i]+'l" describing the l-coordinates of bundle "'$
               + C_bundID[i] + '" in the collection optics setting file "'+collfile+'" !'
      print, errmsg
      print, ''
      return
    endif
    ; store the l-coordinates in the structure and convert them from mm to m
    C_bund_l = create_struct(C_bund_l, C_bundID[i],collset.coll.(idx)*1e-3)
    ; if the number of fibres in this channel is larger than previously, increase the maximum number of fibres per channel
    if n_elements(C_bund_l.(0)) * nslice gt gp_max then gp_max   = n_elements(C_bund_l.(0)) * nslice

    ; find the field with the m-coordinates
    idx = where(strmatch(names,C_bundID[i]+'m', /fold_case))
    ; some error handling if the field doesn't exist
    if idx eq -1 then begin
      print, ''
      errmsg = ' ERROR! Could not find the field "'+C_bundID[i]+'m" describing the m-coordinates of bundle "'$
               + C_bundID[i] + '" in the collection optics setting file "'+collfile+'" !'
      print, errmsg
      print, ''
      return
    endif
    ; store the m-coordinates in the structure and convert them from mm to m
    C_bund_m = create_struct(C_bund_m, C_bundID[i],collset.coll.(idx)*1e-3)
  endfor
endif


;----------------------------------------------------------------------------------
; Calculate the emission positions on the beam and the emission vectors
;----------------------------------------------------------------------------------
; the total number of channels ('nchan') is not necessary the number of channels 
; we want to calculate (listed by 'channels'). So we rename nchan into totalchan
; and then replace nchan by the actual numbers of channels to calculate
totalchan = nchan
nchan=n_elements(channels)

xyz0=dblarr(3,nchan)			; the coordinates of the central points where the beam emission originates
vec0=dblarr(3,nchan)			; the (unit) vector in the viewing direction (from the central points)

; now lets get these emission points and vectors for each channel to calculate
; first we need to know the coordinates of the machine x- and y-vectors in
; the (k,l,m)-collection lens coordinate system:
X1   = coordtrans([1.,0.,0.],[[0,0,0],[C_k],[C_l]])
Y1   = coordtrans([0.,1.,0.],[[0,0,0],[C_k],[C_l]])
; now go through each channel
for i=0,nchan-1 do begin
  ; get the k-,l- and m-coordinates (in the collection lens coordinate system)
  ; of the channel centre from the fibre coordinates
  chanidx = channels[i]
  l0      = mean(C_bund_l.(chanidx))
  m0      = mean(C_bund_m.(chanidx))
  k0      = C_efl
  klm0    = [k0,l0,m0]
  klm0    = klm0/norm(klm0)
  ; this klm0 vector is the emission vector in the (k,l,m)-coordinate system of the
  ; collection lens. In the machine-coordinate system this becomes:
  vec0[*,i] = coordtrans(klm0,[[0,0,0],[X1],[Y1]])
  ; the emission point is determined by the intersection of the beam axis (determined
  ; by B_xyz and B_vec) and the line-of-sight (determined by C_xyz and vec0)
  xyztmp    = intersection(C_xyz,vec0[*,i],B_xyz,B_vec)
  xyz0[*,i] = xyztmp[*,0]	; because the 2 lines do no necessarily touch (the line-of-sight could just
				; just miss the beam axis), there are 2 'intersection' points returned: these
				; are the points on both lines (one the l.o.s. one on the NB) where the 2 lines
				; come closest to one and other. We've chosen the point on the l.o.s. as 'intersection'
endfor

;----------------------------------------------------------------------------------
; The etendue of 1 fibre
;----------------------------------------------------------------------------------
; we assume that the lens design is such that each fibre has the same viewing-opening angle
; (this is the case for the MAST MSE system). That opening angle is then given by: 
; - fib. angle = 2.0 * atan(radius optical core/effective focal length)
; at a distance X from the collection lens
; - the surface of the emission collected by 1 fibre is : [ X*tan(0.5 * fib. angle) ]^2 * pi
; - the solid angle under which the lens is seen is     : (radius lens)^2 * pi / X^2
; hence the etendue of 1 fibre is:
etendue = ( (C_fibre[0] * C_w * !pi)/(2.0*C_efl) )^2 


;----------------------------------------------------------------------------------
; Calculate the points on the collection lens that collect the emitted light
; and over which we will integrate
;----------------------------------------------------------------------------------
; We use the function calc_dp to calculate the coordinates of a number of 
; points homogeneously distributed over the lens.
cp1 = calc_dp(C_w,an_C,seg_C)
; the number of points:
n_C = n_elements(cp1[0,*])

; The coordinates of these points of given in the plane of the lens. In the (k,l,m)-coordinate
; system of the lens, with C_xyz as the origin, cp1 contains the l- and 
; m-coordinates and the k-coordinate is zero. In such a coordinate system
; the orignal x-axis, y-axis and origin are given by:
X1  = coordtrans([1,0,0],[[0,0,0],[C_k],[C_l]])
Y1  = coordtrans([0,1,0],[[0,0,0],[C_k],[C_l]])
O1  = coordtrans([0,0,0],[[C_xyz],[C_k],[C_l]])

; With O1,X1 and Y1 known we can transform the coordinates of the points on the
; collection lens to the original coordinate system:
cp_xyz = coordtrans([replicate(0.0,1,n_C),cp1],[[O1],[X1],[Y1]])


;----------------------------------------------------------------------------------
; Calculate the coordinates of the PINI-beamlets which will be used for the
; determination of the beam velocity distribution (the 'divergence')
;----------------------------------------------------------------------------------
; read in the coordinates of all beamlets
; REMARK: these are (y,z)-coordinates in the beam coordinate system with B_vec as
; as x-axis.
pini  = read_txt(pinifile)
npini = round(0.5*n_elements(pini))
pini  = reform(pini,2,npini)
if (h_div*v_div gt n_elements(pini[0,*])) then begin
  ; if h_div*v_div is greater than the number of pini-beamlets, we just use the all pini-beamlets
  n_div = n_elements(pini[0,*])
endif else begin
  ; if not, we resample the pini-beamlets to [h_div x v_div] 'virtual' beamlets

  ; find the maximum and minimum y-coordinate (i.e. in the horizontal direction)
  hmin = min(pini[0,*])
  hmax = max(pini[0,*])
  ; find the maximum and minimum z-coordinate (i.e. in the vertical direction)
  vmin = min(pini[1,*])
  vmax = max(pini[1,*])
  ; get a reduced set of pini-beamlets based on h_div and v_div
  if h_div eq 1 then begin
    hpini = 0.5*(hmax+hmin)
  endif else begin
    hpini = hmin + findgen(h_div)*(hmax-hmin)/(h_div-1)
  endelse
  if v_div eq 1 then begin
    vpini = 0.5*(vmax+vmin)
  endif else begin
    vpini = vmin + findgen(v_div)*(vmax-vmin)/(v_div-1)
  endelse
  pini = fltarr(2,h_div*v_div)
  m=0
  for i=0,h_div-1 do begin
    for j=0,v_div-1 do begin
      pini[*,m] = [hpini[i],vpini[j]]
      m++
    endfor
  endfor
  ; the number of divergence iteration is
  n_div = h_div*v_div
endelse

;**********************************************************************************
;* Initialisation of the spectrum                                                 *
;**********************************************************************************

;----------------------------------------------------------------------------------
; Calculate the transition wavelength
;----------------------------------------------------------------------------------
; first get the Rydberg constant
case B_M of
 1: B_mass = pm		; Hydrogen beam
 2: B_mass = pm+nm	; Deuterium beam
 3: B_mass = pm+2*nm	; Tritium beam
endcase
RydBrg    = RydInf/(1 + em/B_mass)
; transition goes from 3 (H-alpha), 4 (H-beta), 5 (H-gamma) to 2
n1 = 2.0
n2 = trans+2.0
; transition wavelength (in Angstrom) (Rydberg formula for Hydrogen)
lambda0 = 1./(RydBrg * (1./n1^2 - 1./n2^2)) * 1e10	; wavelength in vacuum
lambda0 = lambda0/nAir					; wavelength in air

;----------------------------------------------------------------------------------
; Get the dipole matrices for the transitions between the unperturbed
; basis states of n=3 (|3 l m>) and n=2 (|2 l m>)
;----------------------------------------------------------------------------------
transition_dipoles, Db, Dc, Dd

;----------------------------------------------------------------------------------
; Simple (just Stark) wavelength shift factor due to an electric field
;----------------------------------------------------------------------------------
ssf = 2.77e-7; [A*m/V]. Wavelength shift = ssf*Efld in [A]


;----------------------------------------------------------------------------------
; Calculate Beam velocity
;----------------------------------------------------------------------------------
B_EJ = B_E * e			; convert beam energy from eV to J
B_v0 = sqrt(2*B_EJ/B_mass)	; beam velocity in m/s (beam energy fluctuations ignored)

;----------------------------------------------------------------------------------
; Calculate the Doppler shifts for the central points of each channel
;----------------------------------------------------------------------------------
Dshift0 =fltarr(nchan)	; the expected Doppler shift for each spectrum
for k=0,nchan-1 do begin
  Dshift0[k]    = -(lambda0/c)*dotp(B_v0*B_vec, vec0[*,k])
endfor
Dmax = max(Dshift0)
Dmin = min(Dshift0)

;----------------------------------------------------------------------------------
; Calculate the maximum Stark shift and the pol. angle for the central 
; points of each channel
;----------------------------------------------------------------------------------
; get the Bfield at the central emission points
if strcmp('none',equifile,/fold_case) then begin	; if there is no equilibrium file than calculate
  equi = calc_equi(xyz0,$				; it using the build-in equilibrium model
                   R0, a, shafr, elong,$
                   Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
endif else begin
  equi  = read_equi(xyz0, equifile)			; else: get the equilibrium from the file
endelse
Bfld0 = equi.Bfld
psi0  = equi.psi
; get the electrostatic field at the central emission points
if strcmp('none',Erfile,/fold_case) then begin		; if there is no Er-file than calculate
  Er0 = calc_Er(xyz0,Bfld0,psi0,Ermax,Eridx)		; Er using the build-in model
endif else begin
  Er0 = read_Er(xyz0, Erfile,$				; else: read Er from file
                equifile,R0,a,shafr,elong,$
                Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
endelse

; calculated the maximum Stark shift and the pol. angle
Efld0   = fltarr(3,nchan)
Sshift0 = fltarr(nchan)
alpha0  = fltarr(nchan)
for k=0,nchan-1 do begin
  Efld0[*,k] = crossp(B_v0*B_vec,Bfld0[*,k]) + Er0[*,k]	; the E-field = vxB + Er (Lorentz + electrostatic)

  Sshift0[k] = 4 * ssf * norm(Efld0[*,k])               ; we calculate the (pure) Stark shift for 4th pi-line

  Efld1     = coordtrans(Efld0[*,k],[[0,0,0],[vec0[*,k]]])	; By expressing the E-field in the coordinate
  alpha0[k] = atan(Efld1[2],Efld1[1])		                ; system defined by the emission direction, we
  if (alpha0[k] lt 0.0) then alpha0[k] += !pi 	                ; find the pol. angle between 0 and pi
endfor
Smax = max(Sshift0)


;----------------------------------------------------------------------------------
; Get the wavelength basis on which to work
;----------------------------------------------------------------------------------
; we will need a large 'overall' wavelength basis and a reduced wavelength basis
; (to increase the speed of the calculation) for each channel:
;   - centred around the Doppler shift for that channel
;   - width determined by the Stark shift for that channel

; overall wavelength basis:
lambdamin  = Dmin - Smax - lambdatol			; minimum wavelength
lambdamax  = Dmax + Smax + lambdatol			; maximum wavelength
nlambda    = ceil((lambdamax-lambdamin)/dlambda)	; number of wavelength bins
lambda     = lambdamin + findgen(nlambda)*dlambda	; wavelength vector

; wavelength basis for each channel
lambdabase = intarr(2,nchan)
for k=0,nchan-1 do begin
  ; the wavelength basis is given as a 'start' and 'stop'-index in the overall wavelength vector
  minidx   = round(((Dshift0[k] - Sshift0[k] - lambdatol) - lambdamin)/dlambda)
  maxidx   = round(((Dshift0[k] + Sshift0[k] + lambdatol) - lambdamin)/dlambda)
  ; roundoff errrors can lead to mindidx<0 of maxidx>nlambda-1. catch this error:
  if (minidx lt 0) then minidx=0
  if (maxidx ge nlambda) then maxidx=nlambda-1
  lambdabase[*,k] = [minidx,maxidx]
endfor


;----------------------------------------------------------------------------------
; Get the psi basis (i.e. the normalised flux coordinate), R basis and Z basis
;----------------------------------------------------------------------------------
dpsi = 1.0/(npsi-1)		; step is psi
psi  = findgen(npsi)*dpsi	; psi vector

; get a minimum and maximum R
Rmin  = min(sqrt(xyz0[0,*]^2 + xyz0[1,*]^2))
Rmax  = max(sqrt(xyz0[0,*]^2 + xyz0[1,*]^2))
; Increase the R-range by 40%
deltaR = Rmax-Rmin
if deltaR lt 0.01 then deltaR = 0.1	; when only one channel is present,Rmin=Rmax! 
Rmin   = Rmin - 0.2*deltaR
Rmax   = Rmax + 0.2*deltaR
deltaR = Rmax-Rmin

; Get a minimum and maximum Z based on the max. (vertical) m-coordinate of the fibres on the focal plane (+ 20 fibre heights to get a margin),
; the distance from the collection lens and the angle between the z-vector and the m-vector of the lens
; coordinate system
dist  = norm(C_xyz-xyz0[*,0])
f_max = 0.0
for k=0,nchan-1 do begin
  if max(C_bund_m.(k)) gt f_max then fmax =  max(C_bund_m.(k))
endfor
Zmax  = 1.2*dist*(f_max + 20.0*C_fibre[1])/C_efl / abs(dotp([0.,0.,1.],C_m)) 
Zmin  = -Zmax

; Get the number R and Z bins. For nR a good number is 2*npsi. nZ is chosen such
; that the step in Z equals the step in R
nR    = 2*npsi				; number of R coordinates
dR    = (Rmax-Rmin)/float(nR-1)		; step in R
dZ    = dR				; approximate step in Z
nZ    = round((Zmax-Zmin)/dZ + 1)	; number of Z coordinates
dZ    = (Zmax-Zmin)/float(nZ-1)		; exact step in Z

; the R and Z vectors are:
R = Rmin + findgen(nR)*dR
Z = Zmin + findgen(nZ)*dZ

; if we want to 'project' psi-coordinates onto the major radius R in the midplane,
; then we need to now the psi-coordinates for R
xyztmp = [transpose(R),replicate(0.0,1,nR),replicate(0.0,1,nR)];
if equifile eq 'none' then begin		; if there is no equilibrium file than calculate
  equi = calc_equi(xyztmp,$			; used the build-in equilibrium model
                   R0, a, shafr, elong,$
                   Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
endif else begin
  equi  = read_equi(xyztmp, equifile)		; else: get the equilibrium from the file
endelse
Rpsi = equi.psi
Rm   = equi.Rm
; there are 2 possibilities: theres a set of Rpsi's for R>=Rm <=> idxpos
; and a set of Rpsi's for R<Rm <=> idxneg
idxpos = where(R ge Rm)
idxneg = where(R lt Rm)

; finally, if we want to overlay the fluxsurfaces over the RZ-emission plot,
; we need to now psi as a function of (R,Z)
xyztmp = fltarr(3,long(nR)*long(nZ))
m=long(0)
for i=0,nR-1 do begin
  for j=0,nZ-1 do begin
    xyztmp[*,m] = [R[i],0.0,Z[j]]
    m++
  endfor
endfor
if equifile eq 'none' then begin		; if there is no equilibrium file than calculate
  equi = calc_equi(xyztmp,$			; used the build-in equilibrium model1
                   R0, a, shafr, elong,$
                   Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
endif else begin
  equi  = read_equi(xyztmp, equifile)		; else: get the equilibrium from the file
endelse
RZpsi = equi.psi
RZpsi = transpose(reform(RZpsi,nZ,nR))


;----------------------------------------------------------------------------------
; Calculate the Beam fluctuations
;----------------------------------------------------------------------------------
; find the energy fluctuation level, depending on the type of modulation
case B_tEf of  
  0: begin				; no modulation
       B_nEf = 1
       B_v   = fltarr(1)
       B_Ef  = [0]			; modulation voltages with respect to the nominal voltage
       B_wf  = 1.0			; weight of the modulation voltages
     end
  1: begin				; square modulation
       B_nEf = 2
       B_Ef  = [-B_aEf/2.,B_aEf/2.]	; modulation voltages with respect to the nominal voltage
       B_wf  = [0.5      ,     0.5]	; weight of the modulation voltages
     end
  2: begin				; triangular modulation
       B_Ef  = B_aEf * (-1./2. + findgen(B_nef)/float(B_nEf-1))	; modulation voltages with respect to the nominal voltage
       B_wf  = 1.0/float(B_nEf) + fltarr(B_nef)			; weight of the modulation voltages
     end
  3: begin				; sine modulation
       B_Ef  = B_aEf/2.*cos( (findgen(B_nEf)*!pi)/(float(B_nEf-1)) )	; modulation voltages with respect to the nominal voltage
       B_wf  = 1.0/float(B_nEf) + fltarr(B_nef)			; weight of the modulation voltages
     end
  4: begin				; Histogram mode for arbitrary beam ripple (read from file)
       data  = read_txt(B_fEf)
       data  = reform(data, 2, n_elements(data)/2.)
       B_nEf = n_elements(data[0,*])	; number of sampling points corresponds with the number of voltages in the file
       B_Ef  = transpose(data[0,*])	; fist column contains the modulation/ripple voltages
       B_wf  = transpose(data[1,*])	; second columns contains their probability (should sum up to 1,
       B_wf  = B_wf/total(B_wf)		; but just in case it doesn't, we normalise it anyway)
     end
endcase
; the total beam energy and the particle velocity:
B_v  = fltarr(B_nEf)
B_Et = fltarr(B_nEf)
B_EJ = fltarr(B_nEf)
B_Et = B_E + B_Ef
B_EJ = B_Et * e			; convert beam energy from eV to J
B_v  = sqrt(2.*B_EJ/B_mass)	; the particle velocity in m/s

; the beam fluctuations give us a Doppler-shift-factor
Dshift_factor = (lambda0*B_v/c)

;----------------------------------------------------------------------------------
; Load the correctionfactor for the stokesvector. This factor takes
; the non-statistical distrubution of the energy levels into account
;----------------------------------------------------------------------------------
; Get the electron density
atomdata_temp = input.spectrum.atomdata  ; Name of the file containing the atomic data.
atomdata_temp = strmid(atomdata_temp, 15, 4) 
endens_tmp = float(atomdata_temp)

; Structure nonstat_data: first row electron density, second raw = sigma0, third raw sigma +/-1, fourth raw is pi+/-2 …
nonstat_data = fltarr(8,6)
OpenR, lun, cwd+'/' + sep + 'equi' + sep+'physics'+sep+'stark'+sep+'stark_plasma_nonstat.csv', /Get_Lun
ReadF, lun, nonstat_data
Free_Lun, lun

; Interpolate: 
nonstat = fltarr(9,4)
for i=0,3 do begin  
nonstat[i,*] = (INTERPOL(nonstat_data[*,5-i], nonstat_data[*,0],endens_tmp)) * make_array(4,1,VALUE=1b)
nonstat[8-i,*] = (INTERPOL(nonstat_data[*,5-i], nonstat_data[*,0],endens_tmp)) * make_array(4,1,VALUE=1b)
end
nonstat[4,*] = (INTERPOL(nonstat_data[*,1], nonstat_data[*,0],endens_tmp)) * make_array(4,1,VALUE=1b) ; sigma 0

;**********************************************************************************
;* declaring the arrays that will contains the final results                      *
;**********************************************************************************
; the data for each grid point:
gp_xyz   = fltarr(3,gp_max,nchan)	; coordinates of each grid point
gp_psi   = fltarr(gp_max,nchan)	; normalised flux coordinate of each grid point
gp_Bfld  = fltarr(3,gp_max,nchan)	; B-field at each grid point
gp_vel   = fltarr(3,gp_max,nchan)	; velocity (unit) vector at each grid point
gp_Efld  = fltarr(3,gp_max,nchan)	; E-field at each grid point
gp_vec   = fltarr(3,gp_max,nchan)	; the emission (unit) vector at each grid point
gp_alpha = fltarr(gp_max,nchan)	; the polarisation angle at each grid point
gp_emis  = fltarr(gp_max,nchan)	; the emission intensity at each grid point

; the emission intensity as function of psi, (R,Z) and R for each channel
psi_emis = fltarr(npsi,nchan)	; array containing the emitted intensity at each psi
RZ_emis  = fltarr(nR,nZ,nchan)	; array containing the emitted intensity at each (R,Z)
R_emis   = fltarr(nR,nchan)	; array containing the emitted intensity at each R

; the spectrum data
pstokes   = fltarr(nlambda,4,nchan)	; pi Stokes vectors (for all channels)
sstokes   = fltarr(nlambda,4,nchan)	; sigma Stokes vectors (for all channels)
stokes    = fltarr(nlambda,4,nchan)	; total Stokes vectors (for all channels)
cwlstokes = fltarr(5,nchan)		; Stokes vectors at the central wavelength (for all channels)
sigstokes = fltarr(5,nchan)             ; Stokes vectors at max. SN for sigma (for all channels)
pirstokes = fltarr(5,nchan)             ; Stokes vectors at max. SN for red shifted pi (for all channels)
pibstokes = fltarr(5,nchan)             ; Stokes vectors at max. SN for blue shifted pi (for all channels)

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen (initialisation done)
;----------------------------------------------------------------------------------
print, "done!"


;**********************************************************************************
;* Now the main loop begins                                                       *
;**********************************************************************************

;----------------------------------------------------------------------------------
; loop through the nchan positions on the beam for which a spectrum is calculated
;----------------------------------------------------------------------------------
for k=0,nchan-1 do begin
  chanidx = channels[k]
  ; get the start time
  tstart = systime(1)

  ;----------------------------------------------------------------------------------
  ; print some info on what's happening to the screen (begin calculation spectrum k)
  ;----------------------------------------------------------------------------------
  print,''
  print,FORMAT='("* Calculating spectrum for channel ",A, " at R=",(F5.2),"m :")',$
                  C_bundID[chanidx], norm(xyz0[0:1,k])

  ;----------------------------------------------------------------------------------
  ; Calculate the coordinates of the grid points and the emission vector
  ;----------------------------------------------------------------------------------
  print,FORMAT='($,A)',"  - calculating the grid point coordinates ............................ "
  ; first find the angle between the beam vector and the main emission vector
  gamma  = acos(abs(dotp(B_vec,vec0[*,k])))   ; where 0 <= gamma <= pi/2
  ; the total beam sampling width (2*B_w) then corresponds with a distance along the line-of-sight of:
  vec0_w = 2.*B_w/sin(gamma)
  ; which is chopped up in nslice slices
  ds     = vec0_w/nslice

  ; calculate the gridpoints and their emission vectors around the central point xyz0, based on the
  ; coordinates of the fibres in the fibre bundle, the lens settings, the neutral beam width, ...
  gp_calc = calc_gp(C_bund_l.(chanidx), C_bund_m.(chanidx), C_efl, C_xyz, C_k, C_l,xyz0[*,k], ds, nslice)
  gp_n    = n_elements(gp_calc[0,*,0])
  gp_xyz[*,0:gp_n-1,k] = gp_calc[*,*,0]
  gp_vec[*,0:gp_n-1,k] = gp_calc[*,*,1]
  print,FORMAT='("done!")'


  ;----------------------------------------------------------------------------------
  ; Get the B-field at the position of the grid points and 
  ; the flux coordinate psi (i.e. the flux surface it is on) 
  ;----------------------------------------------------------------------------------
  print,FORMAT='($,A)',"  - calculating the magnetic field vector at the grid points .......... "
  if strcmp('none',equifile,/fold_case) then begin	; if there is no EFIT-file then calculate our own
    equi = calc_equi(gp_xyz[*,0:gp_n-1,k],$			; equilibrium based on the input parameters
                     R0, a, shafr, elong,$
                     Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
  endif else begin
    equi  = read_equi(gp_xyz[*,0:gp_n-1,k], equifile)		; else: get the equilibrium from the EFIT-file
  endelse
  ; psi of the grid point
  gp_psi[0:gp_n-1,k]    = equi.psi
  ; the B-field at position gp is:
  gp_Bfld[*,0:gp_n-1,k] = equi.Bfld

  print,FORMAT='("done!")'

  ;----------------------------------------------------------------------------------
  ; Get the Er-field at the position of the grid points 
  ;----------------------------------------------------------------------------------
  print,FORMAT='($,A)',"  - calculating the electrostatic field vector at the grid points ..... "
  if strcmp('none',Erfile,/fold_case) then begin	; if there is no file with Er(R)-data
    Er = calc_Er(gp_xyz[*,0:gp_n-1,k],$			; than calculate Er on the input parameters
                 gp_Bfld[*,0:gp_n-1,k],gp_psi[0:gp_n-1,k],Ermax,Eridx)
  endif else begin
    Er = read_Er(gp_xyz[*,0:gp_n-1,k], Erfile,$		; else: get Er from the Er-file
                 equifile,R0,a,shafr,elong,$
                 Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
  endelse
  print,FORMAT='("done!")'


  ;----------------------------------------------------------------------------------
  ; The H-emission depends on the beam density, the electron density and the emission
  ; rate at the position of gp. For this we need a beam model. This can either be
  ; obtained from a file, or we can used our own modest beam model
  ;----------------------------------------------------------------------------------
  print, FORMAT='($,A)',"  - calculating the neutral beam parameters at the grid points ........ "
  if strcmp('none',beamfile,/fold_case) then begin	; if there is no BEAM-file than calculate our own beam model based
    beam = calc_beam(gp_xyz[*,0:gp_n-1,k],B_xyz,$		; on the input parameters.
                     B_vec,B_div,B_w/2,$
                     Bdens0,edens0,Qion,Qemit,$
                     equifile,R0, a, shafr, elong,$
                     Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)

  endif else begin
    beam = read_beam(gp_xyz[*,0:gp_n-1,k],$		; else: get the equilibrium from the BEAM-file
                     B_vec,$
                     beamfile)
  endelse
  gp_emis[0:gp_n-1,k] = beam.emission/(4.0*!pi) ; the 4*pi converts the total emission/volume element
                                                ; to the emission/volume element/solid angle
  ; the collected emission is the emission/volume element/solid angle
  ; times the length of the volume element along the line-of-sight - ds -
  ; times the etendue of 1 fibre - etendue - (Because etendue is preserved the etendue of the
  ; fibre is equal to the etendue of the volume element = solid angle under which the collection
  ; lens is seen times the base surface of the volume element. Because both the solid angle
  ; and the base surface are different for each volume element it is easier to work with the
  ; "fibre" etendue).
  gp_emis[0:gp_n-1,k] *= ds * etendue

  ; we also need to know the velocity distribution of the beam particles at each point
  beamdir = calc_beamdir(gp_xyz[*,0:gp_n-1,k], B_src, B_vec, B_hf, B_vf, B_div, pini)
  gp_vel[*,0:gp_n-1,k] = beamdir.avvel
  n_div   = n_elements(beamdir[0].vel[0,*])
  print, 'done!'


  ;----------------------------------------------------------------------------------
  ; With the B-field, the Er-field and the velocity vector and the emission vector known, we can
  ; now calculate the E-field and polarisation angle for each grid point
  ;----------------------------------------------------------------------------------
  for gp_c =0l,gp_n-1 do begin
    ; the E-field:
    gp_Efld[*,gp_c,k] = crossp(B_v0*gp_vel[*,gp_c,k],gp_Bfld[*,gp_c,k]) + Er[*,gp_c]

    ; By expressing the E-field in the coordinate system defined by the emission direction,
    ; we find the pol. angle between 0 and pi
    Efld1 = coordtrans(gp_Efld[*,gp_c,k], [[0,0,0],[gp_vec[*,gp_c,k]]])
    gp_alpha[gp_c,k] = atan(Efld1[2],Efld1[1])
    if (gp_alpha[gp_c,k] lt 0.0) then gp_alpha[gp_c,k] += !pi
  endfor



  ;----------------------------------------------------------------------------------
  ; Initialise the matrices in which the intensity will be binned as function of
  ; wavelength and pol. angle
  ;----------------------------------------------------------------------------------
  ; to increase the speed, we will not work on the complete wavelength basis, but on
  ; one centered around the expected Doppler shift Dshift0[k].
  nlambdabase = lambdabase[1,k] - lambdabase[0,k] + 1

  ; now initialise the pi and sigma Stokes vector arrays (for this channel)
  pstokes0 = fltarr(nlambdabase,4)	; pi-Stokes vectors as function of wavelength
  sstokes0 = fltarr(nlambdabase,4)	; sigma-Stokes vectors as function of wavelength

  ;----------------------------------------------------------------------------------
  ; Loop trough the grid points
  ;----------------------------------------------------------------------------------
  ; inform the user how many iterations he can expect
  ntot  = long(gp_n)*long(B_nEf)*long(n_div)*long(n_C)	; total number of iterations needed for 1 spectrum
  print,FORMAT='($,A,(I8),A)',"  - calculating emission from each grid point (",ntot," iterations) "
  ; set 'percentage done' to zero
  pctold = 0
  for gp_c = 0l,gp_n-1 do begin
    ;----------------------------------------------------------------------------------
    ; show the progress
    ;----------------------------------------------------------------------------------
    pct = round(50*float(gp_c)/gp_n)
    if pct gt pctold then begin
      print,FORMAT='($,A)','.'
      pctold = pct
    endif

    ;----------------------------------------------------------------------------------
    ; Bin the total emission intensity of the gird points to the emission matrix as 
    ; function of (R,Z), the emission vector as a function of psi,  and the
    ; emission vector as a function of R in the midplane
    ;----------------------------------------------------------------------------------
    ; emission as function of (R,Z)
    Ridx = round((norm(gp_xyz[0:1,gp_c,k])-Rmin)/dR) >0 <(n_elements(R)-1)
    Zidx = round((gp_xyz[2,gp_c,k]-Zmin)/dZ) >0 <(n_elements(Z)-1)
    if   (Ridx ge 0.0) && (Ridx lt nR)$
      && (Zidx ge 0.0) && (Zidx lt nZ) then begin
      RZ_emis[Ridx,Zidx,k] += gp_emis[gp_c,k] / (dR*dZ)
    endif

    ; emission as function of psi
    psi_idx  = round(gp_psi[gp_c,k]/dpsi)
    if (psi_idx lt npsi) then psi_emis[psi_idx,k] += gp_emis[gp_c,k] / dpsi

    ; emission as function of R in the midplane
    ; option 1: get the midplane R that has the same value for psi (i.e. projection using 'positive' and 'negative' Rpsi)
;    if (norm(gp_xyz[0:1,gp_c,k]) ge Rm) then begin		; R falls in the 'positive' Rpsi bins in the midplane.
;      mntmp = min(abs(gp_psi[gp_c,k] - Rpsi[idxpos]), Rpsi_idx)	; In the total R-vector the 'positive' R's start just
;      if ~(idxneg[0] eq -1) then begin                          ; after the 'negative' ones, so we need to add the
;        Rpsi_idx += n_elements(idxneg)				; number of 'negative' R's to the R-index.
;      endif
;    endif else begin						; R falls in the 'negative' Rpsi bins
;      mntmp = min(abs(gp_psi[gp_c,k] - Rpsi[idxneg]), Rpsi_idx)	; in the midplane.
;    endelse
;    R_emis[Rpsi_idx,k] += gp_emis[gp_c,k] / dR
    ; option 2: vertical projection
    Ridx = round((norm(gp_xyz[0:1,gp_c,k])-Rmin)/dR)
    if (Ridx ge 0) && (Ridx lt nR) then R_emis[Ridx,k] += gp_emis[gp_c,k] / dR 


    ;----------------------------------------------------------------------------------
    ; if only the intensity as a function of psi is requested, then we immediately go 
    ; the next iteration and skip over the rest of the integration loop.
    ;----------------------------------------------------------------------------------
    if int_only then continue
    ;----------------------------------------------------------------------------------


    ;----------------------------------------------------------------------------------
    ; Paschen-Back strength gamma (depends on Bfld only).
    ;----------------------------------------------------------------------------------
    ; calc. gamma
    ;----------------------------------------------------------------------------------
    Bfld  = norm(gp_Bfld[*,gp_c,k])     ; Strength B-field at the central points 
                                        ; for each channel and gridpoint
    
    b_b  = gp_Bfld[*,gp_c,k]/Bfld       ; unit vector bb (i.e. bb-axis of (bb,cc,dd)-transition coordinate system)
    
    gamma = ((e*h)/(4*!pi*em))*Bfld     ; Strength of the Zeeman/Paschen-Back part


    ;----------------------------------------------------------------------------------
    ; Loop through the beam divergence
    ;----------------------------------------------------------------------------------
    for i=0,n_div-1 do begin
      ; the diverged velocity direction is:
      gp_div = beamdir[gp_c].vel[*,i]

      ;----------------------------------------------------------------------------------
      ; Loop through the beam energy fluctuations
      ;----------------------------------------------------------------------------------
      for l=0,B_nEf-1 do begin
        ; all information about the Stark E-field is now known:
        Efld    = crossp(B_v[l]*gp_div,gp_Bfld[*,gp_c,k]) + Er[*,gp_c]	; the E-field is v x B + Er

        ;----------------------------------------------------------------------------------
        ; Stark strength and energy levels (depend on Efld and Bfld)
        ;----------------------------------------------------------------------------------
        ; calc. eps, q0 and q1
        ;----------------------------------------------------------------------------------
        eps   = float((3.*eps0*h^2)/(em*!pi) * norm(Efld))	    ; strength of the Stark part
        q0    = sqrt(eps^2 + gamma^2)                         	; energy quanta n = 2 level splitting
        q1    = sqrt(9.*eps^2 + 4.*gamma^2)                   	; energy quanta n = 3 level splitting

        ;----------------------------------------------------------------------------------
        ; Loop through the collection lens points
        ;----------------------------------------------------------------------------------
        for m=0,n_C-1 do begin
          ; the weight factor for the emitted intensity is:
          ; the probability form the velocity distribution * the probability of the beam modulation voltage
          ; * (1/number of collection lens points)
          wf = beamdir[gp_c].p[i] * B_wf[l] /float(n_C)

          ; the emission vector goes from gp to cp:
          gp_ems = (cp_xyz[*,m] - gp_xyz[*,gp_c,k])/norm(cp_xyz[*,m] - gp_xyz[*,gp_c,k])

          ; the observed Doppler shift
          Dshift  = -Dshift_factor[l]*dotp(gp_div, gp_ems)

          ; observation coordinate system
          kk = gp_ems/norm(gp_ems)          
          ll = crossp([0,0,1],kk)
          ll = ll/norm(ll)          
          mm = crossp(kk,ll)
          mm = mm/norm(mm)    

      	  ; (b_b,c_c,d_d)-transition coordinate system
          c_c = kk - dotp(kk,b_b)*b_b
      	  c_c = c_c/norm(c_c)
      	  d_d = crossp(b_b,c_c)
      	  l_bcd = [dotp(ll,b_b), dotp(ll,c_c), dotp(ll,d_d)]
      	  m_bcd = [dotp(mm,b_b), dotp(mm,c_c), dotp(mm,d_d)]      

          ; calculation of  tau
      	  cosphi  = dotp(Efld/norm(Efld),c_c)      >(-1.) <(1.)
      	  sinphi  = dotp(Efld/norm(Efld),d_d)      >(-1.) <(1.)
      	  tau     = complex(-cosphi,sinphi)      ; direction of E in the (c,d)-plane, scalar [rad]
  
          ; run starkpascheback
          starkpaschenback, gamma, eps, tau, $
                            Db, Dc, Dd,      $
                            l_bcd, m_bcd,    $ ; input
                            dE_1, S_1          ; output

          frq_1 = dE_1/h						           ; The frequency shifts of the 9 transitions [Hz]

         lambdashift = ((frq_1/c)*1e-10)*((Dshift+lambda0)^2)		; Wavelength shift Stark-Paschenback effect and conversion of [m] to Angstrom

          ; Stokes vectors pi and sigma lines 
          pS      = S_1
          pS[3,*] = transpose([0,0,0,0]) ;=Sigma -1
          pS[4,*] = transpose([0,0,0,0]) ;=Sigma 0
          pS[5,*] = transpose([0,0,0,0]) ;=Sigma +1
          sS      = S_1
          sS[0,*] = transpose([0,0,0,0]) ;=Pi -4
          sS[1,*] = transpose([0,0,0,0]) ;=Pi -3
          sS[2,*] = transpose([0,0,0,0]) ;=Pi -2
          sS[6,*] = transpose([0,0,0,0]) ;=Pi 2
          sS[7,*] = transpose([0,0,0,0]) ;=Pi 3
          sS[8,*] = transpose([0,0,0,0]) ;=Pi 4
	  
          lidx  = round((lambdashift+Dshift-lambdamin)/dlambda)	         ;length = 9
          lidx -= lambdabase[0,k]
          gidx  = where((lidx le nlambdabase-1) and (lidx ge 0),cnt)
          if cnt gt 0 then begin
          lidx = lidx[gidx]
          sstokes0[lidx,*] += nonstat[gidx,*] * (wf/dlambda) * gp_emis[gp_c, k] * sS[gidx,*]
          pstokes0[lidx,*] += nonstat[gidx,*] * (wf/dlambda) * gp_emis[gp_c, k] * pS[gidx,*]
       	  endif


        endfor
      endfor
    endfor
  endfor
  print, 'done!'
  
 
    
  ;----------------------------------------------------------------------------------
  ; Save the Stokes vectors for this channel to the Stokes vector arrays for all
  ; channels (only if a spectral calculation was performed)
  ;----------------------------------------------------------------------------------
  if ~(int_only) then begin
    ; pi, sigma and total Stokes vectors as function of wavelength
    ;-------------------------------------------------------------
    stokes0     = pstokes0 + sstokes0
    nlambda0    = n_elements(stokes0[*,0])
    pstokes[lambdabase[0,k]:lambdabase[1,k],*,k] = pstokes0
    sstokes[lambdabase[0,k]:lambdabase[1,k],*,k] = sstokes0
    stokes[lambdabase[0,k]:lambdabase[1,k],*,k]  = stokes0

    ; the central wavelength (at the peak of the spectrum) and the corresponding Stokes vector
    ;-----------------------------------------------------------------------------------------
    ; position of the peak of the spectrum
    maxint           = max(stokes0[*,0],maxidx)
    cwlstokes[0:3,k] = reform(stokes0[maxidx,*])                     ;
    cwlstokes[4,k]   = lambda0 + lambda[lambdabase[0,k]+maxidx]      ;  

    ; find the positions of max pi-blue, max. sigma, and max.pi-red and the corresponding Stokes vectors
    ;---------------------------------------------------------------------------------------------------
    ; get SN ratio
    nonzero     = where(stokes0[*,0] gt 0.0,count)   ; where the total intensity isn't zero
    SN          = fltarr(nlambda0)
    if count ne 0 then SN[nonzero] = sqrt( (stokes0[nonzero,1]^2 + stokes0[nonzero,2]^2)/stokes0[nonzero,0] )
    ; distinguish between pi and sigma regions of the spectrum
    pmins = pstokes0[*,0] - sstokes0[*,0]
    ; check if sigma ever exceeds pi (not the case for very wide filters for example):
    if min(pmins) ge 0 then begin
      ; choose sidx to be at the center of the peak of the sigma spectrum
      mx = max(sstokes0[*,0],sidx)
      ; choose pridx and pbidx to be at the peak of the S/N
      mx = max(SN,pridx)
      pbidx=pridx
    endif else begin
      ; approximate positions of pi and sigma
      mx = min(pmins,sidx)            ; approximate position of the sigma peak
      mx = max(pmins[0:sidx], pbidx)  ; approximate position of the blue pi peak
      mx = max(pmins[sidx:*], pridx)  ; approximate position of the red pi peak
      pridx = pridx+sidx
      ; then find actual positions of the maxima in SN
      range    = 4.             ; find maximum with a +/-range number of pixels around the approx. positions
      ; SN maximum around blue shifted pi
      startidx = pbidx-range
      stopidx  = pbidx+range
      if startidx lt 0        then startidx=0
      if stopidx  ge nlambda0 then stopidx=nlambda0-1
      mx       = max(SN[startidx:stopidx], pbidx)
      pbidx    = startidx + pbidx
      ; SN maximum around sigma
      startidx = sidx-range
      stopidx  = sidx+range
      if startidx lt 0        then startidx=0
      if stopidx  ge nlambda0 then stopidx=nlambda0-1
      mx = max(SN[startidx:stopidx], sidx)
      sidx = startidx + sidx
      ; SN maximum around red shifted pi
      startidx = pridx-range
      stopidx  = pridx+range
      if startidx lt 0        then startidx=0
      if stopidx  ge nlambda0 then stopidx=nlambda0-1
      mx = max(SN[startidx:stopidx], pridx)
      pridx = startidx + pridx
    endelse
    ; sigma Stokes vector
    sigstokes[0:3,k] = reform(stokes0[sidx,*])
    sigstokes[4,k]   = lambda0 + lambda[lambdabase[0,k]+sidx]
    ; red shifted pi Stokes vector
    pirstokes[0:3,k] = reform(stokes0[pridx,*])
    pirstokes[4,k]   = lambda0 + lambda[lambdabase[0,k]+pridx]
    ; blue shifted pi Stokes vector
    pibstokes[0:3,k] = reform(stokes0[pbidx,*])
    pibstokes[4,k]   = lambda0 + lambda[lambdabase[0,k]+pbidx]    
  endif
  
  ;-------------------
  ; get the stop time
  ;-------------------
  tstop  = systime(1)
  dt     = tstop-tstart
  th     = floor(dt/3600.)
  tm     = floor((dt-3600.*th)/60.)
  ts     = dt-3600.*th-60.*tm
  print,FORMAT='("                                                    Elapsed time = ",I4.2,":",I2.2,":",I2.2)',th,tm,ts

endfor
;**********************************************************************************
;* Main loop ends here.                                                           *
;**********************************************************************************



;**********************************************************************************
;* Get the spatial resolution                                                     *
;**********************************************************************************
psi_res = fltarr(7,nchan)
R_res   = fltarr(7,nchan)
for k =0,nchan-1 do begin
  ; spatial resolution as function of psi
  mx           = max(psi_emis[*,k],maxidx)
  if mx gt 0.0 then begin
    cumemis      = total(psi_emis[*,k]/mx,/cum)
    COMidx       = value_locate(cumemis, 0.5*max(cumemis))
    loweridx     = value_locate(cumemis, 0.5*(1.0-resfrac)*max(cumemis))
    upperidx     = value_locate(cumemis, 0.5*(1.0+resfrac)*max(cumemis))
    if COMidx   lt 0    then COMidx=0
    if COMidx   ge npsi then COMidx=npsi-1
    if loweridx lt 0    then loweridx=0
    if loweridx ge npsi then loweridx=npsi-1
    if upperidx lt 0    then upperidx=0
    if upperidx ge npsi then upperidx=npsi-1
    psi_res[0,k] = psi0[k]
    psi_res[1,k] = psi[maxidx]
    psi_res[2,k] = psi[COMidx]
    psi_res[3,k] = psi[loweridx]
    psi_res[4,k] = psi[upperidx]
    psi_res[5,k] = psi[upperidx]-psi[loweridx]
    psi_res[6,k] = resfrac
  endif else begin
    psi_res[0,k] = psi0[k]
    psi_res[1,k] = psi0[k]
    psi_res[2,k] = psi0[k]
    psi_res[3,k] = psi0[k]
    psi_res[4,k] = psi0[k]
    psi_res[5,k] = 0.0
    psi_res[6,k] = resfrac
  endelse

  ; spatial resolution as function of R
  mx           = max(R_emis[*,k],maxidx)
  if mx gt 0.0 then begin
    cumemis      = total(R_emis[*,k]/mx,/cum)
    COMidx       = value_locate(cumemis, 0.5*max(cumemis))
    loweridx     = value_locate(cumemis, 0.5*(1.0-resfrac)*max(cumemis))
    upperidx     = value_locate(cumemis, 0.5*(1.0+resfrac)*max(cumemis))
    if COMidx   lt 0  then COMidx=0
    if COMidx   ge nR then COMidx=nR-1
    if loweridx lt 0  then loweridx=0
    if loweridx ge nR then loweridx=nR-1
    if upperidx lt 0  then upperidx=0
    if upperidx ge nR then upperidx=nR-1
    R_res[0,k]   = norm(xyz0[0:1,k])
    R_res[1,k]   = R[maxidx]
    R_res[2,k]   = R[COMidx]
    R_res[3,k]   = R[loweridx]
    R_res[4,k]   = R[upperidx]
    R_res[5,k]   = R[upperidx]-R[loweridx]
    R_res[6,k]   = resfrac
  endif else begin
    Rtmp         = norm(xyz0[0:1,k])
    R_res[0,k]   = Rtmp
    R_res[1,k]   = Rtmp
    R_res[2,k]   = Rtmp
    R_res[3,k]   = Rtmp
    R_res[4,k]   = Rtmp
    R_res[5,k]   = 0.0
    R_res[6,k]   = resfrac
  endelse
endfor

; convert the 'wavelength shift' into a real wavelength:
lambda += lambda0 
; make channel-ID (bundle-ID) string array for the calculated spectra
chanID = C_bundID[channels]
; set the filterflag to zero (because this is a raw spectrum calculation)
filterflag = 0

;**********************************************************************************
;* Save the result                                                                *
;**********************************************************************************
save, channels, chanID, xyz0, B_v0, B_xyz, B_w, B_vec, C_xyz, C_k,     $; central points data
      vec0, Bfld0, Efld0, Dshift0, Sshift0, alpha0, psi0, lambda0,     $
      gp_xyz, gp_vel, gp_vec, gp_Bfld, gp_Efld,                        $; grid point data
      gp_alpha, gp_psi, gp_emis,                                       $
      R, Z, psi, RZpsi, Rm, RZ_emis, R_emis, psi_emis, psi_res, R_res, $; resolution data
      lambda, pstokes, sstokes, stokes,                                $; spectral data
      cwlstokes, pirstokes, pibstokes, sigstokes,                      $
      filterflag,                                                      $; filterflag
      filename=outputfile  ; output file


end
