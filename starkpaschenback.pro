pro starkpaschenback, gamma, eps, tau,   $ ; input : Stark-Paschen-Back terms
                      Db, Dc, Dd, l, m,  $ ; input : transition dipole matrices and l and m unit vectors (in polarisation plane) in (b,c,d) coordinates
                      dE_out, S_out        ; output: energy of the transitions and Stokes vectors
;+
;   STARKPASCHENBACK, gamma, eps, tau,   $ ; input : Stark-Paschen-Back terms
;                     Db, Dc, Dd, l, m,  $ ; input : transition dipole matrices and corresponding unit vectors
;                     dE_out, S_out        ; output: energy of the transitions and Stokes vectors
;
;   calculates the Stark-Paschen-Back spectrum. It returns the energy shift for each of the 15 Stark-Paschen-Back transitions
;   and the corresponding normalized Stokes vector. The normalization of the Stokes vectors is such that total(S[0]) = 1
;
;   REMARK: For the magnetic field, the Paschen-Back rather than the Zeeman effect is used:
;           Whereas Zeeman works on the total angular momentum J of spin-orbit coupled systems, the
;           Paschen-Back effect works on the orbital momentum J and spin momentum S seperately.
;           When the external magnetic field B is larger than 0.3T (i.e. the maximum intrinsic orbital field),
;           the spin will be decoupled from the orbital angular momentum, hence this is when the Paschen-Back effect
;           it to be used. Here we only look at the Paschen-Back effect on the orbital momentum L. The effect on the
;           spin momentum S is neglected (this causes a fine-splitting).
;           For most tokamaks the magnetic field B is larger than 0.3T, so this approach is valid.
;           For spherical tokamaks, however, and for Reversed Field Pinches (RFP) the intrisic orbital field and
;           the external field are of the same order. We consider it still OK to use the Paschen-Back description even
;           in this case, because, although the magnetic field will work on a hybric angular momentum that is a 
;           linear combination of L and J, the dominating electric field will still work on L. 
;
;   the (b,c,d) coordinate system is defined such that b is parallel to the B-field and c is in the plane spanned by
;   the B-field and the emission wavevector k:
;     b = B/norm(B)
;     c = k - dotp(k,b)*b
;     d = crossp(b,c)
;
;   This routine was based on the theoretical work found in:
;   * R.C. Isler, Physical Review A, vol.14 no.3 p.1015 (1976)
;   * Chapter 2 of the PhD thesis of Howard Yuh (PPPL, 2005)
;
; :Input:
;    gamma     : required, type=float, scalar
;               the strength of the Paschen-Back term: gamma = (e*h)/(4*!pi*me)      * Bfld    [eV]
;    eps       : required, type=float, scalar
;               the strength of the Stark term       : eps   = (3*eps0*h^2)/(me*!pi) * Efld    [eV]
;    tau       : required, type=complex, scalar
;               -exp(-i*phi) with phi the angle w.r.t. c-axis of the E-field vector [rad]
;    Db       : required, type=complex, [9 x 4] array
;               Dipole matrix describing the x-component of the polarization
;    Dc       : required, type=complex, [9 x 4] array
;               Dipole matrix describing the y-component of the polarization
;    Dd       : required, type=complex, [9 x 4] array
;               Dipole matrix describing the z-component of the polarization
;    l        : required, type=float, [3 x 1] vector
;               The l unit vector on the polarisaiton plane in (b,c,d)-coordinates
;    m         : required, type=float, [3 x 1] vector
;               The m unit vector on the polarisaiton plane in (b,c,d)-coordinates

;
; :Output:
;    frq_out   : type=float, 15-element vector
;               the frequency shifts of the 15 transitions [Hz]
;    S_out     : type=float, [15 x 4]-array
;               the normalized Stokes vectors for the 15 transitions
; :Keywords:
;    only_stark: optional, type=byte, scalar
;               if set a pure Stark calculation is done (hence no Zeeman part)
;
; :History:
;   29/07/2011 - v1.0 - Maarten De Bock (m.f.m.d.bock@tue.nl)
;                1st version
;-


  ; energy quanta
  ;--------------
  q0    = sqrt(eps^2 + gamma^2)                   ; energy quanta n=2 level splitting
  q1    = sqrt(9.*eps^2 + 4.*gamma^2)             ; energy quanta n=3 level splitting


  ; calculate hybridization matrices for the Stark-Zeeman states
  ;-------------------------------------------------------------
  ; for the n=3-level
  E3 = q1*[ 1.0,-1.0, 0.5,-0.5, 0.0, 0.0, 0.0, 0.5,-0.5]    ; energy shifts Stark-Zeeman states
  H3 = complexarr(9,9)
  ; k = 0 Stark-Zeeman state as function of ...
  H3[0,0] = 3.*sqrt(3.)*eps^2                                  /q1^2;  3 0 0
  H3[0,2] = 1.5*eps*(2.*gamma+q1)*tau                          /q1^2;  3 1 1
  H3[0,3] = 1.5*eps*(2.*gamma-q1)*conj(tau)                    /q1^2;  3 1-1
  H3[0,4] = -3.*sqrt(3.)*eps^2/(2.*sqrt(2.))                   /q1^2;  3 2 0
  H3[0,7] = 0.25*(8.*gamma^2+9.*eps^2+4.*gamma*q1)*tau^2       /q1^2;  3 2 2
  H3[0,8] = 0.25*(8.*gamma^2+9.*eps^2-4.*gamma*q1)*conj(tau)^2 /q1^2;  3 2-2
  ; k = 1 Stark-Zeeman state as function of ...
  H3[1,0] = 3.*sqrt(3.)*eps^2                                  /q1^2;  3 0 0
  H3[1,2] = 1.5*eps*(2.*gamma-q1)*tau                          /q1^2;  3 1 1
  H3[1,3] = 1.5*eps*(2.*gamma+q1)*conj(tau)                    /q1^2;  3 1-1
  H3[1,4] = -3.*sqrt(3.)*eps^2/(2.*sqrt(2.))                   /q1^2;  3 2 0
  H3[1,7] = 0.25*(8.*gamma^2+9.*eps^2-4.*gamma*q1)*tau^2       /q1^2;  3 2 2
  H3[1,8] = 0.25*(8.*gamma^2+9.*eps^2+4.*gamma*q1)*conj(tau)^2 /q1^2;  3 2-2
  ; k = 2 Stark-Zeeman state as function of ...
  H3[2,0] = 4.*sqrt(3.)*eps*gamma                              /q1^2;  3 0 0
  H3[2,2] = 0.5*(4.*gamma^2-9.*eps^2+2.*gamma*q1)*tau          /q1^2;  3 1 1
  H3[2,3] = 0.5*(4.*gamma^2-9.*eps^2-2.*gamma*q1)*conj(tau)    /q1^2;  3 1-1
  H3[2,4] = -sqrt(6.)*eps*gamma                                /q1^2;  3 2 0
  H3[2,7] = -1.5*eps*(2.*gamma+q1)*tau^2                       /q1^2;  3 2 2
  H3[2,8] = -1.5*eps*(2.*gamma-q1)*conj(tau)^2                 /q1^2;  3 2-2
  ; k = 3 Stark-Zeeman state as function of ...
  H3[3,0] = 4.*sqrt(3.)*eps*gamma                              /q1^2;  3 0 0
  H3[3,2] = 0.5*(4.*gamma^2-9.*eps^2-2.*gamma*q1)*tau          /q1^2;  3 1 1
  H3[3,3] = 0.5*(4.*gamma^2-9.*eps^2+2.*gamma*q1)*conj(tau)    /q1^2;  3 1-1
  H3[3,4] = -sqrt(6.)*eps*gamma                                /q1^2;  3 2 0
  H3[3,7] = -1.5*eps*(2.*gamma-q1)*tau^2                       /q1^2;  3 2 2
  H3[3,8] = -1.5*eps*(2.*gamma+q1)*conj(tau)^2                 /q1^2;  3 2-2
  ; k = 4 Stark-Zeeman state as function of ...
  H3[4,0] = sqrt(2.)/3.*(8.*gamma^2-9.*eps^2)                  /q1^2;  3 0 0
  H3[4,2] = -3.*sqrt(6.)*gamma*eps*tau                         /q1^2;  3 1 1
  H3[4,3] = -3.*sqrt(6.)*gamma*eps*conj(tau)                   /q1^2;  3 1-1
  H3[4,4] = -1/6.*(8.*gamma^2-9.*eps^2)                        /q1^2;  3 2 0
  H3[4,7] = 9.*sqrt(3.)/(2.*sqrt(2.)) * eps^2*tau^2            /q1^2;  3 2 2
  H3[4,8] = 9.*sqrt(3.)/(2.*sqrt(2.)) * eps^2*conj(tau)^2      /q1^2;  3 2-2
  ; k = 5 Stark-Zeeman state as function of ...
  H3[5,0] = 1./3.*q1^2                                         /q1^2;  3 0 0
  H3[5,4] = 2./3.*sqrt(2.)*q1^2                                /q1^2;  3 2 0
  ; k = 6 Stark-Zeeman state as function of ...
  H3[6,1] = 2.*gamma                                           /q1  ;  3 1 0
  H3[6,5] = -3./sqrt(2.)*eps*tau                               /q1  ;  3 2 1
  H3[6,6] = -3./sqrt(2.)*eps*conj(tau)                         /q1  ;  3 2-1
  ; k = 7 Stark-Zeeman state as function of ...
  H3[7,1] = 3./sqrt(2.)*eps                                    /q1  ;  3 1 0
  H3[7,5] = 0.5*(2.*gamma+q1)*tau                              /q1  ;  3 2 1
  H3[7,6] = 0.5*(2.*gamma-q1)*conj(tau)                        /q1  ;  3 2-1
  ; k = 8 Stark-Zeeman state as function of ...
  H3[8,1] = 3./sqrt(2.)*eps                                    /q1  ;  3 1 0
  H3[8,5] = 0.5*(2.*gamma-q1)*tau                              /q1  ;  3 2 1
  H3[8,6] = 0.5*(2.*gamma+q1)*conj(tau)                        /q1  ;  3 2-1

  ; for the n=2-level
  E2 = q0*[ 0.0, 1.0,-1.0, 0.0]                             ; energy shifts Stark-Zeeman states
  H2 = complexarr(4,4)
  ; l = 0 Stark-Zeeman state as function of ...
  H2[0,0] = gamma                                              /q0  ;  2 0 0
  H2[0,1] = -1.0/sqrt(2.)*eps*tau                              /q0  ;  2 1 1
  H2[0,2] = -1.0/sqrt(2.)*eps*conj(tau)                        /q0  ;  2 1-1
  ; l = 1 Stark-Zeeman state as function of ...
  H2[1,0] = 1.0/sqrt(2.)*eps                                   /q0  ;  2 0 0
  H2[1,1] = 0.5*(gamma+q0)*tau                                 /q0  ;  2 1 1
  H2[1,2] = 0.5*(gamma-q0)*conj(tau)                           /q0  ;  2 1-1
  ; l = 2 Stark-Zeeman state as function of ...
  H2[2,0] = 1.0/sqrt(2.)*eps                                   /q0  ;  2 0 0
  H2[2,1] = 0.5*(gamma-q0)*tau                                 /q0  ;  2 1 1
  H2[2,2] = 0.5*(gamma+q0)*conj(tau)                           /q0  ;  2 1-1
  ; l = 3 Stark-Zeeman state as function of ...
  H2[3,3] = q0                                                 /q0  ;  2 1 0

  ; calculate the dipole matrices and enery shifts for the
  ; transitions between the Stark-Zeeman states
  ;-------------------------------------------------------
  ; full dipole matrices
  fDb = conj(transpose(H2)) ## Db ## H3
  fDc = conj(transpose(H2)) ## Dc ## H3
  fDd = conj(transpose(H2)) ## Dd ## H3
  fDb = reform(fDb,36,1)
  fDc = reform(fDc,36,1)
  fDd = reform(fDd,36,1)

  ; frequency shifts
  E3  = rebin(E3,9,4)
  E2  = rebin(transpose(E2),9,4)
  dE  = reform((E3-E2),36,1)
  ; electric field vectors
  El = fDb*l[0] + fDc*l[1] + fDd*l[2]
  Em = fDb*m[0] + fDc*m[1] + fDd*m[2]

  ; stokes vectors
  S  = fltarr(36,4)
  S[*,0] = abs(El)^2 + abs(Em)^2
  S[*,1] = abs(El)^2 - abs(Em)^2
  S[*,2] = 2.*real_part(El * conj(Em))
  S[*,3] = 2.*imaginary(El * conj(Em))

;  dE_out      = dE[[10,1,12,19,3,13,21,4,11,22,2,9,20,0,18]]
  dE_out      = dE[[19,3,13,21,4,11,22,2,9]]
;  S_out       = fltarr(15,4)
  S_out       = fltarr(9,4)
;  S_out[0,*]  = S[10,*]
;  S_out[1,*]  = total(S[[1,28],*],1)
;  S_out[2,*]  = total(S[[12,17],*],1)
  S_out[0,*]  = S[19,*]
  S_out[1,*]  = total(S[[3,8,30,35],*],1)
  S_out[2,*]  = total(S[[13,14,15],*],1)
  S_out[3,*]  = total(S[[21,26],*],1)
  S_out[4,*]  = total(S[[4,5,6,31,32,33],*],1)
  S_out[5,*]  = total(S[[11,16],*],1)
  S_out[6,*]  = total(S[[22,23,24],*],1)
  S_out[7,*]  = total(S[[2,7,29,34],*],1)
  S_out[8,*]  = S[9,*]
;  S_out[12,*] = total(S[[20,25],*],1)
;  S_out[13,*] = total(S[[0,27],*],1)
;  S_out[14,*] = S[18,*]

  ; sum over degenerate transitions
;  sidx = sort(dE)
;  dE   = dE[sidx]
;  S    = S[sidx,*]
;  uidx = uniq(dE)
;  dE_out = dE[uidx]
;  S_out   = fltarr(15,4)
;  idx0    = 0
;  for i=0,n_elements(uidx)-1 do begin
;    idx1 = uidx[i]
;    S_out[i,*] = total(S[idx0:idx1,*],1)
;    idx0 = idx1+1
;  endfor

  ; normalize S
  S_out = S_out/total(S_out[*,0])
end
