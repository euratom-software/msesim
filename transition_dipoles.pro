pro transition_dipoles, Db, Dc, Dd
;+
;   TRANSITION_DIPOLES, Db, Dc, Dd
;
;   calculates the dipole matrices for the transitions between the basis states |3 l m> of the n=3 level and
;   the basis states |2 l m> of the n=2 level.
;   The used coordinate system is (b,c,d) - rather than (x,y,z), this to avoid confusion with the 
;   machine coordinate system that is already indicated by (x,y,z). 
;   Here the b-axis will later on - with the definition of the Stark-Paschen-Back interaction -
;   be defined as parallel to the magnetic field.
;   The dipole matrices along the c- and d-axes are constructed from the left and right-handed
;   polarization matrices.
;
;   This routine was based on the theoretical work found in:
;   * Chapter 2, section 2.2.1, p.43-50 of the PhD thesis of Howard Yuh (PPPL, 2005)
;
; :Output:
;    Db       : type=complex, [9 x 4] array
;               Dipole matrix describing the component of the polarization along the b-axis
;    Dc       : type=complex, [9 x 4] array
;               Dipole matrix describing the component of the polarization along the c-axis
;    Dd       : type=complex, [9 x 4] array
;               Dipole matrix describing the component of the polarization along the d-axis
;
; :History:
;   30/08/2011 - v1.0 - Maarten De Bock (m.f.m.d.bock@tue.nl)
;                1st version
;-


  ; Non-zero radial overlap integrals
  R31_20 = 3.0648154     ; equation (2.34) on page 48 in PhD thesis of Howard Yuh (PPPL, 2005)
  R32_21 = 4.7479917     ; equation (2.35) on page 48 in PhD thesis of Howard Yuh (PPPL, 2005)
  R30_21 = 0.9384042     ; equation (2.36) on page 48 in PhD thesis of Howard Yuh (PPPL, 2005)


  ; dipole matrix for the b-component - equation (2.37) on page 48 in PhD thesis of Howard Yuh (PPPL, 2005) => there called z
  ;                  2 0 0       ,        2 1 1       ,        2 1-1       ,        2 1 0
  Db     = [[                  0.,                  0.,                  0.,  sqrt(1./3.)*R30_21],$ ; 3 0 0
            [  sqrt(1./3.)*R31_20,                  0.,                  0.,                  0.],$ ; 3 1 0
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 1 1
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 1-1
            [                  0.,                  0.,                  0., sqrt(4./15.)*R32_21],$ ; 3 2 0
            [                  0.,  sqrt(1./5.)*R32_21,                  0.,                  0.],$ ; 3 2 1
            [                  0.,                  0.,  sqrt(1./5.)*R32_21,                  0.],$ ; 3 2-1
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 2 2
            [                  0.,                  0.,                  0.,                  0.] ] ; 3 2-2


  ; dipole matrix for the left-handed polarization - equation (2.38) on page 49 in PhD thesis of Howard Yuh (PPPL, 2005)
  ;                  2 0 0       ,        2 1 1       ,        2 1-1       ,        2 1 0
  Dp     = [[                  0.,  sqrt(2./3.)*R30_21,                  0.,                  0.],$ ; 3 0 0
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 1 0
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 1 1
            [ -sqrt(2./3.)*R31_20,                  0.,                  0.,                  0.],$ ; 3 1-1
            [                  0.,-sqrt(2./15.)*R32_21,                  0.,                  0.],$ ; 3 2 0
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 2 1
            [                  0.,                  0.,                  0., -sqrt(2./5.)*R32_21],$ ; 3 2-1
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 2 2
            [                  0.,                  0., -sqrt(4./5.)*R32_21,                  0.] ] ; 3 2-2
  
  ; dipole matrix for the right-handed polarization - equation (2.39) on page 49 in PhD thesis of Howard Yuh (PPPL, 2005)
  ;                  2 0 0       ,        2 1 1       ,        2 1-1       ,        2 1 0
  Dm     = [[                  0.,                  0., -sqrt(2./3.)*R30_21,                  0.],$ ; 3 0 0
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 1 0
            [  sqrt(2./3.)*R31_20,                  0.,                  0.,                  0.],$ ; 3 1 1
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 1-1
            [                  0.,                  0., sqrt(2./15.)*R32_21,                  0.],$ ; 3 2 0
            [                  0.,                  0.,                  0.,  sqrt(2./5.)*R32_21],$ ; 3 2 1
            [                  0.,                  0.,                  0.,                  0.],$ ; 3 2-1
            [                  0.,  sqrt(4./5.)*R32_21,                  0.,                  0.],$ ; 3 2 2
            [                  0.,                  0.,                  0.,                  0.] ] ; 3 2-2

  ; converting left- and right-handed dipole matrices in x- and y- matrices,
  ; transposing all to [9x4] and make them complex
  Dc    = transpose(0.5*( Dp + Dm) * complex(1, 0))
  Dd    = transpose(0.5*( Dp - Dm) * complex(0,-1))
  Db    = transpose(Db * complex(1, 0))
end
