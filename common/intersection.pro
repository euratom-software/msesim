
function intersection, A0, Av, B0, Bv
; C = INTERSECTION(A0, AV, B0, BV)
;
; function calculates the 'intersection' of a 2 lines in 3D 
; (A and B, defined by the origins A0/B0 and the vectors Av/Bv).
; These lines do not necessarily lie in one plane: it could be skew lines.
; Therefore, the 'intersection' points are defined as the points on lines
; (one on A, one on B) that have the closest distance to one and other.
; If the line do lie in the same plane the 2 'intersection' points overlap.
; If the lines are parallel, -1 is returned
;
; input:  A0    - [3x1] vector giving a point on line A
;         Av    - [3x1] vector describing the direction of line A
;         B0    - [3x1] vector giving a point on line B
;         Bv    - [3x1] vector describing the direction of line B
;
; output: C     - [3x2] array,
;                 the first row giving the coordinates of the intersection point on line A,
;                 the second row giving the coordinates of the intersection point on line B.
;
; v1.0, mdebock 09/07/2008
;

; points on A and B are given by:
;   As = A0 + s*Av
;   Bt = B0 + t*Bv
; with s and t real numbers.
;
; the distance squared between these points is:
;   D2 = (As-Bt).(As-Bt)
;      = As.As - 2 As.Bt + Bt.Bt
;      =   A0.A0 + 2s A0.Av + s^2 Av.Av
;        -2 (A0.B0 + t A0.Bv + s Av.B0 + st Av.Bv)
;        + B0.B0 + 2t B0.Bv + t^2 Bv.Bv
;      =   s^2 (Av.Av) + t^2 (Bv.Bv) -2st (Av.Bv)
;        + 2s (A0.Av - Av.B0) + 2t (B0.Bv - A0.Bv)
;        + (A0.A0 + B0.B0 - 2 A0.B0)
;
; the closest distance is where the partial derivatives to s and t are both zero:
;   s (Av.Av) - t (Av.Bv) + Av.(A0-B0) = 0
;   t (Bv.Bv) - s (Av.Bv) + Bv.(B0-A0) = 0
; or
;   [ Av.Av,-Av.Bv]   *   [s]    = [(B0-A0).Av]
;   [-Av.Bv, Bv.Bv]       [t]      [(A0-B0).Bv]
; or, inverting the 2x2 matrix:
;   [s] =  1/((Av.Av)(Bv.Bv)-(Av.Bv)^2) [Bv.Bv, Av.Bv] * [(B0-A0).Av]
;   [t]                                 [Av.Bv, Av.Av]   [(A0-B0).Bv]
; or
;    s = (-(Bv.Bv)(A0-B0).Av + (Av.Bv)(A0-B0).Bv)/((Av.Av)(Bv.Bv)-(Av.Bv)^2)
; and
;    t = (-(Av.Bv)(A0-B0).Av + (Av.Av)(A0-B0).Bv)/((Av.Av)(Bv.Bv)-(Av.Bv)^2)
;
; Only if N=(Av.Av)(Bv.Bv)-(Av.Bv)^2=0, there is no solution.
; Which is normal, because it means Av=Bv and the lines are parallel!
;
; now, enough math-talk, lets implement

; calculate dot products
AvAv = total(Av*Av)
BvBv = total(Bv*Bv)
AvBv = total(Av*Bv)
AB0  = A0-B0
AB0Av= total(AB0*Av)
AB0Bv= total(AB0*Bv)

; calculate denumenator
N = AvAv*BvBv - AvBv^2
if N eq 0.0 then return, -1 ; parallel lines

; calculate s and t
s = (-BvBv*AB0Av + AvBv*AB0Bv)/N
t = (-AvBv*AB0Av + AvAv*AB0Bv)/N

; and from that the points
C      = fltarr(3,2)
C[*,0] = A0+s*Av
C[*,1] = B0+t*Bv

return, C

end