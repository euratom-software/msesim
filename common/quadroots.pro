function quadroots, coeff
; Function returns the roots of a (real) quadratic equation (i.e. polynomial of the second order)
;
; - intput: coeff  [3xn] vector-array with in each row the coefficients [c,b,a]:  a*x^2 + b*x + c = 0
;
; -returns: [2xn] vector x:  x[0,n] = (-b + sqrt(b^2 - 4*a*c))/(2*a)
;                            x[1,n] = (-b - sqrt(b^2 - 4*a*c))/(2*a)
;
;  v1.0 mdebock, 04/06/2007
;

eps = 1e-10

a = real_part(coeff[2,*])
b = real_part(coeff[1,*])
c = real_part(coeff[0,*])

n = n_elements(a)
x =  complexarr(2,n)


D = b^2 - 4*a*c
idx = where((D ge 0) AND (abs(a) gt eps),count)
;idx = where(D ge 0,count)
if (count ne 0) then begin
  x[0,idx] = complex( (-b[idx] + sqrt(D[idx]))/(2*a[idx]), 0)
  x[1,idx] = complex( (-b[idx] - sqrt(D[idx]))/(2*a[idx]), 0)
endif

idx = where((D lt 0) AND (abs(a) gt eps),count)
;idx = where(D lt 0,count)
if (count ne 0) then begin
  x[0,idx] = complex(-b[idx]/(2*a[idx]) , sqrt(-D[idx])/(2*a[idx]) )
  x[1,idx] = complex(-b[idx]/(2*a[idx]) , -sqrt(-D[idx])/(2*a[idx]) )
endif

idx = where(abs(a) le eps,count)	; equation is not quadratic but linear!
if (count ne 0) then begin
  x[0,idx] = complex(-c[idx]/b[idx] , 0)
  x[1,idx] = complex(-c[idx]/b[idx] , 0)
endif

return, x

end