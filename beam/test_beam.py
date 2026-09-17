import numpy as np
import matplotlib.pyplot as plt
from jqc import jqc_plot
from scipy.io import readsav

from scipy.misc import rotate

bdens = readsav('/work/sgibson/msesim/bdens.sav')['bdens3']

R0    = 0.83 #;0.88 ;0.83 for MAST
a     = 0.5 #;0.65 ;0.6 for MAST
Bphi  = -0.5
q0    = 1.0
qa    = 3.0
qidx  = 3.0
Bp0   = 0.0
Bpa   = 0.0
Bpidx = 1.0
shafr = 0.1
elong = 1.66

# beam parameters
B0   =[0,-2,0]
w0   =0.1
xi   = 70
delta= 90
chi  = 1.2

Bdens0=5e15
edens0= 3e19 #;5e19
Qion  =5e-20
Qemit =2.5e-12

Brange =[0.4,1.90]
nl     = 20
nw     = 20
ntheta = 18

# degrees to radians
d2r = np.pi/180
dl     = (Brange[1]-Brange[0])/(nl-1) #steps along beam
dw     = (3*w0)/(nw-1) #steps across beam axis
dtheta = 2*np.pi/ntheta #angle segments around beam axis

l = Brange[0] + np.arange(0,nl,1) * dl
d = np.arange(0,nw,1) * dw
d2 = np.arange(0, 2*nw-1, 1) * dw - (nw - 1) * dw
lddens = np.zeros((nl, 2 * nw - 1))
idxpi = ntheta / 2

ldens = bdens[:,0,0]

lddens[0:nl - 1, 0:nw - 1] = rotate(Bdens3[*, *, idxpi[0]], 7)
lddens[0:nl - 1, nw:2 * nw - 2] = Bdens3[*, 1: *, 0]


# surface, lddens, l, d2, color = 0, charsize = 2.0, charthick = 1.2,
# az = -30, xtitle = 'Length along the beam (m)',
# zs = 1, zr = [0, 1.1 * Bdensmax],$
# ytitle = 'Distance from beam axis (m)', ztitle = 'Beam density (m^-3)'








