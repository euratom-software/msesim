from scipy.io import readsav

import matplotlib.pyplot as plt
import numpy as np

file = '/work/sgibson/msesim/equi/equi_27193_0.2.sav'
eq = readsav(file)

fluxcoord = eq['fluxcoord']
r = eq['r']
z = eq['z']
bfld = eq['bfld']
rm = eq['rm']

rr,zz = np.meshgrid(r,z)
levels=np.arange(0,1.1,0.1)

#plot the fluxcoordinates

plt.figure()
plt.contour(rr,zz, fluxcoord, levels=levels)
plt.colorbar()
plt.show()

#plot each B field component

br_lvls = np.arange(-0.6,0.6,0.05)
bz_lvls = np.arange(-2,0,0.05)
bphi_lvls = np.arange(-6,1,0.05)

plt.figure()

plt.subplot(331)
plt.title('br')
plt.contour(rr,zz,bfld[:,:,0], levels=br_lvls)
plt.colorbar()

plt.subplot(332)
plt.title('bz')
plt.contour(rr,zz,bfld[:,:,1], levels=bz_lvls)
plt.colorbar()


plt.subplot(333)
plt.title('Bphi')
plt.contour(rr,zz,bfld[:,:,2], levels=bphi_lvls)
plt.colorbar()
plt.show()