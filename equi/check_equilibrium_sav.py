from scipy.io import readsav

import matplotlib.pyplot as plt
import numpy as np

#file = '/work/sgibson/msesim/equi/MAST_equilibrium/equi_MAST_19522_t0.070.xdr'
#file = '/work/sgibson/msesim/equi/equi_MASTU_1MA_P4_CATIA.sav'
#file = '/work/sgibson/msesim/equi/equi_MASTU_k25_scenario_centre.sav'
file = '/work/sgibson/msesim/equi/equi_MASTU_mastlike.sav'

# file = '/work/sgibson/msesim/equi/equi_JET_87123_49.634s_eftm.sav'

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
plt.contour(rr,zz, fluxcoord.T, levels=levels)
plt.colorbar()
plt.show()

plt.figure()
plt.plot(r, bfld[65,:,1])
plt.show()

#plot each B field component

br_lvls = np.arange(-0.6,0.6,0.05)
bz_lvls = np.arange(-2,0.5,0.05)
bphi_lvls = np.arange(-6,1,0.05)

plt.figure()

plt.subplot(331)
plt.title('br')
plt.contourf(rr,zz,bfld[:,:,0], levels=br_lvls)
plt.colorbar()

plt.subplot(332)
plt.title('bz')
plt.contourf(rr,zz,bfld[:,:,1], levels=bz_lvls)
plt.colorbar()

plt.subplot(333)
plt.title('Bphi')
plt.contourf(rr,zz,bfld[:,:,2], levels=bphi_lvls)
plt.colorbar()
plt.show()