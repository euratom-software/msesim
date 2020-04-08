from pyEquilibrium.equilibrium import equilibrium
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

efit = equilibrium(gfile="/work/sgibson/msesim/equi/novel_edge_p_only.eqdsk")

psi_grid = np.linspace(efit.psi_axis, efit.psi_bnd, len(efit.R))
psi_gridn = ( psi_grid - efit.psi_axis )/(efit.psi_bnd - efit.psi_axis)
ffprime = efit.ffprime(psi_grid)
pprime = efit.pprime(psi_grid)

j_phi = (1/1.256) * ((ffprime/efit.R) + (pprime*efit.R))

plt.figure()
plt.subplot(221)
plt.plot(psi_gridn, ffprime)
plt.xlabel('$\psi_{n}$')
plt.ylabel('ff prime')

plt.subplot(222)
plt.plot(psi_gridn, pprime)
plt.xlabel('$\psi_{n}$')
plt.ylabel('p prime')

plt.subplot(223)
plt.plot(psi_gridn, j_phi)
plt.xlabel('$\psi_{n}$')
plt.ylabel('$j_{\phi}$ (MA/m^2)')

plt.show()