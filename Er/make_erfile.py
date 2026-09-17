import numpy as np
import matplotlib.pyplot as plt

from pyEquilibrium.equilibrium import equilibrium

from scipy.interpolate import interp1d

#text file must have first row as R values, second row as Er values

data = np.loadtxt('/work/sgibson/msesim/Er/30166_200.dat')

data_psi = data[:,0]
data_Er = data[:,1]

eq = equilibrium(gfile='/work/sgibson/msesim/equi/MASTU_equilibrium/k25_scenario_centre_conventional_for_Sam.eqdsk')

psi_n = eq.psiN(eq.R, 0)

psi_interp = interp1d(psi_n[0,50:], eq.R[50:])

r_new = psi_interp(data_psi)

Er_data = np.array([r_new[::-1], data_Er[::-1]])

Er_file = np.savetxt('Er_30166.dat', Er_data.T, delimiter=' ')