# from equi.calculate_bfieldhdf5 import equilibrium

from pyEquilibrium.equilibrium import equilibrium

import idlbridge as idl
import numpy as np
import matplotlib.pyplot as plt


efit = equilibrium(device='MAST', shot=27193, time=0.2)

#Gather the relevant variables that msesim wants

print(efit.R)

R = efit.R
Z = efit.Z
nR = len(R)
nZ = len(Z)

#Turn B field components into 3D array

Bfld = np.zeros((3,nR,nZ))
Bfld[0,:,:] = efit.BR(R,Z)
Bfld[1,:,:] = efit.BZ(R,Z)
Bfld[2,:,:] = efit.Bt(R,Z)

#Get normalised magnetic flux co-ordinates and radial co-ordinate of magnetic axis
fluxcoord = efit.psiN(R,Z)
Rm = np.array([efit.axis[0]])

rr,zz = np.meshgrid(R,Z)

plt.figure()
plt.contour(rr,zz,fluxcoord)
plt.colorbar()
plt.show()

#Need to put our variables into idl using the idlbridge:
idl.put('R', R)
idl.put('Z', Z)
idl.put('Bfld', Bfld)
idl.put('fluxcoord', fluxcoord)
idl.put('Rm', Rm)

idl.execute("save, R, Z, Bfld, fluxcoord, Rm, filename='test_equifile.sav'")
