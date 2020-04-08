
from pyEquilibrium.equilibrium import equilibrium

import idlbridge as idl
import numpy as np
import matplotlib.pyplot as plt

plt.ion()


#efit = equilibrium(device='MAST', shot=27193, time=0.2)
#efit = equilibrium(gfile="/work/sgibson/msesim/equi/JET_equilibrium/87123/run_87123_EFTM_JETPPF_74/gx087123.49634")
# efit=equilibrium(device='JET', shot=87123, time=49.634, dda='EFTM', with_bfield=True, verbose=True)
# efit = equilibrium(gfile="/work/sgibson/msesim/equi/Lowbeta_12.geqdsk")

efit = equilibrium(gfile="/work/sgibson/msesim/equi/novel_core_p_and_j.eqdsk")

#Gather the relevant variables that msesim wants

R = efit.R
Z = efit.Z
nR = len(R)
nZ = len(Z)

#Turn B field components into 3D array

Bfld = np.zeros((nZ,nR,3))
Bfld[:,:,0] = efit.BR(R,Z)
Bfld[:,:,1] = efit.BZ(R,Z)
Bfld[:,:,2] = efit.Bt(R,Z)

#Get normalised magnetic flux co-ordinates and radial co-ordinate of magnetic axis

fluxcoord = efit.psiN(R,Z)
Rm = np.array([efit.axis[0]])

rr,zz = np.meshgrid(R,Z)
levels=np.arange(0,1.1,0.1)

#plot the fluxcoordinates

# plt.figure()
# plt.contourf(rr,zz, fluxcoord, levels=levels)
# plt.colorbar()
# plt.show()

#plot each B field component

br_lvls = np.arange(-0.6,0.6,0.05)
bz_lvls = np.arange(-2,0.5,0.05)
bphi_lvls = np.arange(-6,1,0.05)

# plt.figure()
#
# plt.subplot(331)
# plt.title('br')
# plt.contourf(rr,zz,Bfld[:,:,0], levels=br_lvls)
# plt.colorbar()
#
# plt.subplot(332)
# plt.title('bz')
# plt.contourf(rr,zz,Bfld[:,:,1], levels=bz_lvls)
# plt.colorbar()
#
# plt.subplot(333)
# plt.title('Bphi')
# plt.contourf(rr,zz,Bfld[:,:,2], levels=bphi_lvls)
# plt.colorbar()
#plt.show()

#Need to put our variables into idl using the idlbridge:
# idl.put('R', R)
# idl.put('Z', Z)
# idl.put('Bfld', Bfld)
# idl.put('fluxcoord', fluxcoord.T)
# idl.put('Rm', Rm)
#
# idl.execute("save, R, Z, Bfld, fluxcoord, Rm, filename='equi_MASTU_novel_core_p_and_j.sav'")
