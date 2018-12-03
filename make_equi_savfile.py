from equi.calculate_bfieldhdf5 import equilibrium
import idlbridge as idl
import numpy as np
from scipy.io import readsav

#Load a local HDF5 equilibrium file, and specify the timeslice

efit = equilibrium(specified_time=0.365, efit_filepath='/work/sgibson/msesim/equi/efitOut.hdf5')

#Gather the relevant variables that msesim wants
R = efit.R
Z = efit.Z
nR = len(R)
nZ = len(Z)

#Turn B field components into 3D array

Bfld = np.zeros((3,nR,nZ))
Bfld[0,:,:] = efit.B_R
Bfld[1,:,:] = efit.B_Z
Bfld[2,:,:] = efit.B_phi

#Get magnetic flux co-ordinates and radial co-ordinate of magnetic axis
fluxcoord = efit.psi_2d
Rm = np.array([efit.r_magaxis])

#Need to put our variables into idl using the idlbridge:
idl.put('R', R)
idl.put('Z', Z)
idl.put('Bfld', Bfld)
idl.put('fluxcoord', fluxcoord)
idl.put('Rm', Rm)

#Get them back, and execute the save function to save these variables in a sav file format
idl.get('R')
idl.get('Z')
idl.get('Bfld')
idl.get('fluxcoord')
idl.get('Rm')

idl.execute("save, R, Z, Bfld, fluxcoord, Rm, filename='equi_MAST_24409_t0.365s.sav'")
