import idlbridge as idl
import numpy as np
import matplotlib.pyplot as plt

import pyuda

client = pyuda.Client()

def find_nearest(array, value):
    """
    Find a value closest to the specified value in an array.
    :param array: Array of values
    :param value: Particular value to find the closest index to
    :return: Index of the requested value in the array.
    """

    if value < array.min() or value > array.max():
        raise IndexError("Requested value is outside the range of the data.")

    index = np.searchsorted(array, value, side="left")

    if (value - array[index]) ** 2 < (value - array[index + 1]) ** 2:
        return index
    else:
        return index + 1

pulse = 53704
run = '01'
requested_time = 0.5
#filepath = '/common/uda-scratch/sgibson/efit_runs/mastu/{}/efit_sgibson_{}/epq0{}.nc'.format(pulse, run, pulse)
filepath=pulse
runtype='epq'

R = client.get(f'/{runtype}/output/profiles2D/r', source=filepath)
Z = client.get(f'/{runtype}/output/profiles2D/z', source=filepath)
Bt = client.get(f'/{runtype}/output/profiles2D/Bphi', source=filepath)
Br = client.get(f'/{runtype}/output/profiles2D/Br', source=filepath)
Bz = client.get(f'/{runtype}/output/profiles2D/Bz', source=filepath)
psin = client.get(f'/{runtype}/output/profiles2D/psiNorm', source=filepath)
time = client.get(f'/{runtype}/time', source=filepath)
rmag = client.get(f'/{runtype}/output/globalParameters/magneticAxis/R', source=filepath)

#Gather the relevant variables that msesim wants

nR = len(R.data)
nZ = len(Z.data)

#Turn B field components into 3D array

tidx = find_nearest(time.data, requested_time)

Bfld = np.zeros((nZ,nR,3))
Bfld[:,:,0] = Br.data[tidx,:,:].T
Bfld[:,:,1] = Bz.data[tidx,:,:].T
Bfld[:,:,2] = Bt.data[tidx,:,:].T

#Get normalised magnetic flux co-ordinates and radial co-ordinate of magnetic axis
fluxcoord = psin.data[tidx,:,:].T
Rm = rmag.data[tidx]

rr,zz = np.meshgrid(R.data,Z.data)
levels=np.arange(0,1.2,0.1)

#plot the fluxcoordinates

plt.figure()
plt.contour(rr,zz, fluxcoord, levels=levels)
plt.colorbar()

#plot each B field component

br_lvls = np.arange(-0.6,0.6,0.05)
bz_lvls = np.arange(-2,0.5,0.05)
bphi_lvls = np.arange(-6,1,0.05)

plt.figure()

plt.subplot(331)
plt.title('br')
plt.contourf(rr,zz,Bfld[:,:,0], levels=br_lvls)
plt.colorbar()

plt.subplot(332)
plt.title('bz')
plt.contourf(rr,zz,Bfld[:,:,1], levels=bz_lvls)
plt.colorbar()

plt.subplot(333)
plt.title('Bphi')
plt.contourf(rr,zz,Bfld[:,:,2], levels=bphi_lvls)
plt.colorbar()
plt.show()

#Need to put our variables into idl using the idlbridge:
idl.put('R', R.data)
idl.put('Z', Z.data)
idl.put('Bfld', Bfld)
idl.put('fluxcoord', fluxcoord)
idl.put('Rm', Rm)
#
idl.execute("save, R, Z, Bfld, fluxcoord, Rm, filename='equi_mu05_53704_epq.sav'")
