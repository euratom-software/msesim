import idlbridge as idl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

idl.execute('restore, "/work/sgibson/msesim/beam/beam_bes_PINI_S_18501_t0.29s_63keV.xdr", /verbose')


def load_beamfile():
    """
    Stores the output of an msesim run (stored in a .dat file) in a dictionary for use in python.
    :return: Dictionary of outputs from msesim. Retrievable using the object_names as given below.
    """

    data = {}

    key_names = ("I_SPACE", "SHOT", "WHATBEAM",
                 "NEUTRAL", "NEL_SPACE", "TEL_SPACE",
                 "XC", "YC", "ZC",
                 "BEAMTIME",
                 "NEL", "TEL", "RSHOT", "ZSHOT",
                 "CURRENTFULL", "CURRENTHALF", "CURRENTTHIRD",
                 "BEAMVOLTAGE", "ABEAM", "RATEBES")

    for key_name, object_name in zip(key_names, key_names):
        data[object_name] = idl.get(key_name)

    return data

data = load_beamfile()

ratebes_full = data['RATEBES'][0,0,:,:,:]*10**6 #as a function of neutral density, electron density, and electron temperature
tel = data['TEL']['data'][58,:,:] # time, rshot, zshot
nel = data['NEL']['data'][58,:,:] # time, rshot, zshot
neutral = data['NEUTRAL'][0,0,:,:,10] #component, x, y, z

plt.figure()
plt.plot(neutral[:,10])
plt.show()