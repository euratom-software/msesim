import adas as adas
import numpy as np
import matplotlib.pyplot as plt

#Get the beam stopping coefficients for a 6% carbon impurity plasma

stopping_files=['/home/adas/adas/adf21/bms10#h/bms10#h_h1.dat',
        '/home/adas/adas/adf21/bms97#h/bms97#h_c6.dat']

fraction=[0.94, 0.06]

beam_energies=np.geomspace(60e3, 75.0e3, 100)

te=1.0e3 #1keV plasma

dens= np.geomspace(1.0e12, 1e14, 100) #cm^3

beam_stopping_coefficient = np.zeros((len(dens), len(beam_energies)))

for d, density in enumerate(dens):
    for e, energy in enumerate(beam_energies):
        beam_stopping_coefficient[d,e] = adas.read_adf21(files=stopping_files,fraction=fraction,energy=energy,te=te,dens=density)

#now get the fractional population of the n=3 level in the beam due to a 6% carbon contaminated plasma

emission_files=['/home/adas/adas/adf22/bme10#h/bme10#h_h1.dat'] #,'/home/adas/adas/adf22/bms97#h/bms97#h_c6.dat'

fraction= [1] #[0.94, 0.06]

beam_emission_coefficients = np.zeros((len(dens), len(beam_energies)))

for d, density in enumerate(dens):
    for e, energy in enumerate(beam_energies):
        beam_emission_coefficients[d,e] = adas.read_adf22(files=emission_files, fraction=fraction, energy=beam_energies[e],te=te,dens=dens[d])


plt.figure()
plt.plot(dens*10**6, beam_stopping_coefficient[:,50]*10**-6)
plt.xscale('log')
plt.xlabel('Density n$_{e}$ (m$^{-3}$')
plt.ylabel('Beam Stopping Coefficient (m$^{3}$/s)')
plt.show()

plt.figure()
plt.plot(dens*10**6, beam_emission_coefficients[:,50]*10**-6, color='C1')
plt.xscale('log')
plt.xlabel('Density n$_{e}$ (m$^{-3}$')
plt.ylabel('Beam Emission Coefficient (ph m$^{3}$/s)')
plt.show()

# plt.figure()
# plt.subplot(221)
# c = plt.pcolormesh(dens, beam_energies*10**-3, beam_stopping_coefficient*10**-6)
# plt.xlabel('n$_{e}$ (cm$^{-3}$)')
# plt.ylabel('Beam Energy (keV)')
# plt.colorbar(c)
# c.set_label('Effective Beam Stopping Coefficient (m$^{3}$/s)')
#
# plt.subplot(222)
# c2 = plt.pcolormesh(dens, beam_energies*10**-3, beam_emission_coefficients*10**-6)
# plt.xlabel('n$_{e}$ (cm$^{-3}$)')
# plt.ylabel('Beam Energy (keV)')
# plt.colorbar(c)
# c2.set_label('Beam Emission Coefficients (ph m$^{3}$/s)')
#
# plt.tight_layout()
# plt.show()
