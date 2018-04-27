import idlbridge as idl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from msesim_analysis.read_output import Channel, CentralPoints, GridPoints, Resolution, SpectralData
from mse_spectrometer_data.mast_spectrum import Spectrum

def get_channel_list(channel, channel_number):
    """
    Finds the index associated with the channel number we want to plot.

    :param channel: (class) channels simulated in MSESIM
    :param channel_number: (int) The channel number we want to look at
    :return: index in list associated with the specified channel number
    """
    channel_list = channel.data['channels']
    if channel_number in channel_list:
        channel_index = np.argwhere(channel_list==channel_number)[0][0]
    else:
        print('MSESIM was not run for this channel number! Please pick from the following channel numbers: {}'.format(channel.data['channels']))
    return channel_index

def total_stokes_vectors(spectral_data, channel_index):

    """

    :param spectral_data: (class) Contains all spectral data output from MSESIM
    :param channel_index: list index of desired channel number
    :return: Total stokes vector components s0, s1, s2, s3 and total linear polarised intensity (sigma + pi)
    """

    s0_total = spectral_data.data["total_stokes_vector"][channel_index,0,:]
    s1_total = spectral_data.data["total_stokes_vector"][channel_index,1,:]
    s2_total = spectral_data.data["total_stokes_vector"][channel_index,2,:]
    s3_total = spectral_data.data["total_stokes_vector"][channel_index,3,:]

    linear_polarised_intensity_total = np.sqrt(s1_total**2 + s2_total**2)

    return s0_total, s1_total, s2_total, s3_total, linear_polarised_intensity_total

def get_wavelength(spectral_data):
    wavelength = spectral_data.data["wavelength_vector"]
    return wavelength

def pi_sigma_stokes_vectors(spectral_data):

    s0_pi = spectral_data.data["pi_stokes_vector"][0,0,:]
    s1_pi = spectral_data.data["pi_stokes_vector"][0,1,:]
    s2_pi = spectral_data.data["pi_stokes_vector"][0,2,:]
    s3_pi = spectral_data.data["pi_stokes_vector"][0,3,:]

    linear_polarisation_intensity_pi = np.sqrt(s1_pi**2 + s2_pi**2)

    s0_sigma = spectral_data.data["sigma_stokes_vector"][0,0,:]
    s1_sigma = spectral_data.data["sigma_stokes_vector"][0,1,:]
    s2_sigma = spectral_data.data["sigma_stokes_vector"][0,2,:]
    s2_sigma = spectral_data.data["sigma_stokes_vector"][0,3,:]

    return

def fit(x,c):
    return c

def remove_background(spectrum):

    center = int(len(spectrum.data[4,:])/2)
    y = spectrum.data[4,center:]
    x = spectrum.wavelengths[center:]

    p0 = [500]
    params = curve_fit(fit,x,y,p0=p0)
    c = params[0]
    background = np.ones((len(spectrum.data[3,:])))*c

    data_nobackground = spectrum.data[3,:] - background

    return data_nobackground

def plot_mse_spectra(data_nobackground, spectrum, wavelength, linear_polarised_intensity_total,channel, channel_index):

    print(spectrum.wavelengths)

    plt.figure()
    plt.plot(spectrum.wavelengths*10, data_nobackground, label='MSE spectrometer data')
    plt.plot(wavelength,linear_polarised_intensity_total, label='Channel {}, fiber {}'.format(channel.data['channels'][channel_index], channel.data['channel_id'][channel_index]))
    plt.legend()
    plt.title('MSE Spectrum for t={}'.format(spectrum.times[0][4]))
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    plt.show()

    return

def beam_emission(resolution):

    emission_RZ = resolution.data["emission_intensity(RZ)"]
    R = resolution.data["R"]
    Z = resolution.data["Z"]
    return R,Z, emission_RZ

def plot_beam_emission(R, Z, emission_RZ):
    plt.figure()
    plt.pcolormesh(R, Z, emission_RZ[1, :, :])
    plt.title('Beam emission intensity in poloidal plane')
    plt.colorbar()
    plt.show()
    return


print('Loading spectrum data...')
#calib_filename = '/home/sgibson/PycharmProjects/msesim/mse_spectrometer_data/wavelength_calibration_28101.sav'
calib_filename = '/projects/diagnostics/MAST/MSE/repos_mse/codes/lib/spec_parameters_28101.sav'
spectrum = Spectrum(28101, 1, calib_filename)

channel = Channel()
central_points = CentralPoints()
grid_points = GridPoints()
resolution = Resolution()
spectral_data = SpectralData()

channel_number = 23
channel_index = get_channel_list(channel, channel_number)
print('Simulated channels {} connected to fibers {}'.format(channel.data['channels'], channel.data['channel_id']))

s0_total, s1_total, s2_total, s3_total, linear_polarised_intensity_total = total_stokes_vectors(spectral_data, channel_index)
wavelength = get_wavelength(spectral_data)

data_nobackground = remove_background(spectrum)
plot_mse_spectra(data_nobackground, spectrum, wavelength, linear_polarised_intensity_total, channel, channel_index)





