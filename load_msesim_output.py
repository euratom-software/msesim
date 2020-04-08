import idlbridge as idl
import numpy as np

from scipy.io import readsav

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import pyuda

client = pyuda.Client()

plt.ion()

class MSESIM(object):

    def __init__(self):

        self.object_names = (
                'intersection_coordinates', 'beam_velocity', 'beam_duct_coordinates', 'half_beam_sampling_width',
                'beam_axis_vector',
                'collection_lens_coordinates', 'optical_axis', 'emission_vector', 'bfield_vector', 'efield_vector',
                'doppler_shift', 'max_stark_shift', 'polarisation_angles', 'psi_at_intersection_coordinates',

                'grid_coordinates', 'grid_beam_velocity_vectors', 'grid_emission_vectors',
                'grid_bfield_vector', 'grid_efield_vector', 'grid_polarisation_angles',

                'grid_psi_normalised', 'grid_beam_emission_intensity',

                'R', 'Z', 'psi', 'psi(R,Z)', 'magnetic_axis',

                'emission_intensity(RZ)', 'emission_intensity(R)', 'emission_intensity(psi)',
                'radial_resolution(psi)', 'radial_resolution(R)',

                'wavelength_vector', 'pi_stokes', 'sigma_stokes', 'total_stokes',
                'cwl_stokes', 'optimal_sigma_wavelength_stokes', 'optimal_blueshift_pi_wavelength_stokes',
                'optimal_redshift_pi_wavelength_stokes')

        self.data = self.load_msesim_spectrum()
        self.parameters_1d()


    @staticmethod
    def _find_nearest(array, value):

        if value < array.min() or value > array.max():
            raise IndexError("Requested value is outside the range of the data. The range of the data is from {}m to {}m".format(array.min(),array.max()))

        index = np.searchsorted(array, value, side="left")

        print(index)

        if (value - array[index])**2 < (value - array[index + 1])**2:
            return index
        else:
            return index + 1

    def load_msesim_spectrum(self):

        """
        Stores the output of an msesim run (stored in a .dat file) in a dictionary for use in python.
        :return: Dictionary of outputs from msesim. Retrievable using the object_names as given below.
        """

        self.data = {}

        key_names = ("xyz0", "B_v0", "B_xyz", "B_w", "B_vec",
                     "C_xyz", "C_k", "vec0", "Bfld0", "Efld0",
                     "Dshift0", "Sshift0", "alpha0", "psi0",

                     "gp_xyz", "gp_vel", "gp_vec", "gp_bfld",
                     "gp_efld", "gp_alpha", "gp_psi", "gp_emis",

                     "R", "Z", "psi", "RZpsi", "Rm",
                     "RZ_emis", "R_emis", "psi_emis", "psi_res",
                     "R_res",

                     "lambda", "pstokes", "sstokes", "stokes",
                     "cwlstokes", "sigstokes", "pibstokes", "pirstokes")

        for key_name, object_name in zip(key_names, self.object_names):
            self.data[object_name] = idl.get(key_name)

        return self.data

    def parameters_1d(self):

        """
        Calculate useful parameters using the Stokes vectors, including:

        - Polarised Fraction
        - Total Unpolarised Fraction
        - Linearly polarised fraction
        - Circularly polarised fraction
        - Polarisation angle at the central wavelengths

        :return:
        """

        self.stokes_total = np.array(self.data["cwl_stokes"])

        self.S0 = self.stokes_total[:,0]
        self.S1 = self.stokes_total[:,1]
        self.S2 = self.stokes_total[:,2]
        self.S3 = self.stokes_total[:,3]
        self.cwl = self.stokes_total[:,4]

        self.wavelength = np.array(self.data["wavelength_vector"])/10 #in nanometers

        self.major_radius = np.array(self.data["resolution_vector(R)"])[:,2]


        self.polarised_fraction = np.sqrt(self.S1**2 + self.S2**2+self.S3**2)/self.S0
        self.total_unpolarised = (self.S0 - np.sqrt(self.S1**2 + self.S2**2 + self.S3**2))/self.S0
        self.total_circular = np.sqrt(self.S3**2 - (self.total_unpolarised**2*self.S0))/self.S0

        self.LPF = np.sqrt(self.S1**2 + self.S2**2)/self.S0
        self.CPF = self.total_circular/self.S0

        self.circular_total = np.sqrt(self.S3**2)/self.S0 #total circular polarisation fraction

        self.gamma = -0.5*np.arctan2(self.S2,self.S1)*(180./np.pi) #polarisation angle

    def plot_spectrum(self, radius):

        # Find the radii for the spectrum you want to plot
        idx = self._find_nearest(self.major_radius, value=radius)

        #manager = plt.get_current_fig_manager()
        #manager.window.showMaximized()

        fig= plt.figure(1)
        gs1 = gridspec.GridSpec(nrows=3, ncols=3)
        ax1 = fig.add_subplot(gs1[:-1, :])
        plt.plot(self.wavelength, self.S0[idx,:].T, color='black', label='$I_{\mathrm{total}}$')
        plt.plot(self.wavelength, np.sqrt(self.S2[idx,:].T**2 + self.S1[idx,:].T**2), label='$I_{\mathrm{linear}}$')
        plt.plot(self.wavelength, self.S3[idx,:].T, label='$I_{\mathrm{circular}}$')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False, useMathText=True)
        plt.xlim(659.6,660.2)
        plt.legend(prop={'size': 40})
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Intensity $I$ (photons/s)', labelpad=5)

        plt.show()