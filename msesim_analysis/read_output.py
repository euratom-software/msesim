import idlbridge as idl

#Directory of msesim output
idl.execute("restore, '/home/sgibson/PycharmProjects/msesim/runs/benchmarkMAST/output/data/MAST_28101_t200ms_filtered.dat' , /VERBOSE")

class MSESIM_Output(object):

    def __init__(self):

        self.key_names = None
        self.object_names = None

    def read_data(self):
        self.data = {}

        for key_name, object_name in zip(self.key_names, self.object_names):
            self.data[object_name] = idl.get(key_name)

class Channel(MSESIM_Output):

    def __init__(self):

        self.key_names = ["CHANNELS", "CHANID"]
        self.object_names = ["channels", "channel_id"]
        self.read_data()

class CentralPoints(MSESIM_Output):

    def __init__(self):

        self.key_names = ["xyz0", "B_v0", "B_xyz", "B_w", "B_vec",
                          "C_xyz", "C_k", "vec0", "Bfld0", "Efld0",
                          "Dshift0", "Sshift0", "alpha0", "psi0"]

        self.object_names = ['central_coordinates', 'beam_velocity_vector', 'beam_duct_coordinates', 'half_beam_sampling_width', 'beam_axis_vector',
                            'collection_lens_coordinates', 'optical_axis', 'emission_vector', 'bfield_vector', 'efield_vector',
                             'doppler_shift', 'max_stark_shift', 'polarisation_angles', 'psi_normalised']
        self.read_data()

class GridPoints(MSESIM_Output):

    def __init__(self):

        self.key_names = ["gp_xyz", "gp_vel", "gp_vec",
                          "gp_bfld", "gp_efld", "gp_alpha",
                          "gp_psi", "gp_emis"]

        self.object_names = ['grid_coordinates', 'beam_velocity_vector', 'emission_vector',
                             'bfield_vector', 'efield_vector', 'polarisation_angles',
                             'psi_normalised', 'emission_intensity' ]

        self.read_data()

class Resolution(MSESIM_Output):

    def __init__(self):

        self.key_names = ["R", "Z", "psi", "RZpsi", "Rm",
                          "RZ_emis", "R_emis", "psi_emis",
                          "psi_res", "R_res"]

        self.object_names = ['R', 'Z', 'psi', 'psi(R,Z)', 'magnetic_axis',
                             'emission_intensity(RZ)', 'emission_intensity(R)', 'emission_intensity(psi)',
                             'resolution_vector(psi)', 'resolution_vector(R)']

        self.read_data()


class SpectralData(MSESIM_Output):

    def __init__(self):

        self.key_names = ["lambda", "pstokes", "sstokes", "stokes",
                          "cwlstokes", "sigstokes", "pibstokes", "pirstokes"]

        self.object_names = ['wavelength_vector', 'pi_stokes_vector', 'sigma_stokes_vector', 'total_stokes_vector',
                             'cwl_stokes', 'optimal_sigma_wavelength_stokes', 'optimal_blueshift_pi_wavelength', 'optimal_redshift_pi_wavelength']

        self.read_data()




