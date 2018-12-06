#External Imports
import numpy as np
from equi.read_efit_hdf5 import read_efit
''' B Field reconstruction from local equilibrium runs for MAST.  '''

class Equilibrium:

    def __init__(self, time_series, Rvector, Zvector, R, Z, psi_2d, psi_n, B_phi, B_Z, B_R, Bvac_magaxis, psi_lcfs, boundary_coordinates, j_phi_2d, jphi_1d, r_magaxis, z_magaxis):
        self.time_series = time_series
        self.Rvector = Rvector
        self.Zvector = Zvector
        self.R = R
        self.Z = Z
        self.psi_2d = psi_2d
        self.psi_n = psi_n
        self.B_phi = B_phi
        self.B_Z = B_Z
        self.B_R = B_R
        self.Bvac_magaxis = Bvac_magaxis
        self.psi_lcfs = psi_lcfs
        self.boundary_coordinates = boundary_coordinates
        self.j_phi_2d = j_phi_2d
        self.jphi_1d = jphi_1d
        self.r_magaxis = r_magaxis
        self.z_magaxis = z_magaxis

def get_slice(efit_headers, specified_time):

    time_series = np.array([efit_headers['times']])

    time_slice = abs(time_series - specified_time).argmin()

    return time_slice, time_series

def vectors(efit_vectors):

    Rvector = np.array([list(efit_vectors['rVector'])])
    Zvector = np.array([list(efit_vectors['zVector'])])

    return Rvector, Zvector

def grid(Rvector, Zvector):
    R,Z = np.meshgrid(Rvector, Zvector, indexing='ij')
    return R,Z

def magnetic_axes(global_parameters, time_slice):
    magnetic_axis = np.array([global_parameters['magneticAxis']])
    R_magaxis = magnetic_axis[0,time_slice,0]
    Z_magaxis = magnetic_axis[0,time_slice,1]

    return R_magaxis, Z_magaxis

def plasmacurrent(global_parameters):
    plasma_current = global_parameters['plasmaCurrent']
    return plasma_current

def boundary_coords(efit_geometry, time_slice):
    boundary_coordinates = np.array([efit_geometry['boundaryCoordinates']])
    return boundary_coordinates[0,time_slice,:,:]

def psi_RZ(efit_2d_profiles,time_slice):
    #psi as a function of R,Z
    psi_2d = np.array([efit_2d_profiles['poloidalFlux']])
    return psi_2d[0,time_slice,:,:]

def psi_magnetic_axis(global_parameters, time_slice):
    #value of psi at the magnetic axis
    psi_axis = np.array([global_parameters['psiAxis']])
    return psi_axis[0,time_slice]

def psi_boundary(global_parameters, time_slice):
    #value of psi at last closed flux surface
    psi_lcfs = np.array([global_parameters['psiBoundary']])
    return psi_lcfs[0,time_slice]

def toroidal_flux_function(efit_outputs, time_slice):
    #current flux function = RBpi/mu0
    f_current = np.array([efit_outputs['fluxFunctionProfiles']['fDia']])
    return f_current[0,time_slice,:]

def BVac_magnetic_axis(global_parameters, time_slice):
    #Value of the vaccuum magnetic field at the magnetic axis
    Bvac_magaxis = np.array([global_parameters['bvacRmag']])
    return Bvac_magaxis[0,time_slice]

def getjphi_2d(efit_2d_profiles, time_slice):
    jphi_2d = np.array([efit_2d_profiles['jphi']])
    return jphi_2d[0,time_slice,:,:]

def getjphi_1d(radial_profiles, time_slice):
    jphi_1d = np.array([radial_profiles['jphi']])
    return jphi_1d[0,time_slice,:]

def psi_normalised(psi_2d, psi_axis, psi_lcfs):
    #normalised poloidal flux, such that psi_n **2 = rho (where rho = 1 at the LCFS)
    psi_n = (psi_2d - psi_axis)/(psi_lcfs-psi_axis)
    return psi_n

def grad_psi_Z(psi_2d, Z):
    #derivative of 2d flux function

    dPsidZ = np.gradient(psi_2d, axis=1)/np.gradient(Z[0,:])

    return dPsidZ

def grad_psi_R(psi_2d, R):
    #derivative of 2d flux function with respect to R

    dPsidR = np.gradient(psi_2d, axis=0)/np.gradient(R[:,0])
    return dPsidR

def current_flux_function_2d(Rvector, Zvector, psi_axis, psi_lcfs, psi_2d, Bvac_magaxis, R_magaxis, f_current):

    psi_grid = np.linspace(psi_axis, psi_lcfs, len(f_current))

    f_RZ = np.zeros((len(Rvector[0,:]),len(Zvector[0,:])))

    for r in range(len(Rvector[0,:])):
        for z in range(len(Zvector[0,:])):
            if psi_2d[r,z] < psi_grid[-1] and psi_2d[r,z] > psi_grid[0]:
                f_RZ[r,z] = f_current(psi_2d[r,z])
            else:
                f_RZ[r,z] = Bvac_magaxis*1.0

    return f_RZ

def B_toroidal(f_RZ, R):
    B_phi = f_RZ/R
    return B_phi

def BZ(R, dPsidR):
    B_Z = (1/R) * dPsidR
    return B_Z

def BR(R, dPsidZ):
    B_R = (-1/R)*dPsidZ
    return B_R

def equilibrium(specified_time, efit_filepath):

    efitOut, efit_headers, efit_inputs, efit_outputs, global_parameters, radial_profiles, efit_geometry, efit_vectors, efit_2d_profiles = read_efit(efit_filepath)

    time_slice, time_series = get_slice(efit_headers, specified_time)
    Rvector, Zvector = vectors(efit_vectors)
    R,Z = grid(Rvector, Zvector)
    f_current = toroidal_flux_function(efit_outputs, time_slice)
    plasma_current = plasmacurrent(global_parameters)
    Bvac_magaxis = BVac_magnetic_axis(global_parameters, time_slice)
    R_magaxis, Z_magaxis = magnetic_axes(global_parameters, time_slice)
    psi_2d = psi_RZ(efit_2d_profiles,time_slice)
    psi_axis = psi_magnetic_axis(global_parameters,time_slice)
    psi_lcfs = psi_boundary(global_parameters,time_slice)
    psi_n = psi_normalised(psi_2d, psi_axis, psi_lcfs)

    dPsidZ = grad_psi_Z(psi_2d, Z)
    dPsidR = grad_psi_R(psi_2d, R)

    j_phi_2d = getjphi_2d(efit_2d_profiles, time_slice)
    jphi_1d = getjphi_1d(radial_profiles, time_slice)

    f_RZ = current_flux_function_2d(Rvector, Zvector, psi_axis, psi_lcfs, psi_2d, Bvac_magaxis, R_magaxis, f_current)

    B_phi = B_toroidal(f_RZ, R)
    B_Z = BZ(R, dPsidR)
    B_R = BR(R, dPsidZ)

    r_magaxis, z_magaxis = magnetic_axes(global_parameters, time_slice)
    boundary_coordinates = boundary_coords(efit_geometry, time_slice)

    efit = Equilibrium(time_series, Rvector, Zvector, R, Z, psi_2d, psi_n, B_phi, B_Z, B_R, Bvac_magaxis, psi_lcfs, boundary_coordinates, j_phi_2d, jphi_1d, r_magaxis, z_magaxis)

    return efit

