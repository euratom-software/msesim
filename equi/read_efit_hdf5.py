import h5py
import numpy as np

def data_from_efitfile(efit_filepath):
    efitOut = h5py.File(efit_filepath, 'r')
    keys = list(efitOut.keys())  # hdf5 returns view like objects instead of a list, so convert it to a list
    return efitOut

def get_efit_headers(efitOut):

    efitOut_headers = ['codeVersion', 'equilibriumStatus', 'pulseNumber', 'times']
    efit_headers = {}

    for k in range(len(efitOut_headers)):
        try:
            efit_headers[efitOut_headers[k]] = efitOut['equilibrium']['header'][efitOut_headers[k]]
        except KeyError:
            print('{} does not exist for this run'.format(efitOut_headers[k]))

    return efit_headers

def get_efit_inputs(efitOut):

    efitOut_inputs = ['BVacRadiusProduct', 'codeControls', 'constraints', 'currentSetup', 'fitdz', 'limiter',
                      'numericalControls', 'regularGrid']
    efit_inputs = {}

    for i in range(len(efitOut_inputs)):
        try:
            efit_inputs[efitOut_inputs[i]] = efitOut['equilibrium']['input'][efitOut_inputs[i]]
        except KeyError:
            print('{} does not exist for this run'.format(efitOut_inputs[i]))

    return efit_inputs

def get_efit_outputs(efitOut):

    efit_outputs = {}

    efitOut_outputs = ['boozerProfiles', 'fitDetails', 'fluxFunctionProfiles', 'geometry', 'globalParameters',
                       'profiles2D', 'radialProfiles']

    for j in range(len(efitOut_outputs)):
        try:
            efit_outputs[efitOut_outputs[j]] = efitOut['equilibrium']['output'][efitOut_outputs[j]]
        except KeyError:
            print('{} is not produced in this run'.format(efitOut_outputs[j]))

    return efit_outputs

def get_efit_global_parameters(efit_outputs):

    global_parameters = {}

    for l in range(len(list(efit_outputs['globalParameters']))):
        key = list(efit_outputs['globalParameters'])[l]
        global_parameters[key] = list(efit_outputs['globalParameters'][key])

    return global_parameters

def get_efit_radial_profiles(efit_outputs):

    radial_profiles = {}

    for l in range(len(list(efit_outputs['radialProfiles']))):
        key = list(efit_outputs['radialProfiles'])[l]
        radial_profiles[key] = list(efit_outputs['radialProfiles'][key])

    return radial_profiles

def get_efit_geometry(efit_outputs):

    efit_geometry = {}

    for n in range(len(list(efit_outputs['geometry']))):
        key = list(efit_outputs['geometry'])[n]
        efit_geometry[key] = list(efit_outputs['geometry'][key])

    return efit_geometry

def get_grid_vectors(efit_inputs):

    efit_vectors = {}

    for l in range(len(list(efit_inputs['regularGrid']))):
        key = list(efit_inputs['regularGrid'])[l]
        efit_vectors[key] = list(efit_inputs['regularGrid'][key])

    return efit_vectors

def profiles_2D(efit_outputs):

    efit_2d_profiles = {}

    for m in range(len(list(efit_outputs['profiles2D']))):
        key = list(efit_outputs['profiles2D'])[m]
        efit_2d_profiles[key] = list(efit_outputs['profiles2D'][key])

    return efit_2d_profiles

def read_efit(efit_filepath):

    efitOut = data_from_efitfile(efit_filepath)
    efit_headers = get_efit_headers(efitOut)
    efit_inputs = get_efit_inputs(efitOut)
    efit_outputs = get_efit_outputs(efitOut)
    global_parameters = get_efit_global_parameters(efit_outputs)
    radial_profiles = get_efit_radial_profiles(efit_outputs)
    efit_geometry = get_efit_geometry(efit_outputs)
    efit_vectors = get_grid_vectors(efit_inputs)
    efit_2d_profiles = profiles_2D(efit_outputs)

    return efitOut, efit_headers, efit_inputs, efit_outputs, global_parameters, radial_profiles, efit_geometry, efit_vectors, efit_2d_profiles
