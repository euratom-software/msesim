import idlbridge as idl
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def create_channel_dataset():

    central_pos = xr.DataArray(data=idl.get('xyz0'), dims=('channel', 'xyz'), name='central_pos')
    central_pos.attrs = {'units': 'm', 'long_name': 'Central positions of the cylindrical volume element in machine co-ordinates'}

    lens_coords = xr.DataArray(data=idl.get('c_xyz'), dims=('xyz'), name='lens_pos')
    lens_coords.attrs = {'units': 'm', 'long_name': 'Collection lens position in machine co-ordinates'}

    lens_vector = xr.DataArray(data=idl.get('c_k'), dims=('xyz'), name='lens_vector')
    lens_vector.attrs = {'units': 'Arb Units', 'long_name': 'Normalised optical axis'}

    beam_pos = xr.DataArray(data=idl.get('b_xyz'), dims=('xyz'), name='beam_pos')
    beam_pos.attrs = {'units': 'm', 'long_name': 'Beam duct position in machine co-ordinates'}

    beam_vector = xr.DataArray(data=idl.get('b_vec'), dims=('xyz'), name='beam_vector')
    beam_vector.attrs = {'units': 'Arb Units', 'long_name': 'Normalised beam axis'}


    beam_halfwidth = xr.DataArray(data=idl.get('b_w'), name='beam_halfwidth')
    beam_halfwidth.attrs = {'units': 'Arb Units', 'long_name': 'Normalised beam axis'}

    beam_velocity = xr.DataArray(data=idl.get('b_v0'), name='beam_velocity')
    beam_halfwidth.attrs = {'units': 'm/s', 'long_name': 'Beam velocity in m/s calculated from beam voltage'}


    E_field = xr.DataArray(data=idl.get('Efld0'), dims=('channel', 'xyz'), name='E_field')
    E_field.attrs = {'units': 'V/m', 'long_name': 'Electric field vector at central points'}

    B_field = xr.DataArray(data=idl.get('Bfld0'), dims=('channel', 'xyz'), name='B_field')
    B_field.attrs = {'units': 'T', 'long_name': 'Magnetic field vector at central points'}

    doppler_shift = xr.DataArray(data=idl.get('dshift0'), dims=('channel'), name='doppler_shift')
    doppler_shift.attrs = {'units': 'angstroms', 'long_name': 'Doppler shift of emission for each channel due to atomic motion'}

    stark_shift = xr.DataArray(data=idl.get('sshift0'), dims=('channel'), name='stark_shift')
    stark_shift.attrs = {'units': 'angstroms', 'long_name': 'Stark shift of emission for each channel due to E field'}

    alpha_central = xr.DataArray(data=idl.get('alpha0'), dims=('channel'), name='alpha_central')
    alpha_central.attrs = {'units': 'Radians', 'long_name': 'Polarisation angle at central points'}

    psi_central = xr.DataArray(data=idl.get('psi0')[:,0], dims=('channel'), name='psi_central')
    psi_central.attrs = {'units': 'Arb Units', 'long_name': 'Normalised poloidal flux at central positions'}

    channel_dataset = xr.merge([central_pos, lens_coords, lens_vector,
                                beam_pos, beam_vector, beam_halfwidth, beam_velocity,
                                E_field, B_field,
                                doppler_shift, stark_shift,
                                alpha_central, psi_central])

    return channel_dataset.assign_coords(xyz=['x', 'y', 'z'], channel=np.arange(40))


def create_grid_dataset():

    gp_pos = xr.DataArray(data=idl.get('gp_xyz'), dims=('channel','grid_point','xyz'), name='grid_pos')
    gp_pos.attrs = {'units': 'm', 'long_name': 'Grid positions of the cylindrical volume element in machine co-ordinates'}

    gp_vel = xr.DataArray(data=idl.get('gp_vel'), dims=('channel', 'grid_point','xyz'), name='gp_vel')
    gp_vel.attrs = {'units': 'm', 'long_name': 'Beam velocity vector at grid points in machine co-ordinates'}

    gp_vec = xr.DataArray(data=idl.get('gp_vec'), dims=('channel', 'grid_point','xyz'), name='gp_vec')
    gp_vec.attrs = {'units': 'm', 'long_name': 'Beam emission vector at grid points in machine co-ordinates'}

    gp_bfld = xr.DataArray(data=idl.get('gp_bfld'), dims=('channel', 'grid_point','xyz'), name='gp_bfld')
    gp_bfld.attrs = {'units': 'V/m', 'long_name': 'B field vector at grid points in machine co-ordinates'}

    gp_efld = xr.DataArray(data=idl.get('gp_efld'), dims=('channel', 'grid_point','xyz'), name='gp_efld')
    gp_efld.attrs = {'units': 'T', 'long_name': 'E field vector at grid points in machine co-ordinates'}

    gp_alpha = xr.DataArray(data=idl.get('gp_alpha'), dims=('channel', 'grid_point'), name='g_alpha')
    gp_alpha.attrs = {'units': 'Radians', 'long_name': 'Polarisation angle at grid points'}

    gp_psi = xr.DataArray(data=idl.get('gp_psi'), dims=('channel', 'grid_point'), name='gp_psi')
    gp_psi.attrs = {'units': 'Arb Units', 'long_name': 'Normalised poloidal flux at grid points'}

    gp_emis = xr.DataArray(data=idl.get('gp_emis'), dims=('channel','grid_point'), name='gp_emis')
    gp_emis.attrs = {'units': 'photons/m$^{3}$/sr', 'long_name': 'Emission intensity (photons/volume element/solid angle) at grid points'}

    gp_dataset = xr.merge([gp_pos, gp_vel, gp_vec,
                           gp_bfld, gp_efld, gp_alpha,
                           gp_psi, gp_emis])

    return gp_dataset.assign_coords(xyz=['x', 'y', 'z'], channel=np.arange(40))

def create_stokes_dataset():

    wavelength = xr.DataArray(data=idl.get('lambda'), dims=('wavelength'), name='lambda')
    wavelength.attrs = {'units': 'Angstrom', 'long_name': 'Wavelength'}

    pi_stokes = xr.DataArray(data=idl.get('pstokes'), dims=('channel', 'stokes', 'wavelength'), name='pi_stokes')

    sigma_stokes = xr.DataArray(data=idl.get('sstokes'), dims=('channel', 'stokes', 'wavelength'), name='sigma_stokes')

    total_stokes = xr.DataArray(data=idl.get('stokes'), dims=('channel', 'stokes', 'wavelength'), name='total_stokes')

    stokes_dataset = xr.merge([wavelength, pi_stokes, sigma_stokes, total_stokes])
    stokes_dataset = stokes_dataset.assign_coords(wavelength=wavelength)

    #Calculate the polarisation angle from the total stokes vector

    polarisation_angle = 0.5*np.arctan(np.sum(total_stokes.isel(stokes=2),axis=1)/np.sum(total_stokes.isel(stokes=1),axis=1))

    plt.figure()
    plt.plot(polarisation_angle*180./np.pi)
    plt.show()

    sigma_polarised = (sigma_stokes.isel(stokes=1)**2 + sigma_stokes.isel(stokes=2)**2)**0.5

    pi_polarised = (pi_stokes.isel(stokes=1)**2 + pi_stokes.isel(stokes=2)**2)**0.5

    gamma_dataset = xr.merge([wavelength, polarisation_angle, sigma_polarised, pi_polarised])

    gamma_dataset = gamma_dataset.assign_coords(channel=np.arange(40), wavelength=wavelength)

    return stokes_dataset, gamma_dataset

def create_cwl_stokes_dataset():

    cwl_stokes = xr.DataArray(data=idl.get('cwlstokes'), dims=('channel', 'cwlstokes'), name='cwl_stokes')
    cwl_stokes.attrs =  {'long_name': 'Stokes Component'}

    sigma_optimal_stokes = xr.DataArray(data=idl.get('sigstokes'), dims=('channel', 'cwlstokes'), name='sigma_optimal_stokes')
    sigma_optimal_stokes.attrs =  {'long_name': 'Stokes Component'}

    pired_optimal_stokes = xr.DataArray(data=idl.get('pirstokes'), dims=('channel', 'cwlstokes'), name='pired_optimal_stokes')
    pired_optimal_stokes.attrs =  {'long_name': 'Stokes Component'}

    piblue_optimal_stokes = xr.DataArray(data=idl.get('pibstokes'), dims=('channel', 'cwlstokes'), name='piblue_optimal_stokes')
    piblue_optimal_stokes.attrs =  {'long_name': 'Stokes Component'}


    cwl_stokes_dataset = xr.merge([cwl_stokes, sigma_optimal_stokes,
                                   pired_optimal_stokes, piblue_optimal_stokes])

    return cwl_stokes_dataset

def create_resolution_dataset():


    R = xr.DataArray(data=idl.get('R'), dims=('R'), name='r')
    R.attrs = {'units': 'm', 'long_name': 'Major Radius'}

    Z = xr.DataArray(data=idl.get('Z'), dims=('Z'), name='z')
    Z.attrs = {'units': 'm', 'long_name': 'Vertical extent'}

    psi_rz = xr.DataArray(data=idl.get('rzpsi'), dims=('Z', 'R'), name='psi_rz')
    psi_rz.attrs = {'units': 'wb/m (?)', 'long_name': 'Poloidal flux $\Psi$(R,Z)'}

    emission_rz = xr.DataArray(data=idl.get('rz_emis'), dims=('channel', 'Z', 'R'), name='emission_rz')
    emission_rz.attrs = {'units': 'photons/m$^{3}$/sr', 'long_name': 'Emission intensity (R,Z)'}

    psi_n = xr.DataArray(data=idl.get('psi'), dims=('psin'), name='psi_n')
    psi_n.attrs = {'units': 'Arb Units', 'long_name': 'Normalised poloidal flux'}

    emission_psi = xr.DataArray(data=idl.get('psi_emis'), dims=('channel', 'psin'), name='emission_psi')
    emission_psi.attrs = {'units': 'photons/m$^{3}$/sr', 'long_name': 'Emission intensity ($\Psi$)'}

    r_resolution = xr.DataArray(data=idl.get('r_res'), dims=('channel', 'nr'), name='r_resolution')
    r_resolution.attrs = {'units': 'm', 'long_name': 'Spatial resolution'}

    emission_dataset = xr.merge([R, Z, psi_rz, emission_rz])

    resolution_dataset = xr.merge([psi_n, emission_psi, r_resolution])

    return emission_dataset, resolution_dataset

def create_datasets(fname):

    idl.execute("restore, '{0}' , /VERBOSE".format(fname))
    channel_ds = create_channel_dataset()
    grid_ds = create_grid_dataset()
    stokes_dataset, gamma_dataset = create_stokes_dataset()
    emission_ds, res_ds = create_resolution_dataset()

    return channel_ds, grid_ds, stokes_dataset, gamma_dataset, res_ds, emission_ds


def spectrum(fname, channel, title):
    channel_ds, grid_ds, stokes_ds, gamma_ds, res_ds, emission_ds = create_datasets(fname)
    try:
        stokes_ds['total_stokes'].isel(stokes=0, channel=channel).plot(label='R={}m, {}'.format(np.round(res_ds['r_resolution'].isel(channel=channel, nr=2).values,3), title))
        # stokes_ds['total_stokes'].isel(stokes=1, channel=channel).plot(label='S1')
        # stokes_ds['total_stokes'].isel(stokes=2, channel=channel).plot(label='S2')
        # stokes_ds['total_stokes'].isel(stokes=3, channel=channel).plot(label='S3')
        plt.ylabel('Intensity [ph/s]', fontsize=18)
        plt.xlabel('Wavelength [Angstrom]', fontsize=18)
        plt.legend(loc=2, prop={'size': 12})
        plt.title('{}'.format(title))
        #plt.xlim(6598,6602.5)
        #plt.ylim(-4500,8500)
        # plt.show()
    except IndexError:
        print('Channel number must be between 0 and 40!')

def create_gamma_dataset():

    wavelength = xr.DataArray(data=idl.get('lambda'), dims=('wavelength'), name='lambda')
    wavelength.attrs = {'units': 'Angstrom', 'long_name': 'Wavelength'}

    r_resolution = xr.DataArray(data=idl.get('r_res'), dims=('channel', 'nr'), name='r_resolution')
    r_resolution.attrs = {'units': 'm', 'long_name': 'Spatial resolution'}

    r_coords = r_resolution.isel(nr=2)

    r_coords = xr.DataArray(data=r_coords, name='r_coords')

    #Calculate the polarisation angle from the total stokes vector
    total_stokes = xr.DataArray(data=idl.get('stokes'), dims=('channel', 'stokes', 'wavelength'), name='total_stokes')

    polarisation_angle = -1*0.5*np.arctan(np.sum(total_stokes.isel(stokes=2),axis=1)/np.sum(total_stokes.isel(stokes=1),axis=1))*180./np.pi

    polarisation_angle = xr.DataArray(data=polarisation_angle, dims=('channel'), name='gamma')

    gamma_dataset = xr.merge([wavelength, polarisation_angle, r_coords])

    gamma_dataset = gamma_dataset.assign_coords(channel=np.arange(40))

    return gamma_dataset


def create_dataset(fname):

    idl.execute("restore, '{0}' , /VERBOSE".format(fname))
    gamma_dataset = create_gamma_dataset()

    return gamma_dataset

mastu_1ma = '/work/sgibson/msesim/runs/conventional_mse_mastu_fiesta1MA/output/data/conventional_mse_mastu.dat'
jet = '/work/sgibson/msesim/runs/JET/KS5_test/output/data/KS5_87123.dat'

mast = '/work/sgibson/msesim/runs/conventional_mse_mast/output/data/MAST_27193_0.2_lpmse.dat'
mast_deltabeam = '/work/sgibson/msesim/runs/mse_mast_test/output/data/MAST_27193_0.2_lpmse.dat'

mastu_k25 = '/work/sgibson/msesim/runs/conventional_mse_mastu_K25_scenario/output/data/conventional_mse_mastu_k25_scenario.dat'
mastu_mastlike = '/work/sgibson/msesim/runs/conventional_mse_mastu_mastlike/output/data/conventional_mse_mastu.dat'

test = '/work/sgibson/msesim/runs/mse_mast_test/output/data/MAST_27193_0.2_lpmse.dat'

ds_1ma = create_dataset(mastu_1ma)
plt.plot(ds_1ma['r_coords'], ds_1ma['gamma'], color='C0', label='Ip = 1MA, MAST-U plasma, 3-4 beam')

ds_k25 = create_dataset(mastu_k25)
plt.plot(ds_k25['r_coords'], ds_k25['gamma'], color='C1', label='Ip = 1MA, MAST-U plasma, 2 beam')

ds_mastlike = create_dataset(mastu_mastlike)
plt.plot(ds_mastlike['r_coords'], ds_mastlike['gamma'], color='C2', label='Mast-like plasma')

plt.ylabel('Polarisation Angle (Degrees)', fontsize=28)
plt.xlabel('R (m)', fontsize=28)
plt.legend(loc=2, prop={'size': 20})
plt.show()





# spectrum(mastu_1ma, channel=10, title='MAST-U 1MA')

#spectrum(test, channel=20, title='test tiny beam width')
# spectrum(mastu_mastlike, channel=20, title='MAST-like plasma, 65kV beam')
# spectrum(mastu_k25, channel=20, title='MAST-U plasma, 70kV beam')
# spectrum(mastu_1ma, channel=20, title='MAST-U 1MA plasma, 75kV beam')
# # plt.show()
#
# channel_ds, grid_ds, stokes_dataset, gamma_dataset, res_ds, emission_ds = create_datasets(mastu_1ma)


