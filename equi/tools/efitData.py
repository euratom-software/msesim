#code from Lucy Kogan/Ivan Lupelli

import os
import time
import importlib

import numpy as np
import xarray as xr
import netCDF4 as nc

from equi.tools.efitUtils import chi2calc, errorsCalc, writeColumns

class EfitData():
    def __init__(self, include_converged_times=True, include_vacuum_times=True,
                 include_non_converged_times=True, include_fatal_times=False,
                 subsample=1, sample_count=None, time_min=None, time_max=None,
                 filename="efitOut", netcdf=False, hdf5=True, shot=None, tokamak='', debug=0):

        self.nanosecond = 1e-9
        self.debug = debug

        if netcdf:

        # Dictionary of xarray objects mirroring hierarchy of the NetCDF file.
        # Initially setup empty dictionary structure
            self.data = {"top": None,
                         "constraints": {"plasmaCurrent": None,
                                         "pfCircuits" : None,
                                         "fluxLoops" : None,
                                         "magneticProbes" : None,
                                         "mse" : None,
                                         "lineIntegratedDensities" : None,
                                         "discreteDensities" : None,
                                         "pressures" : None,
                                         "faradayRotationChannels" : None,
                                         "diamagneticFlux" : None,
                                         "boundaries" : None},
                         "ironModel": None,
                         "regularGrid": None,
                         "geometry": {"limiter": None, "pfSystem": None},
                         "profiles2D": None,
                         "globalParameters": None,
                         "radialProfiles": None,
                         "fluxFunctionProfiles": None,
                         "separatrixGeometry": None,
                         "numerics": {"numericalControls": None,
                                      "pp": None,
                                      "ff": None,
                                      "ww": None,
                                      "ne": None,
                                      "numericalDetails": None}}

            self.time_resolution = 0

        if hdf5:

            # Dictionary of xarray objects mirroring hierarchy of the old hdf5 file.
            # Initially setup empty dictionary structure
            self.data = { "equilibrium": {
                            "header": {"codeVersion":None,
                                     "equilibriumStatus": None,
                                     "pulseNumber": None,
                                     "times":None
                                     },
                            "input": {'BVacRadiusProduct': None,
                                    'codeControls': None,
                                    'constraints': None,
                                    'currentSetup': None,
                                    'fitdz': None,
                                    'limiter': None,
                                    'numericalControls': None,
                                    'regularGrid': None
                                    },
                            "output": {'boozerProfiles': None,
                                     'fitDetails': None,
                                     'fluxFunctionProfiles': None,
                                     'geometry': None,
                                     'globalParameters': None,
                                     'profiles2D': None,
                                     'radialProfiles': None
                                      }}}

            self.time_resolution = 0

        if shot is None:
            self._read_data(filename, include_converged_times=include_converged_times, include_vacuum_times=include_vacuum_times,
                            include_non_converged_times=include_non_converged_times, include_fatal_times=include_fatal_times,
                            subsample=subsample, sample_count=sample_count, time_min=time_min, time_max=time_max)
            self.calc_chi2()
            self.calc_errors()
        elif shot is not None and tokamak.lower() != '':
            tok_module = importlib.import_module('efitplot_{}'.format(tokamak.lower()))
            all_data = tok_module.read_tok_data(shot, include_converged_times=include_converged_times, include_vacuum_times=include_vacuum_times,
                                                 include_non_converged_times=include_non_converged_times, include_fatal_times=include_fatal_times,
                                                 subsample=subsample, sample_count=sample_count, time_min=time_min, time_max=time_max)
            self.data = all_data

    def _read_group(self, filename, group, itimes=None):
        try:
            xrdata = xr.open_dataset(filename, group=group)

            if itimes is not None:
                xrdata = xrdata.isel(time=itimes)

            # This is needed since when reading in at group level the coordinate values are not picked up from higher up levels
            if "time" in xrdata.dims:
                xrdata = xrdata.assign(time = self.data["top"].time.values)


        except OSError:
            if self.debug > 0:
                print("<<WARNING>> Could not find group {} in netcdf file {}".format(group, filename))

            xrdata = None

        return xrdata

    def _read_data(self, filename, include_converged_times=True, include_vacuum_times=True,
                   include_non_converged_times=True, include_fatal_times=False,
                   subsample=1, sample_count=None, time_min=None, time_max=None):

        # --------------------------
        # Read in top level status and time
        # And select time slices that were requested
        # --------------------------
        top = xr.open_dataset(filename)

        # Change time to floats, in seconds
        top.time.values = top.time.values.astype(np.float64) * self.nanosecond

        equilibrium_status = top.equilibriumStatusInteger.values

        flag = np.zeros(equilibrium_status.size, dtype="bool")
        if include_converged_times:
            flag = np.logical_or(flag, equilibrium_status == 1)
        if include_vacuum_times:
            flag = np.logical_or(flag, equilibrium_status == 2)
        if include_non_converged_times:
            flag = np.logical_or(flag, equilibrium_status == -1)
        if include_fatal_times:
            flag = np.logical_or(flag, equilibrium_status == -2)
        if time_min is not None:
            flag = np.logical_and(flag, top.time.values.astype(np.float64) >= time_min)
        if time_max is not None:
            flag = np.logical_and(flag, top.time.values.astype(np.float64) <= time_max)

        efit_time_indices = np.arange(top.time.size)[flag]

        if efit_time_indices.size == 0:
            if self.debug > 0:
                print("<<ERROR>> No timeslices were found in the file")
                print("          You requested to include converged? {} vacuum? {} non-converged? {} fatal? {}".format(include_converged_times,
                                                                                                                       include_vacuum_times,
                                                                                                                       include_non_converged_times,
                                                                                                                       include_fatal_times))
                if time_min is not None and time_max is not None:
                    print("          For times between {} and {} seconds.".format(time_min, time_max))
                elif time_min is not None:
                    print("          For times >= {} seconds.".format(time_min))
                elif time_max is not None:
                    print("          For times <= {} seconds.".format(time_max))
            return

        # Subsample the time-slices
        efit_time_indices = efit_time_indices[::subsample]

        # specify the number of time-slices (optional)
        if sample_count is not None:
            efit_time_indices = efit_time_indices[:sample_count]

        top = top.isel(time=efit_time_indices)

        if top.time.size > 1:
            self.time_resolution = (top.time.values[1] - top.time.values[0])

        mintime = top.time.min().values
        maxtime = top.time.max().values

        if len(efit_time_indices) == 1:
            print("<<INFO>> EFIT++ data available for a single time slice: {} seconds".format(mintime))
        else:
            efitTimeResolution = np.diff(top.time.values)
            averageEfitTimeResolution = np.mean(efitTimeResolution)
            print("<<INFO>> EFIT++ data available from {} to {} seconds (av. tiime resol. {} ms)".format(mintime,
                                                                                                         maxtime,
                                                                                                         averageEfitTimeResolution * 1000))

        print("<<INFO>> efitplot will process: {}  time-slices in the time range: {}<t<{}".format(len(efit_time_indices), mintime, maxtime) )

        self.data["top"] = top

        # --------------------------
        # Constraints
        # --------------------------
        self.data["constraints"] = {}
        self.data["constraints"]["plasmaCurrent"] = self._read_group(filename, "/input/constraints/plasmaCurrent", itimes = efit_time_indices)
        self.data["constraints"]["pfCircuits"] = self._read_group(filename, "/input/constraints/pfCircuits", itimes = efit_time_indices)
        self.data["constraints"]["fluxLoops"] = self._read_group(filename, "/input/constraints/fluxLoops", itimes = efit_time_indices)
        self.data["constraints"]["magneticProbes"] = self._read_group(filename, "/input/constraints/magneticProbes", itimes = efit_time_indices)
        self.data["constraints"]["mse"] = self._read_group(filename, "/input/constraints/mse", itimes = efit_time_indices)
        self.data["constraints"]["lineIntegratedDensities"] = self._read_group(filename, "/input/constraints/lineIntegratedDensities", itimes = efit_time_indices)
        self.data["constraints"]["discreteDensities"] = self._read_group(filename, "/input/constraints/discreteDensities", itimes = efit_time_indices)
        self.data["constraints"]["pressures"] = self._read_group(filename, "/input/constraints/pressures", itimes = efit_time_indices)
        self.data["constraints"]["faradayRotationChannels"] = self._read_group(filename, "/input/constraints/faradayRotationChannels", itimes = efit_time_indices)
        self.data["constraints"]["diamagneticFlux"] = self._read_group(filename, "/input/constraints/diamagneticFlux", itimes = efit_time_indices)
        self.data["constraints"]["boundaries"] = self._read_group(filename, "/input/constraints/boundaries", itimes = efit_time_indices)

        # --------------------------
        # Iron Model
        # --------------------------
        self.data["ironModel"] = self._read_group(filename, "/input/constraints/ironModel/geometry")

        # --------------------------
        # Grid
        # --------------------------
        self.data["regularGrid"] = self._read_group(filename, "/input/regularGrid")

        # --------------------------
        # Geometry
        # --------------------------
        self.data["geometry"] = {}
        self.data["geometry"]["limiter"] = self._read_group(filename, "/input/limiter")
        self.data["geometry"]["pfSystem"] = self._read_group(filename, "/input/pfSystem")

        # --------------------------
        # Variables
        # --------------------------
        self.data["profiles2D"] =  self._read_group(filename, "/output/profiles2D", itimes = efit_time_indices)
        self.data["globalParameters"] = self._read_group(filename, "/output/globalParameters", itimes = efit_time_indices)
        self.data["radialProfiles"] = self._read_group(filename, "/output/radialProfiles", itimes = efit_time_indices)
        self.data["fluxFunctionProfiles"] = self._read_group(filename, "/output/fluxFunctionProfiles", itimes=efit_time_indices)

        # Can't read this in directly using xarray, since pandas can't cope with nan's for arrays of objects (ie. compound types...)
        self.data["separatrixGeometry"] = self._read_separatrix(filename, efit_time_indices)
        #self.data["separatrixGeometry"] = self._read_group(filename, "/output/separatrixGeometry", itimes=efit_time_indices)

        # --------------------------
        # Numerics
        # --------------------------
        self.data["numerics"] = {}
        self.data["numerics"]["numericalControls"] = self._read_group(filename, "/input/numericalControls", itimes=efit_time_indices)
        self.data["numerics"]["pp"] = self._read_group(filename, "/input/numericalControls/pp", itimes=efit_time_indices)
        self.data["numerics"]["ff"] = self._read_group(filename, "/input/numericalControls/ff", itimes=efit_time_indices)
        self.data["numerics"]["ww"] = self._read_group(filename, "/input/numericalControls/ww", itimes=efit_time_indices)
        self.data["numerics"]["ne"] = self._read_group(filename, "/input/numericalControls/ne", itimes=efit_time_indices)
        self.data["numerics"]["numericalDetails"] = self._read_group(filename, "/output/numericalDetails", itimes=efit_time_indices)

        # Computed boundary at midplane
        if self.data["constraints"]["boundaries"] is not None:
            bound_computed = self.data["separatrixGeometry"].rmidplaneOut
            self.data["constraints"]["boundaries"] = self.data["constraints"]["boundaries"].assign(computedSep = (['time', 'boundariesDim'], bound_computed.values.reshape((self.data["top"].time.size, 1))))

    def _read_separatrix(self, filename, itimes):
        if self.data["top"] is None:
            if self.debug > 0:
                print("<<ERROR>> Read main data before separatrix data.")
            return None

        rootgrp = nc.Dataset(filename)

        try:
            group = rootgrp.groups['output'].groups['separatrixGeometry']
            rbound = group.variables['boundaryCoords'][itimes]['R']
            zbound = group.variables['boundaryCoords'][itimes]['Z']
            xpointCount = group.variables['xpointCount'][itimes]
            xpointR = group.variables['xpointCoords'][itimes]['R']
            xpointZ = group.variables['xpointCoords'][itimes]['Z']
            strikepointR = group.variables['strikepointCoords'][itimes]['R']
            strikepointZ = group.variables['strikepointCoords'][itimes]['Z']
            limiterCoordsR = group.variables['limiterCoords'][itimes]['R']
            limiterCoordsZ = group.variables['limiterCoords'][itimes]['Z']
            rGeom = group.variables['geometricAxis'][itimes]['R']
            zGeom = group.variables['geometricAxis'][itimes]['Z']

            minorRadius = group.variables['minorRadius'][itimes]
            elongation = group.variables['elongation'][itimes]
            upperTriangularity = group.variables['upperTriangularity'][itimes]
            lowerTriangularity = group.variables['lowerTriangularity'][itimes]
            rmidplaneIn = group.variables['rmidplaneIn'][itimes]
            rmidplaneOut = group.variables['rmidplaneOut'][itimes]
        except KeyError as err:
            if self.debug > 0:
                print("<<WARNING>> Could not find output group or variable in output group: {} in netcdf file {}".format(err, filename))

            return None


        separatrix_ds = xr.Dataset({'rBoundary': (['time', 'boundaryCoordsDim'], rbound),
                                    'zBoundary': (['time', 'boundaryCoordsDim'], zbound),
                                    'minorRadius': (['time'], minorRadius),
                                    'elongation': (['time'], elongation),
                                    'upperTriangularity': (['time'], upperTriangularity),
                                    'lowerTriangularity': (['time'], lowerTriangularity),
                                    'rmidplaneIn': (['time'], rmidplaneIn),
                                    'rmidplaneOut': (['time'], rmidplaneOut),
                                    'rGeom': (['time'], rGeom),
                                    'zGeom': (['time'], zGeom),
                                    'xpointCount': (['time'], xpointCount),
                                    'xpointR': (['time', 'xpointDim'], xpointR),
                                    'xpointZ': (['time', 'xpointDim'], xpointZ),
                                    'strikepointR': (['time', 'strikepointDim'], strikepointR),
                                    'strikepointZ': (['time', 'strikepointDim'], strikepointZ),
                                    'limiterR': (['time'], limiterCoordsR),
                                    'limiterZ': (['time'], limiterCoordsZ)},
                                    coords={'time': self.data["top"].time.values})

        return separatrix_ds

    def calc_chi2(self):
        for constraint in self.data["constraints"].keys():
            if self.data["constraints"][constraint] is None:
                continue

            if constraint != "boundaries":
                self.data["constraints"][constraint] = chi2calc(self.data["constraints"][constraint])
            else:
                self.data["constraints"][constraint] = chi2calc(self.data["constraints"][constraint], computedname='computedSep', targetname='rCoords', sigmaname='rSigmas')

    def calc_errors(self):
        for constraint in self.data["constraints"].keys():
            if self.data["constraints"][constraint] is None:
                continue

            if constraint != "boundaries":
                self.data["constraints"][constraint] = errorsCalc(self.data["constraints"][constraint])
            else:
                self.data["constraints"][constraint] = errorsCalc(self.data["constraints"][constraint], computedname='computedSep', targetname='rCoords')

    def calc_bfield(self):
        if self.data["profiles2D"] is None:
            return

        # Assuming grid doesn't change with time
        Refit = self.data["profiles2D"].r.values[0, :]
        Zefit = self.data["profiles2D"].z.values[0, :]

        # Calculate psiNorm
        psi = np.transpose(self.data["profiles2D"].poloidalFlux.values,
                           axes=[0, 2, 1])

        psiNorm = ((self.data["profiles2D"].poloidalFlux - self.data["globalParameters"].psiAxis)
               / (self.data["globalParameters"].psiBoundary - self.data["globalParameters"].psiAxis))
        if (psiNorm.values.any() < 0):
            if self.debug > 0:
                print("<<WARNING>> Negative values for psinorm!!!")
            np.place(psiNorm.values, psiNorm.values < 0, 0.0)

        psiNorm.values = (psiNorm - psiNorm.min(axis=(2,1))).values

        # Calculate gradients
        dR = Refit[1] - Refit[0]
        dZ = Zefit[1] - Zefit[0]
        try:
            psiDz, psiDr = np.gradient(psi, dZ, dR, edge_order=2, axis=(1,2))
        except ValueError:
            psiDz, psiDr = np.gradient(psi, dZ, dR, axis=(1, 2))

        # Bz, Br 2D
        Bz2d = 1.0 * psiDr / Refit
        Br2d = -1.0 * psiDz / Refit

        self.data["profiles2D"] = self.data["profiles2D"].assign(psiNorm = (['time', 'rGrid', 'zGrid'], psiNorm))
        self.data["profiles2D"] = self.data["profiles2D"].assign(Bz = (['time', 'rGrid', 'zGrid'], Bz2d))
        self.data["profiles2D"] = self.data["profiles2D"].assign(Br = (['time', 'rGrid', 'zGrid'], Br2d))

        # B poloidal
        Bpol2d = np.sign(self.data["globalParameters"].plasmaCurrent) \
                 * (self.data["profiles2D"].Bz ** 2 + self.data["profiles2D"].Br ** 2) ** 0.5

        self.data["profiles2D"] = self.data["profiles2D"].assign(Bpol = Bpol2d)

    def save_eqdsk(self, time_index, savedir):
        if self.data is None:
            if self.debug > 0:
                print("<<ERROR>> EfitData: No data available. Read data first.")
                return

        dirpath=os.path.abspath(savedir)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

        pulseNumber = self.data["top"].pulseNumber
        time_value = self.data["top"].time.values[time_index]
        filename = "{}/g_p{}_t{:.5f}".format(dirpath, pulseNumber, time_value)

        # Create or update index.dat file
        if time_index == 0:
            idexfile=open(dirpath+'/index.dat','w')
            idexfile.write('ntimes={}\n'.format(self.data["top"].time.size))
            idexfile.write('path="{}"\n'.format(dirpath))
        else:
            idexfile=open(dirpath+'/index.dat','a')
            idexfile.write('time={:.5f}   filename={}\n'.format(time_value, filename))

        # Close index file for eqdsk
        idexfile.close()

        if self.debug > 0:
            print('----------------------------------------------------------')
            print('Writing EQDSK file')
            print('----------------------------------------------------------')

        with open(filename, 'w') as fp:
            title='EFIT++   '+str(time.strftime("%x"))+'#'+str(pulseNumber)+'-'+str(int(time_value*1000))+'ms'
            print("%40s%4s%8s%4d%4d"%(title,'    ',0,
                                      self.data["regularGrid"].nr.values[0],
                                      self.data["regularGrid"].nz.values[0]),
                  file=fp)

            rMax = self.data["regularGrid"].rMax.values[0]
            rMin = self.data["regularGrid"].rMin.values[0]
            zMax = self.data["regularGrid"].zMax.values[0]
            zMin = self.data["regularGrid"].zMin.values[0]
            rGeom = self.data["separatrixGeometry"].rGeom.values[time_index]
            magneticAxis = self.data["globalParameters"].magneticAxis.values[time_index]

            print("% .9e% .9e% .9e% .9e% .9e"% ( rMax-rMin, zMax-zMin,
                                                 rGeom, rMin, 0.5*(zMax+zMin)), file=fp)
            print("% .9e% .9e% .9e% .9e% .9e"% ( magneticAxis[0], magneticAxis[1],
                                                 self.data["globalParameters"].psiAxis.values[time_index],
                                                 self.data["globalParameters"].psiBoundary.values[time_index],
                                                 self.data["globalParameters"].bvacRgeom.values[time_index] ),
                                                 file=fp)
            print("% .9e% .9e% .9e% .9e% .9e"% ( self.data["globalParameters"].plasmaCurrent.values[time_index],
                                                 self.data["globalParameters"].psiAxis.values[time_index],
                                                 0.0,  magneticAxis[0], 0.0 ), file=fp)
            print("% .9e% .9e% .9e% .9e% .9e"% ( magneticAxis[1], 0.0,
                                                 self.data["globalParameters"].psiBoundary.values[time_index],
                                                 0.0, 0.0), file=fp)
            writeColumns(fp,self.data["fluxFunctionProfiles"].rBphi.values[time_index][:],5,"% .9e")
            writeColumns(fp,self.data["fluxFunctionProfiles"].staticPressure.values[time_index][:],5,"% .9e")
            writeColumns(fp, self.data["fluxFunctionProfiles"].ffPrime.values[time_index][:], 5, "% .9e")
            writeColumns(fp, self.data["fluxFunctionProfiles"].staticPPrime.values[time_index][:], 5, "% .9e")

            psiefit_t=np.transpose(self.data["profiles2D"].poloidalFlux.values[time_index][:])
            writeColumns(fp,psiefit_t,5,"% .9e")
            writeColumns(fp, self.data["fluxFunctionProfiles"].q.values[time_index][:], 5, "% .9e")

            rlimsize = self.data["geometry"]["limiter"].rValues.size
            print("%5d%5d"% ( self.data["separatrixGeometry"].boundaryCoordsDim.size,
                              rlimsize ), file=fp)

            rBoundary = self.data["separatrixGeometry"].rBoundary.values[time_index][:]
            zBoundary = self.data["separatrixGeometry"].zBoundary.values[time_index][:]
            boundaryCoordsDim = self.data["separatrixGeometry"].zBoundary.boundaryCoordsDim.size
            # boundaryCoords.resize(boundaryCoordsDim,1)
            i=0
            bc = np.array( np.zeros(2*boundaryCoordsDim) )
            for v0, v1 in zip(rBoundary, zBoundary):
                 (bc[i], bc[i+1]) = (v0, v1)
                 i += 2
            writeColumns(fp, bc,5,"% .9e")
            i=0
            lm = np.array( np.zeros(2*rlimsize) )
            for v in range(rlimsize):
                 lm[i] = self.data["geometry"]["limiter"].rValues.values[v]
                 lm[i+1] = self.data["geometry"]["limiter"].zValues.values[v]
                 i += 2
            writeColumns(fp,lm,5,"% .9e")

        if self.debug > 0:
            print('----------------------------------------------------------')