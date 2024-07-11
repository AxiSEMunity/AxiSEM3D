import os
import pyvtk
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


# The data structure in element-wise output is too complicated for xarray.open_mfdataset.
# Here we open the files as individual datasets and concatenate them on the variable level.
# This code is compatible with parallel netcdf build (single file output)

# load_wave_data=True:  read wave data and return numpy.ndarray
# load_wave_data=False: do not read wave data and return xarray.DataArray (use False if data is big)

def read_element_output(data_dir, load_wave_data=True):
    ################ open files ################
    # filenames
    nc_fnames = [f for f in os.listdir(data_dir) if 'axisem3d_synthetics.nc' in f]
    print('files to open: ', nc_fnames)

    # open files
    nc_files = []
    for nc_fname in nc_fnames:
        nc_files.append(xr.open_dataset(data_dir + '/' + nc_fname))

    ################ variables that are the same in the datasets ################
    # read Na grid (all azimuthal dimensions)
    na_grid = nc_files[0].data_vars['list_na_grid'].values.astype(int)

    # read time
    data_time = nc_files[0].data_vars['data_time'].values

    ################ variables to be concatenated over the datasets ################
    # define empty lists of xarray.DataArray objects
    xda_list_element_na = []
    xda_list_element_coords = []
    dict_xda_list_element = {}
    dict_xda_data_wave = {}
    for nag in na_grid:
        dict_xda_list_element[nag] = []
        dict_xda_data_wave[nag] = []

    # loop over nc files
    for nc_file in nc_files:
        # append DataArrays
        xda_list_element_na.append(nc_file.data_vars['list_element_na'])
        xda_list_element_coords.append(nc_file.data_vars['list_element_coords'])
        for nag in na_grid:
            dict_xda_list_element[nag].append(nc_file.data_vars['list_element__NaG=%d' % nag])
            dict_xda_data_wave[nag].append(nc_file.data_vars['data_wave__NaG=%d' % nag])

    # concat xarray.DataArray
    xda_list_element_na = xr.concat(xda_list_element_na, dim='dim_element')
    xda_list_element_coords = xr.concat(xda_list_element_coords, dim='dim_element')
    for nag in na_grid:
        dict_xda_list_element[nag] = xr.concat(dict_xda_list_element[nag], dim='dim_element__NaG=%d' % nag)
        dict_xda_data_wave[nag] = xr.concat(dict_xda_data_wave[nag], dim='dim_element__NaG=%d' % nag)

    # read data to numpy.ndarray
    list_element_na = xda_list_element_na.values.astype(int)
    list_element_coords = xda_list_element_coords.values
    dict_list_element = {}
    dict_data_wave = {}
    for nag in na_grid:
        dict_list_element[nag] = dict_xda_list_element[nag].values.astype(int)
        if load_wave_data:
            dict_data_wave[nag] = dict_xda_data_wave[nag].values

    ############### return ################
    if load_wave_data:
        return na_grid, data_time, list_element_na, list_element_coords, dict_list_element, dict_data_wave
    else:
        return na_grid, data_time, list_element_na, list_element_coords, dict_list_element, dict_xda_data_wave








# data dir
data_dir = './output/elements/orthogonal_azimuthal_slices'

# read
na_grid, data_time, list_element_na, list_element_coords, \
dict_list_element, dict_data_wave = read_element_output(data_dir)

# wave dimension to animation
wave_channel = 'Z'
wave_dim = 'RTZ'.index(wave_channel)

# time steps
ntime = len(data_time)

# phi of the slices
phi_slices = np.radians(np.arange(0, 360, 90))
nslice = len(phi_slices)

# GLL coords on elements
nelem = list_element_coords.shape[0]
ngll = list_element_coords.shape[1]
# flattened coords, (s, z)
element_coords_sz = list_element_coords.reshape((nelem * ngll), 2)

# connectivity list, shared by all slices
# with GLL_points_one_edge = [0, 2, 4] in the inparam.output.yaml,
# the GLL point layout on each element is
# 0--1--2
# |  |  |
# 3--4--5
# |  |  |
# 6--7--8
connectivity = []
for ielem in np.arange(nelem):
    start = ielem * 9
    connectivity.append([start + 0, start + 1, start + 4, start + 3])
    connectivity.append([start + 1, start + 2, start + 5, start + 4])
    connectivity.append([start + 3, start + 4, start + 7, start + 6])
    connectivity.append([start + 4, start + 5, start + 8, start + 7])

# loop over slices
for islice, phi in enumerate(phi_slices):
    # create vtk folder
    vtk_dir = data_dir + '/vtk/slice%d' % islice
    os.makedirs(vtk_dir, exist_ok=True)

    # vtk mesh
    xyz = np.ndarray((nelem * ngll, 3))
    xyz[:, 0] = element_coords_sz[:, 0] * np.cos(phi)
    xyz[:, 1] = element_coords_sz[:, 0] * np.sin(phi)
    xyz[:, 2] = element_coords_sz[:, 1]
    vtk_mesh = pyvtk.UnstructuredGrid(list(zip(xyz[:, 0], xyz[:, 1], xyz[:, 2])),
                                      quad=connectivity)

    # loop over elements to read wave data
    wave = np.ndarray((nelem * ngll, ntime))
    for ielem in np.arange(nelem):
        wave[(ielem * ngll):(ielem * ngll + ngll), :] = dict_data_wave[nslice][ielem, islice, :, wave_dim, :]

    # loop over time to write vtk
    for itime in np.arange(ntime):
        vtk = pyvtk.VtkData(vtk_mesh, pyvtk.PointData(pyvtk.Scalars(wave[:, itime], name='U' + wave_channel)))
        vtk.tofile(vtk_dir + '/' + 'wave%d.vtk' % itime, 'binary')
        print('Done time step %d / %d' % (itime + 1, ntime), end='\r')
    print('\nDone slice %d' % (islice + 1))