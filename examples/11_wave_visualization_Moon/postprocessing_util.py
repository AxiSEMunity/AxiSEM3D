#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Util code for the jupyter notebooks
Author: Ceri Nunn, JPL
"""

import os
import csv
import yaml
import numpy as np
import netCDF4 as nc4
from netCDF4 import chartostring
from collections import OrderedDict
import math

# import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
from matplotlib.dates import date2num
from matplotlib.patches import Rectangle
from matplotlib import patheffects

from obspy.core import Stream, Trace, UTCDateTime, Stats
from obspy.core.event import read_events
# from obspy.io.sac import SACTrace
from obspy.signal.detrend import spline
from obspy import read
from obspy import read_inventory
from obspy.signal.filter import envelope
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

from obspy.taup.taup_create import build_taup_model
from obspy.taup import TauPyModel
from obspy.taup.tau import plot_ray_paths
from obspy.taup.utils import parse_phase_list, split_ray_path
from obspy.taup.seismic_phase import SeismicPhase

from irfpy.moon import moon_map

from pdart.util import linear_interpolation, remove_negative_ones

import seaborn

from matplotlib.pyplot import get_cmap

from postprocessing_util_observations import get_station_details


plt.rcParams['figure.figsize'] = [16, 7]

MOON_RADIUS_IN_KM = 1737.1
# MOON_FLATTENING = 1/825
MOON_FLATTENING = 0.0 # flattening for the spherical model

# import pyvtk
# import xarray as xr
from obspy.clients.fdsn.client import Client, FDSNNoDataException

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
# warnings.filterwarnings("ignore", category=FutureWarning)

phase_name_dict = {'kmps' : 'km/s',
        "Sv1s" : 'S (trapped 0-1 km)',
        "Sv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1sSv1s" : 'S (trapped 0-1 km)',
    "Sv12s^1Sv12s": 'S (trapped 1-12 km)',
        "Sv12s^1Sv12s^1Sv12s^1Sv12s^1Sv12s^1Sv12s^1Sv12s^1Sv12s": 'S (trapped 1-12 km)',
    "Sv12s" : 'S (trapped 0-12 km)',
    "Sv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12sSv12s" : 'S (trapped 0-12 km)',
}


### Combined netcdf file
def combined_file(top_dir=None,
        run=None,station_group=None,station_file=None, folder=None, delete=True):

    rank_file = os.path.join(top_dir, run, folder, 'output', 'stations', station_group,'rank_station.info')

    # read rank-station info
    rank_station_info = np.loadtxt(rank_file, dtype=str, skiprows=1)

    # deal with slightly weird case when there is only one row
    if len(rank_station_info.shape) == 1:
        rank_station_info = [[rank_station_info[0], rank_station_info[1], rank_station_info[2]]]

    # dict: mpi-rank -> [station keys]
    rank_station_dict = {}

    # (lat, lon) of stations re-ordered by data
#     stlatlon_in_data_order = []

    ranks = []
    for item in rank_station_info:
        rank = item[0]
        ranks.append(int(rank))


    ordered_ranks = list(OrderedDict.fromkeys(ranks))


    dim_station_str_length = 0

    list_of_ds = []
    station_numbers = []
    for mpi_rank in ordered_ranks:
        fn_file1 = '{}{}'.format('axisem3d_synthetics.nc.rank',str(mpi_rank))
        fn_file = os.path.join(top_dir, run, folder, 'output', 'stations', station_group, fn_file1)

        ds = nc4.Dataset(fn_file)

        list_of_ds.append(ds)
        data_wave = ds.variables['data_wave']
        list_channel = ds.variables['list_channel']
        list_station = ds.variables['list_station']
        if dim_station_str_length < list_station.shape[1]:
            dim_station_str_length = list_station.shape[1]

#         print('This is list_station',list_station.shape )
        station_numbers.append(data_wave.shape[0])

    # get the dimensions (from the last file in the list)

    _, dim_channel, dim_time = data_wave.shape
#
#     print(dim_channel, dim_time)
    _, dim_channel_str_length = list_channel.shape

# #     print(ds['dim_station'].dimensions)
#     print(dim_channel_str_length,dim_station_str_length)
#     return
    i_station = 0
#     print(str(chartostring(list_station[i_station])))


#     print('data shape ', ds['data_wave'].shape)

#     for i in range(0,len(list_channel)):

#         print('Channel ', str(chartostring(list_channel[i])))
#         print(data_wave[i_station, i, 0:10])
#         print(data_wave[i_station, i, -10:])
#         print(abs(data_wave[i_station, i, :]).max())


#     42, 3, 51764

    # make a combined file -
    file1 = 'axisem3d_synthetics.nc.rank_all.nc'
    file = os.path.join(top_dir, run, folder, 'output', 'stations', station_group, file1)


    if os.path.exists(file):
        os.remove(file)

    f = nc4.Dataset(file,'w', format='NETCDF4', clobber=True) #'w' stands for write

    total_stations = sum(station_numbers)

    f.createDimension('dim_time', dim_time)
    f.createDimension('dim_station', total_stations)
    f.createDimension('dim_channel', dim_channel)
    f.createDimension('dim_channel_str_length', dim_channel_str_length)
    f.createDimension('dim_station_str_length', dim_station_str_length)

    data_time = f.createVariable('data_time', 'float64', 'dim_time')
    data_wave = f.createVariable('data_wave', 'float32', ('dim_station', 'dim_channel', 'dim_time'))
    list_channel = f.createVariable('list_channel', 'S1', ('dim_channel', 'dim_channel_str_length'))
    list_station = f.createVariable('list_station', 'S1', ('dim_station', 'dim_station_str_length'))
#     print('dimensions', total_stations, dim_station_str_length)

#              original_stream.append(Trace(data_wave[rank,ich,:], header=stats))

    i_start = 0
    for i, ds in enumerate(list_of_ds):

        i_end = i_start + station_numbers[i]

        data_time[:] = ds.variables['data_time'][:]
        data_wave[i_start:i_end,:,:] = ds.variables['data_wave'][:,:,:]
        list_channel[:,:] = ds.variables['list_channel'][:,:]
#         print(ds.variables['list_station'].shape)
#         print(list_station.shape)
#         print('before')
        list_station_shape = ds.variables['list_station'].shape[1]
#         list_station[i_start:i_end,0:list_station_shape] = ds.variables['list_station'][:,:]
        list_station[i_start:i_end,0:list_station_shape] = ds.variables['list_station'][:,:]
#         print(list_station)
#         print(i_start,i_end)

        i_start = i_end


    # read rank-station info
    rank_station_info = np.loadtxt(rank_file, dtype=str, skiprows=1)

    station_file1 = os.path.join(top_dir, run, folder, 'input',station_file)
    print(station_file1)

    list_lat = f.createVariable('list_lat', 'float32', ('dim_station'))
    list_lon = f.createVariable('list_lon', 'float32', ('dim_station'))

#     print(str(chartostring(list_station[i_station])))

    # read station locations
    sf = np.loadtxt(station_file1, dtype=str, skiprows=3)


    # deal with slightly weird case when there is only one row
    if len(sf.shape) == 1:
        sf_list = sf[0], sf[1], sf[2], sf[3], sf[5]
        sf1 = np.array(sf_list)
        sf = np.array([sf1])

    # dict: station key -> [lat, lon]
    stlatlon_dict = {}
    for ist in np.arange(total_stations):
        key = sf[ist, 1] + '.' + sf[ist, 0]
        stlatlon_dict[key] = np.array([float(sf[ist, 2]), float(sf[ist, 3])])

#         print(key, stlatlon_dict[key])



    for ist in np.arange(total_stations):
        key_nc = str(chartostring(list_station[ist]))
#         print(key_nc)

        list_lat[ist] = stlatlon_dict[key_nc][0]
        list_lon[ist] = stlatlon_dict[key_nc][1]

    f.close()

    print('made combined file ', file)

    # read it back in, to check it is valid
    f_read = nc4.Dataset(file,'r', format='NETCDF4') #'r' stands for read
    for mpi_rank in ordered_ranks:
        fn_file1 = '{}{}'.format('axisem3d_synthetics.nc.rank',str(mpi_rank))
        fn_file = os.path.join(top_dir, run, folder, 'output', 'stations', station_group, fn_file1)
        if delete:
            os.remove(fn_file)
    fn_file = os.path.join(top_dir, run, folder, 'output', 'stations', station_group, 'rank_station.info')
    if delete:
        os.remove(fn_file)
        print('deleted original files ')


def combined_file_slices(top_dir=None,
        run=None,element_name=None,folder=None,delete=True):

    data_dir = os.path.join(top_dir,run,folder,'output','elements',element_name)

    # filenames
    nc_fnames = [f for f in os.listdir(data_dir) if 'axisem3d_synthetics.nc' in f and 'rank_all.nc' not in f]

    element_numbers = []
    for i, nc_frame in enumerate(nc_fnames):
        nc_file = os.path.join(data_dir,nc_frame)
#         print(nc_file)

        ds = nc4.Dataset(nc_file)
        list_element_na = ds.variables['list_element_na']
        dim_element, _ = list_element_na.shape
        element_numbers.append(dim_element)

#         print(ds.variables['list_na_grid'][:])
#         print(ds.variables['list_element_na'][:])

    total_elements = sum(element_numbers)

    print(total_elements)

    i_start = 0
    for i, nc_frame in enumerate(nc_fnames):
        nc_file = os.path.join(data_dir,nc_frame)

        source_ds = nc4.Dataset(nc_file)

#         print(nc_file)
#         for name, dimension in source_ds.dimensions.items():

#             print(name, dimension)

        if i == 0:

            # make a combined file -
            file1 = 'axisem3d_synthetics.nc.rank_all.nc'
            file = os.path.join(top_dir,run,folder,'output','elements',element_name, file1)

            if os.path.exists(file):
                os.remove(file)

            dest_ds = nc4.Dataset(file,'w', format='NETCDF4', clobber=True) #'w' stands for write
#             # copy global attributes all at once via dictionary
#             dest_ds.setncatts(source_ds.__dict__)

            for name, dimension in source_ds.dimensions.items():
                if 'dim_na__NaG' in name:
                    dim_na_number = str(dimension.size)

            dim_element__NaG_input = 'dim_element__NaG=' + dim_na_number
            list_element__NaG_input = 'list_element__NaG=' + dim_na_number
            data_wave__NaG_input = 'data_wave__NaG=' + dim_na_number
            dim_na__NaG_input = 'dim_na__NaG=' + dim_na_number

            # copy dimensions
            for name, dimension in source_ds.dimensions.items():
                if name in ['dim_element']:
                    dest_ds.createDimension(name,total_elements)
                elif name in [dim_element__NaG_input]:
                    dest_ds.createDimension('dim_element__NaG',total_elements)
                elif name in [dim_na__NaG_input]:
                    dest_ds.createDimension(
                        'dim_na__NaG', (len(dimension) if not dimension.isunlimited() else None))
                else:
                    dest_ds.createDimension(
                        name, (len(dimension) if not dimension.isunlimited() else None))

#
            # Copy variables
            exclude = [list_element__NaG_input, 'list_element_na', 'list_element_coords', data_wave__NaG_input]
            for v_name, varin in source_ds.variables.items():
                if v_name not in exclude:
                    outVar = dest_ds.createVariable(v_name, varin.datatype, varin.dimensions)

                    # Copy variable attributes
                    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
                    outVar[:] = varin[:]

            # create the variables which will be concatenated
            list_element__NaG = dest_ds.createVariable('list_element__NaG', 'int32', 'dim_element__NaG')
            list_element_na = dest_ds.createVariable('list_element_na', 'int32', ('dim_element', 'dim_5'))
            list_element_coords = dest_ds.createVariable('list_element_coords', 'float64', ('dim_element', 'dim_GLL', 'dim_2'))
            data_wave__NaG = dest_ds.createVariable('data_wave__NaG', 'float32', (
                  'dim_element__NaG', 'dim_na__NaG', 'dim_GLL', 'dim_channel', 'dim_time'))



        # for all files, copy the data for the variables that need concatenating

        # concatenate the data
        i_end = i_start + element_numbers[i]
        print(i_start,i_end)
        list_element__NaG[i_start:i_end] = source_ds.variables[list_element__NaG_input][:]
        list_element_na[i_start:i_end,:] = source_ds.variables['list_element_na'][:,:]
        list_element_coords[i_start:i_end,:,:] = source_ds.variables['list_element_coords'][:,:,:]
#         print(source_ds.variables['data_wave__NaG=1'][:,:,:,:,:])
#         print(len(source_ds.variables['data_wave__NaG=1'][:,:,:,:,:]))
        data_wave__NaG[i_start:i_end,:,:,:,:] = source_ds.variables[data_wave__NaG_input][:,:,:,:,:]
        i_start = i_end


    print('Source')
    print(source_ds)
    print('Destination')
    print(dest_ds)

    print('made combined file')
    print(file)


    dest_ds.close()

    # read it in again to check it has worked
    dest_ds_read = nc4.Dataset(file,'r', format='NETCDF4') #'r' stands for read
    for i, nc_frame in enumerate(nc_fnames):
        nc_file = os.path.join(data_dir,nc_frame)
        if delete:
            os.remove(nc_file)
            print('deleted original files')

### Plot locations

def plot_locations(top_dir=None,
        run=None,station_group=None,station_file=None,folder=None,source_name=None):

    station_file = os.path.join(top_dir, run, folder, 'input', station_file)

    # read station info
    station_file1 = np.loadtxt(station_file, dtype=str, skiprows=3)

    ####################################
    # draw a map of event and stations #
    ####################################

    print('preparing the map ...')

    plt.figure(dpi=150)
    # draw map
    map1 = Basemap(projection='cyl', resolution='l', lon_0=0)

    # prepare the map
    map = moon_map.MoonMapSmall()
    map_on_sphere = map.gridsphere()
    blon, blat = map_on_sphere.get_bgrid()
    level = map_on_sphere.get_average()

    # pcolor on the map
    map1.pcolormesh(blon, blat, level, latlon=True, cmap='gray')

    # read event location
    location_file = os.path.join(top_dir, run, folder, 'input','inparam.source.yaml')
    with open(location_file, 'r') as file:
        source_yaml = yaml.load(file, Loader=yaml.FullLoader)
    loc_leaf = source_yaml['list_of_sources'][0][source_name]['location']
    event_latlon = loc_leaf['latitude_longitude']
    event_depth = loc_leaf['depth']

    # draw event
    map1.scatter(event_latlon[1], event_latlon[0], latlon=True,
                s=150, c='r', marker='*', lw=0, zorder=100)
    # draw stations
    map1.scatter(station_file1[:, 3].astype(float), station_file1[:, 2].astype(float), latlon=True,
                s=30, c='b', marker=7, lw=0, zorder=10)
#     map1.scatter(station_file1[-5:, 3].astype(float), station_file1[-5:, 2].astype(float), latlon=True,
#                 s=30, c='r', marker=7, lw=0, zorder=10)
    plt.show()




### Get all the streams
def get_all_streams_from_netcdf(top_dir,run,short_title,station_group,station_file,
  source_name,show_seismogram=False,folder=None,model_taup=None,model_type='MOON'):
    rank_file = os.path.join(top_dir, run, folder, 'output', 'stations', station_group,'rank_station.info')

    # read a combined file
    fn_file = os.path.join(top_dir, run, folder, 'output', 'stations', station_group, 'axisem3d_synthetics.nc.rank_all.nc')

    ds = nc4.Dataset(fn_file)

    list_station = []
    list_station1 = ds.variables['list_station']
    for n in range(list_station1.shape[0]):
        station = str(chartostring(list_station1[n]))
        list_station.append(station)


    data_time = ds.variables['data_time']
    data_wave = ds.variables['data_wave']
    list_lat = ds.variables['list_lat']
    list_lon = ds.variables['list_lon']
    # print(data_wave.dimensions)
    # float32 data_wave(dim_station, dim_channel, dim_time)
    # current shape = (27, 3, 32635)
#     print(ds.variables['data_wave'])

    st_rank_all = {}
    list_station = []
    list_station1 = ds.variables['list_station']
    for n in range(list_station1.shape[0]):
        key = str(chartostring(list_station1[n]))
        list_station.append(key)
        st_rank_all[key] = n

    # read event location
    file_source = os.path.join(top_dir,run,folder,'input','inparam.source.yaml')
    with open(file_source, 'r') as file:
        source_yaml = yaml.load(file, Loader=yaml.FullLoader)
    loc_leaf = source_yaml['list_of_sources'][0][source_name]['location']
    event_latlon = loc_leaf['latitude_longitude']
    event_latitude = event_latlon[0]
    event_longitude = event_latlon[1]
    event_depth = loc_leaf['depth']
    event_depth = float(event_depth)

    # trace header
    stats = Stats()
    stats.starttime = UTCDateTime(data_time[0])
    stats.delta = UTCDateTime(data_time[1] - data_time[0])
    stats.npts = len(data_time)

    if model_type == 'MOON':
        radius_in_km = MOON_RADIUS_IN_KM
        flattening = MOON_FLATTENING
    elif model_type == 'EARTH':
        radius_in_km = EARTH_RADIUS_IN_KM
        flattening = EARTH_FLATTENING
    else:
        print('Model type not found')

    original_stream = Stream()

    for key in st_rank_all.keys():

        # get the rank for each station
        rank = st_rank_all[key]

        station_latitude = list_lat[rank]
        station_longitude = list_lon[rank]

        distance, azimuth_A_B, azimuth_B_A =  gps2dist_azimuth(
          event_latitude, event_longitude, station_latitude, station_longitude, radius_in_km*1000., flattening)

        distance_in_km = distance/1000.
        event_depth = event_depth/1000.

        distance_in_degree = kilometers2degrees(distance_in_km,radius=radius_in_km)

        stats.distance_in_km = distance_in_km
        stats.distance_in_degree = distance_in_degree

        stats.run = run
        stats.short_title = short_title

        if model_taup is not None:
            arrivals = model_taup.get_travel_times(source_depth_in_km=0, distance_in_degree=distance_in_degree, phase_list=(['P']))
            if len(arrivals) > 0:
                stats.P_arrival = arrivals[0].time
            arrivals = model_taup.get_travel_times(source_depth_in_km=0, distance_in_degree=distance_in_degree, phase_list=(['S']))
            if len(arrivals) > 0:
                stats.S_arrival = arrivals[0].time
            arrivals = model_taup.get_travel_times(source_depth_in_km=0, distance_in_degree=distance_in_degree, phase_list=(['PP']))
            if len(arrivals) > 0:
                stats.PP_arrival = arrivals[0].time


        list_channel = []
        list_channel1 = ds.variables['list_channel']
        for n in range(list_channel1.shape[0]):
            channel = str(chartostring(list_channel1[n]))
            list_channel.append(channel)

        list_channel1 = []
        for ch in list_channel:
            list_channel1.append(channel_dict[ch])

#         if list_channel == ['U1','U2','U3']:
#             list_channel = ['R','T','Z']

        stats.network, stats.station = key.split('.')

#         print(list_channel)

    #     numpy.where(condition[, x, y])
        # note that rank is using the combined rank

        for ich, ch in enumerate(list_channel1):
#             print(ich, ch)
            stats.channel = ch
            # latitude/longitude are extra attributes which will not be
            # saved if the stream is saved to the drive
            stats.latitude = list_lat[rank]
            stats.longitude = list_lon[rank]

            original_stream.append(Trace(data_wave[rank,ich,:], header=stats))
#             if rank == 10:
#                 print('Data Wave ', ch, data_wave[rank,ich,0:10])
#                 print('Data Wave ', ch, data_wave[rank,ich,-10:])
#             if ch == 'T':
#                 print('Min/Max', rank, ch, data_wave[:,ich,:].min(),data_wave[:,ich,:].max())

    for ich, ch in enumerate(list_channel1):
        print('Min/Max', ch, data_wave[:,ich,:].min(),data_wave[:,ich,:].max())


    return original_stream, st_rank_all


phase_list_colors = {"P" : '#1f78b4', "PP": '#33a02c', "PcP": '#e31a1c',
                     "Pdiff" : '#ff7f00', "PvmP": '#6a3d9a', "Pvmp": '#6a3d9a',
                     "S" : '#710193', 'SS' : 'cyan', "PS" : '#BE9394',
                     "Pv9P" : 'green', "Pv9S" : 'blue',  "Pv20P" : 'purple', "Pv20S": 'orange',  "Pv25P": 'black', "Pv25S": 'yellow',
                    "2kmps": 'midnightblue', "Pv50P" : 'green', "Pv50S" : 'blue',
                    }


### Plot epicentral distance with taup
def plot_spectrogram(original_trace1,title=None,startsecond=0,endsecond=1800,clip=[0,1]):

    original_trace1_copy = original_trace1.copy()


    # print('{:} Distance from Source {:.1f} km'.format(station, original_trace1_copy.stats.distance_in_km))
#     original_trace1_copy.resample(lowpass)



    if endsecond is not None:
        original_trace1_copy = original_trace1_copy.slice(endtime=UTCDateTime(endsecond))
    if startsecond is not None:
        original_trace1_copy = original_trace1_copy.slice(starttime=UTCDateTime(startsecond))

    log = True
    fig = original_trace1_copy.spectrogram(show=False, log=log, title=title,wlen=10, samp_rate=original_trace1_copy.stats.sampling_rate,clip=clip)
    ax = fig.axes[0]
    if log:
        mappable = ax.collections[0]
    else:
        mappable = ax.images[0]
    plt.colorbar(mappable=mappable, ax=ax,
#                  label='??'
                )
    plt.show()

    original_trace1_copy.plot()
    plt.show()

def plot_different_filtering(original_trace1,original_trace2=None,
        taup_show=True,model_taup=None,title=None,
        freqs=[('lowpass',0.5), (0.01,0.02),(0.02,0.03),(0.03,0.04),
        (0.04,0.05),(0.05,0.06),(0.06,0.07),(0.07,0.08),(0.08,0.09),
        (0.09,0.1),(0.1,0.2),(0.2,0.3),(0.3,0.4),(0.4,0.5)],
        startsecond=0,endsecond=1800,normalize='relative',scale=1,figsize=(11, 17)):

    # The normalization factor can be 'all' (the trace is lowpassed, and the max amplitude is calculated)
    # or 'relative', where each filtered trace is normalized.
    # Relative is better for seeing small amplitude phases, but 'all' is better for getting a realistic
    # sense of which frequencies have power.
    # The default is 'relative'
    if normalize not in ['relative','all']:
        print("Please select 'relative' or 'all' for the normalization factor.")
        return

    taper_len = 10

    normalize_text = None

    print('event_depth=0')
    event_depth=0

    if taup_show:
        distance_in_degree = original_trace1.stats.distance_in_degree
        arrivals = model_taup.get_ray_paths(source_depth_in_km=event_depth,
          distance_in_degree=original_trace1.stats.distance_in_degree,
          phase_list=['P'])
        print('arrivals ', original_trace1.stats.distance_in_km, arrivals)

    # if the normalization factor is all, get it first from the lowpass filter version
    if normalize == 'all':
        for i, freq in enumerate(freqs):


            if freq[0] == 'lowpass':
                process_trace_freq1 = original_trace1.copy()

                # taper the right side of the trace with a 10 second taper
                process_trace_freq1 = process_trace_freq1.taper(max_percentage=None,max_length=taper_len,side='right')

                process_trace_freq1.filter('lowpass', freq=freq[1],zerophase=True)

                # slice the trace if required
                orig_endtime = process_trace_freq1.stats.endtime
                print('orig ', orig_endtime)
                if endsecond is not None:
                    process_trace_freq1 = process_trace_freq1.slice(endtime=UTCDateTime(endsecond))
                    print('end ', process_trace_freq1.stats.endtime)
                if process_trace_freq1.stats.endtime > orig_endtime - taper_len:
                    process_trace_freq1 = process_trace_freq1.slice(endtime=orig_endtime - taper_len)
                    print('end ', process_trace_freq1.stats.endtime)


                if startsecond is not None:
                    process_trace_freq1 = process_trace_freq1.slice(starttime=UTCDateTime(startsecond))
                print('start ', process_trace_freq1.stats.starttime)
#                 return

                # find the normalization factor
                max_data = abs(process_trace_freq1.data).max()

     # a second, optional trace

    fig, ax = plt.subplots(figsize=figsize,dpi=300)
    for i, freq in enumerate(freqs):
        process_trace_freq1 = original_trace1.copy()

        # taper the right side of the trace with a 10 second taper
        process_trace_freq1 = process_trace_freq1.taper(max_percentage=None,max_length=taper_len,side='right')

        if freq[0] == 'lowpass':
            process_trace_freq1.filter('lowpass', freq=freq[1],zerophase=True)
            annotation = '< {:.3f} Hz'.format(freq[1])
        else:
            process_trace_freq1.filter('bandpass', freqmin=freq[0], freqmax=freq[1],zerophase=False)
            annotation = '{:.3f}-{:.3f} Hz'.format(freq[0],freq[1])

        # resample introduces periodic signals to these data - if required, try decimation
        # process_trace_freq1.resample(sampling_rate=sampling_rate)
        # process_trace_freq1.decimate(2,no_filter=True)

        # slice the trace if required
        orig_endtime = process_trace_freq1.stats.endtime
        if endsecond is not None:
            process_trace_freq1 = process_trace_freq1.slice(endtime=UTCDateTime(endsecond))
        if process_trace_freq1.stats.endtime > orig_endtime - taper_len:
            process_trace_freq1 = process_trace_freq1.slice(endtime=orig_endtime - taper_len)

        if startsecond is not None:
            process_trace_freq1 = process_trace_freq1.slice(starttime=UTCDateTime(startsecond))

        if normalize != 'all':
            # find the normalization factor
            max_data = abs(process_trace_freq1.data).max()

        # normalize the trace
        process_trace_freq1.data = (process_trace_freq1.data/max_data*scale) + i
#         print(process_trace_freq1.times())
#         print(process_trace_freq1.times(reftime=UTCDateTime(0)))
#         print('trace ', process_trace_freq1.times(reftime=UTCDateTime(0)[0]))
        plt.plot(process_trace_freq1.times(reftime=UTCDateTime(0)), process_trace_freq1.data,color='b')

        if taup_show:
            for arrival in arrivals:
                color = phase_list_colors.get(arrival.name,'r')

                plt.vlines(x=arrival.time, ymin=i-scale/4, ymax=i+scale/4, color=color, label=arrival.name)
#                     if i == 0:
#                         ax.annotate(arrival.name, (arrival.time,ax.get_ylim()[1]), color=color)
                ax.annotate(text=annotation, xy=(endsecond*1.01,i), xycoords='data',
                            horizontalalignment='left', verticalalignment='center',fontsize=12, color='k',annotation_clip=False)

    ax.annotate(text=normalize_text, xy=(0.99,0.01), xycoords='axes fraction',
                horizontalalignment='right', verticalalignment='bottom',fontsize=10, color='k')

    plt.xlabel('Time [s]')
    legend_without_duplicate_labels(ax)
#     plt.xlim(process_trace_freq1.times(reftime=UTCDateTime(0))[0],process_trace_freq1.times(reftime=UTCDateTime(0))[-1])
    plt.xlim(startsecond,endsecond)
#     ax.set_yticklabels([])
#     print(freq)
    if title is not None:
        plt.title(title)
#     plt.title('Propagation of an Explosive Source\nLowpass filtered below {:.1f} Hz'.format(freq))
    plt.show()
    plt.close()

# def get_cmap_phaselist(n, name='rainbow'):
#     # See get_phase_colors(phase_list) - it's more consistent
#     '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
#     RGB color; the keyword argument name must be a standard mpl colormap name.'''
#     return plt.cm.get_cmap(name, n)

# phase_list_dict = {'P': 'r', 'S': 'b'}

## Plot the epicentral distance with Taup
# plot taup and the simulations or observations - the plot shows the seismograms vertically by default
def plot_epicentral_distance_taup(original_stream,inv=None,seismogram_show=True,
    model_taup=None,model_taup_label=None,taup_show=True,just_taup=False,title=None,
    freqmin=None,freqmax=None,channel='Z',obs_start=-1800,startsecond=0,endsecond=1800,normalize='relative',scale=20,scale_list=None,taup_height=100,
    observations=True,taper_len=10,phase_list=["P", "PP", "PcP", "Pdiff", "PvmP"],source_depth=0.0,degree_min=0,degree_max=180,
    catalogs=None,
    raw=False,trace_title_show=True,show_impact_type=False,
#     pre_filt_env = [[0.1,0.25,0.75,1],[0.5,0.75,1.25,1.5],[1,1.25,1.75,2.25]]):
#     pre_filt_env = [[0.3,0.4,0.5,0.6]]
    pre_filt_env=None,smooth_periods=10,fill_envelope=False, save_fig=False,figsize=(11, 17),legend_loc='best', ncol=2, show_legend=True, seismograms_vertical=True):

    if model_taup_label is not None:
        print('Calculated using model : ', model_taup_label)
        taup_label = model_taup_label.replace('_taup','')
        taup_label = "Seismic phases: {}".format(taup_label)
    else:
        taup_label = ''

    epicentral_stream = Stream()

    print(channel)
    print(original_stream)

# select tshe streams to work on
# working on the epicentral stream, remove the instrument response, remove the mean if necessary
# and taper if necessary
# plot the epicentral stream


    for tr in original_stream.select(channel=channel):
#         if tr.stats.latitude == 0.0 and tr.stats.longitude >= 0.0 and tr.stats.longitude < 60:
        tr1 = tr.copy()
        'Lunar Model M1 with heterogeneity '

#     remove negative_ones

#     remove the mean

#     can interpolate across gaps of only one sample - but only relevant for 'SHZ'

#     add the linear interpolation  - but it needs to not be a mask

    # make some basic preparations for the traces (remove the -1 values, remove the single sample gaps, taper, remove the mean)
    if observations:
        remove_negative_ones(epicentral_stream)
        for tr in epicentral_stream:
            # interpolate across the gaps of one sample
            linear_interpolation(tr,interpolation_limit=1)
            # if it's an observational trace, and raw, remove the mean
            if raw:
                av = tr.data.mean()
                tr.data = tr.data - av


    # taper
    if observations:
        epicentral_stream.taper(max_percentage=None,max_length=taper_len,side='both')
    else:
        epicentral_stream.trim(starttime=UTCDateTime(-200.0),pad=True,fill_value=0.0)
        epicentral_stream.taper(max_percentage=None,max_length=taper_len,side='right')


    if pre_filt_env is not None:
        epicentral_stream_env = epicentral_stream.copy()

    if seismogram_show:
        # process the trace
        if observations:
            if raw == False:

                if inv is None:
                    client = Client("IRIS")
                    # get the response file (wildcards allowed)
                    inv = client.get_stations(starttime=UTCDateTime('1969-01-01'),endtime=UTCDateTime('1977-09-30T23:59:59'),
                        network='XA', sta='*', loc='*', channel='*',
                        level="response")

                for tr in epicentral_stream:
                # for removing the instrument response from a seimogram,
                # it is useful to get a mask, then interpolate across the gaps,
                # then mask the trace again.

                    if tr.stats.channel in ['MH1', 'MH2', 'MHZ','SHZ']:
                        # add linear interpolation but keep the original mask
                        original_mask = linear_interpolation(tr,interpolation_limit=None)

#                         # apply a bandpass with a zerophase filter (AFTER interpolating across the gaps)
#                         tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)

                        # remove the instrument response
                        pre_filt = [freqmin/1.5,freqmin,freqmax,freqmax*1.5]
                        tr.remove_response(inventory=inv, pre_filt=pre_filt, zero_mean=True, taper=True, output="DISP",
                                   water_level=None, plot=False)
                        # apply the mask back to the trace
                        tr.data = np.ma.masked_array(tr, mask=original_mask)

        else:
            # bandpass the simulated trace if required
            if raw == False:
#                 epicentral_stream.plot()
                epicentral_stream.filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=False)

#                 epicentral_stream.filter('lowpass', freq=freqmax,zerophase=False)
#                 epicentral_stream.trim(starttime=UTCDateTime(-40.0),endtime=UTCDateTime(40.0)).plot()
#                 epicentral_stream.plot()
#                 return

#                 from obspy.signal.invsim import (cosine_taper, cosine_sac_taper,
#                                                  invert_spectrum)
#                 # remove the instrument response
#                 pre_filt = [freqmin/1.5,freqmin,freqmax,freqmax*1.5]
#                 freq_domain_taper = cosine_sac_taper(freqs, flimit=pre_filt)
#                 data *= freq_domain_taper

    for i, tr in enumerate(epicentral_stream):
        if observations:
            tr.trim(starttime=tr.stats.impact_time + startsecond, endtime=tr.stats.impact_time + endsecond)
        else:
            tr.trim(starttime=UTCDateTime(startsecond), endtime=UTCDateTime(endsecond))
            print(tr.stats.station, tr.stats.distance_in_km)

        if normalize == 'relative':
            # find the normalization factor
            max_data = abs(tr.data).max()
            # normalize the trace
            tr.data = tr.data/max_data
        elif normalize == 'P_arrival':
            if 'P_arrival' in tr.stats.keys():
                tr_local = tr.copy()
                tr_local.trim(starttime=UTCDateTime(tr.stats.P_arrival)-10.,endtime=UTCDateTime(tr.stats.P_arrival)+30.)
                # find the normalization factor based on -10 and +60 seconds after the P arrival
                if len(tr_local.data) > 0:
                    # calculate the normalization factor on the short section, but apply to longer trace
                    max_data = abs(tr_local.data).max()
                    tr.data = tr.data/max_data
                else:
                    print('P_arrival trace too small - removing the trace', tr)
                    epicentral_stream.remove(tr)
            else:
                print('P_arrival not found - removing the trace', tr)
                epicentral_stream.remove(tr)
        elif normalize == 'S_arrival':
            if 'S_arrival' in tr.stats.keys():
#                 print('Before ', tr.stats.starttime, tr.stats.endtime)
                tr.trim(starttime=UTCDateTime(tr.stats.S_arrival)-30.,endtime=UTCDateTime(tr.stats.S_arrival)+70.)
#                 tr.trim(starttime=tr.stats.S_arrival-40.)
#                 print('After', tr.stats.starttime, tr.stats.endtime, tr.stats.endtime-tr.stats.starttime )
                # find the normalization factor based on -60 and +60 seconds after the P arrival
                if len(tr.data) > 0:
                    max_data = abs(tr.data).max()
                    tr.data = tr.data/max_data
                else:
                    print('S_arrival trace too small - removing the trace', tr)
                    epicentral_stream.remove(tr)
            else:
                print('S_arrival not found - removing the trace', tr)
                epicentral_stream.remove(tr)

        elif normalize == 'PP_arrival':
            # this might be the best way to normalize
            if 'PP_arrival' in tr.stats.keys():
#                 print('Before ', tr.stats.starttime, tr.stats.endtime)
                tr_local = tr.copy()
                tr_local.trim(starttime=UTCDateTime(tr.stats.PP_arrival)-30.,endtime=UTCDateTime(tr.stats.PP_arrival)+70.)
#                 tr.trim(starttime=tr.stats.S_arrival-40.)
#                 print('After', tr.stats.starttime, tr.stats.endtime, tr.stats.endtime-tr.stats.starttime )
                # find the normalization factor based on -30 and +70 seconds after the P arrival
                if len(tr_local.data) > 0:
                    # calculate the normalization factor on the short section, but apply to longer trace
                    max_data = abs(tr_local.data).max()
                    tr.data = tr.data/max_data
                else:
                    print('PP_arrival trace too small - removing the trace', tr)
                    epicentral_stream.remove(tr)
            else:
                print('PP_arrival not found - removing the trace', tr)
                epicentral_stream.remove(tr)
        elif normalize == 'first_arrivals':
            if 'P_arrival' in tr.stats.keys() and 'S_arrival' in tr.stats.keys():
#                 print('Before ', tr.stats.starttime, tr.stats.endtime)
                tr.trim(starttime=UTCDateTime(tr.stats.P_arrival)-10.,endtime=UTCDateTime(tr.stats.S_arrival)+30.)
#                 tr.trim(starttime=tr.stats.S_arrival-40.)
#                 print('After', tr.stats.starttime, tr.stats.endtime, tr.stats.endtime-tr.stats.starttime )
                # find the normalization factor based on -60 and +60 seconds after the P arrival
                if len(tr.data) > 0:
                    max_data = abs(tr.data).max()
                    tr.data = tr.data/max_data
                else:
                    print('first_arrivals trace too small - removing the trace', tr)
                    epicentral_stream.remove(tr)
            else:
                print('S_arrival not found - removing the trace', tr)
                epicentral_stream.remove(tr)
        elif normalize == 'first_arrivals_longer':
            if 'P_arrival' in tr.stats.keys() and 'S_arrival' in tr.stats.keys():
#                 print('Before ', tr.stats.starttime, tr.stats.endtime)
                tr_local = tr.copy()
                tr_local.trim(starttime=UTCDateTime(tr.stats.P_arrival)-10.,endtime=UTCDateTime(tr.stats.S_arrival)+30.)
#                 tr.trim(starttime=tr.stats.S_arrival-40.)
#                 print('After', tr_local.stats.starttime, tr_local.stats.endtime, tr_local.stats.endtime-tr_local.stats.starttime )
                # also trim the trace to a slightly longer
                tr.trim(starttime=UTCDateTime(tr.stats.P_arrival)-10.,endtime=UTCDateTime(tr.stats.S_arrival)+30.+30.)
                # find the normalization factor based on -60 and +60 seconds after the P arrival
                if len(tr_local.data) > 0:
                    # calculate the normalization factor on the short section, but apply to longer trace
                    max_data = abs(tr_local.data).max()
                    tr.data = tr.data/max_data
                else:
                    print('first_arrivals_all trace too small - removing the trace', tr)
                    epicentral_stream.remove(tr)
            else:
                print('S_arrival not found - removing the trace', tr)
                epicentral_stream.remove(tr)
        if scale_list is not None:
            tr.data = tr.data*scale_list[i]*scale
        else:
            tr.data = tr.data*scale


    # cmap_phaselist = get_cmap_phaselist(len(phase_list))




    print('event_depth=0')
    event_depth=0

    annotate_time = (endsecond - 1.07*(endsecond-startsecond))
    annotate_time2 = (endsecond + 0.02*(endsecond-startsecond))

#     fig = plt.figure(dpi=150)
    fig, ax = plt.subplots(figsize=figsize,dpi=300)




    nasa_red = '#FC3D21'
    nasa_blue = '#1B3C8C'
    bright_blue = '#1B3CF8'
    blue_grey =  '#454756'


    alpha = 0.6
    if seismogram_show:
#     plot seismograms
        for i in range(len(epicentral_stream)):
            if i == 0:
                print('CURRENT RUN:', epicentral_stream[i].stats.run)

    #     plt.gca().axis('off')



            if observations:
                if epicentral_stream[i].stats.channel == 'MH1':
                    color = blue_grey
                else:
                    color = bright_blue
#                 elif (i % 2) == 0:
#                     color = nasa_blue
#                 else:
#                     color = bright_blue
            else:
                if (i % 2) == 0:
                    color = blue_grey
                else:
                    color = 'black'


            if observations:
                ref_time = epicentral_stream[i].stats.impact_time
            else:
                ref_time = UTCDateTime(0)

            time = epicentral_stream[i].times(reftime = ref_time)
            data = epicentral_stream[i].data + epicentral_stream[i].stats.distance_in_degree

            observation_label = '{} ({:.2f}-{:.2f} Hz)'.format(epicentral_stream[i].stats.channel,freqmin,freqmax)

            # plot the data
            if seismograms_vertical:
                plt.plot(data, time,color=color,alpha=alpha,label=observation_label,zorder=8,linewidth=0.5)
            else:
                plt.plot(time, data,color=color,alpha=alpha,label=observation_label,zorder=8,linewidth=0.5)


        # plot tmax from Gillet et al, 2017
#     for tr in epicentral_stream:
#         if 'tmax' in tr.stats.keys():
#             print('tmax found ', tmax, tr.stats.title, tr.stats.station)
#             ax.hlines(y=tr.stats.tmax, xmin=tr.stats.distance_in_degree, xmax=tr.stats.distance_in_degree+taup_height, linewidth=2, color='r',zorder=9)

    if observations and trace_title_show:
        for i in range(len(epicentral_stream)):
            try:
                trace_title = epicentral_stream[i].stats.impact
            except AttributeError:
                trace_title = ''

            trace_title = '{} {}'.format(trace_title, epicentral_stream[i].stats.station)

            print('{} {:.01f}\n'.format(trace_title, epicentral_stream[i].stats.distance_in_degree))


            if seismograms_vertical:
                xy=(epicentral_stream[i].stats.distance_in_degree,annotate_time)
                rotation = -90
            else:
                xy=(annotate_time,epicentral_stream[i].stats.distance_in_degree)
                rotation = 0

            ax.annotate(text=trace_title, xy=xy,
                        xycoords='data', horizontalalignment='center',
                        verticalalignment='center',fontsize=10, color='k',
                        rotation=rotation,annotation_clip=False)

    if observations and show_impact_type:
        for i in range(len(epicentral_stream)):
            impact_type = ''
            try:
                trace_title = epicentral_stream[i].stats.impact
                if 'S-IVB' in trace_title:
                    impact_type = 'S-IVB'
                    impact_color = '#FF69B4'
                elif 'LM' in trace_title:
                    impact_type = 'LM'
                    impact_color = '#800080'
            except AttributeError:
                impact_type = ''


            plt.scatter(annotate_time2,epicentral_stream[i].stats.distance_in_degree,s=30,color=impact_color,label=impact_type,alpha=0.7,clip_on=False)

    filter_colors = ['#710193', '#9E7BB5', '#BE9394']
    if pre_filt_env is not None:

        # taper
        if observations:
            epicentral_stream_env.taper(max_percentage=None,max_length=taper_len,side='both')
        else:
            epicentral_stream_env.trim(starttime=UTCDateTime(-200.0),pad=True,fill_value=0.0)
            epicentral_stream_env.taper(max_percentage=None,max_length=taper_len,side='right')



        for ii, pre_filt in enumerate(pre_filt_env):

            # find a good smoothing size
            mid_period  = 1/((pre_filt[1] + pre_filt[2]) / 2)
            smooth_kernel_size = mid_period * smooth_periods / epicentral_stream_env[0].stats.delta
            smooth_kernel_size = int(smooth_kernel_size) + 1

            print('smooth_periods={}, smooth_kernel_size={} smoothing length={:.1f} s'.format(smooth_periods,smooth_kernel_size, smooth_kernel_size*tr.stats.delta))

            for i, tr1 in enumerate(epicentral_stream_env):

                tr = tr1.copy()
#                 tr.plot()
#                 print(pre_filt)
#                 return

                # process the traces for the envelopes
                if observations:

                    # add linear interpolation but keep the original mask
                    original_mask = linear_interpolation(tr,interpolation_limit=None)

                    if tr.stats.channel in ['MH1', 'MH2', 'MHZ','SHZ']:
                        # remove the instrument response
#                             pre_filt = [freqmin/1.5,freqmin,freqmax,freqmax*1.5]
                        tr.remove_response(inventory=inv, pre_filt=pre_filt, zero_mean=True, taper=True, output="DISP",
                                   water_level=None, plot=False)

                else:
                    # bandpass the simulated trace
                    tr.filter('bandpass', freqmin=pre_filt[1], freqmax=pre_filt[2],zerophase=False)

                # find the envelope
                tr.data=envelope(tr.data)

                if smooth_kernel_size > 1:
                    kernel = np.ones(smooth_kernel_size) / smooth_kernel_size
                    data_convolved = np.convolve(tr.data, kernel, mode='same')

                    tr.data = data_convolved

#                     # y_avg = np.zeros((1, len(tr)))
#                     # avg_mask = np.ones(smooth_length)/smooth_length
#                     # y_avg = np.convolve(tr.data, avg_mask, 'same')

#                     # tr.data = y_avg



                if observations:
                    # apply the mask back to the trace
                    tr.data = np.ma.masked_array(tr, mask=original_mask)

                if observations:
                    tr.trim(starttime=tr.stats.impact_time + startsecond, endtime=tr.stats.impact_time + endsecond)
                else:
                    tr.trim(starttime=UTCDateTime(startsecond), endtime=UTCDateTime(endsecond))

                if normalize == 'relative':
                    # find the normalization factor
                    max_data = abs(tr.data).max()
                    # normalize the trace
                    tr.data = tr.data/max_data
                    print('trace should be normalized')

                elif normalize == 'PP_arrival':
                    # this might be the best way to normalize
                    if 'PP_arrival' in tr.stats.keys():
            #                 print('Before ', tr.stats.starttime, tr.stats.endtime)
                        tr_local = tr.copy()
                        tr_local.trim(starttime=UTCDateTime(tr.stats.PP_arrival)-30.,endtime=UTCDateTime(tr.stats.PP_arrival)+70.)
            #                 tr.trim(starttime=tr.stats.S_arrival-40.)
            #                 print('After', tr.stats.starttime, tr.stats.endtime, tr.stats.endtime-tr.stats.starttime )
                        # find the normalization factor based on -30 and +70 seconds after the P arrival
                        if len(tr_local.data) > 0:
                            # calculate the normalization factor on the short section, but apply to longer trace
                            max_data = abs(tr_local.data).max()
                            tr.data = tr.data/max_data
                        else:
                            print('PP_arrival trace too small - removing the trace', tr)
                            epicentral_stream.remove(tr)
                    else:
                        print('PP_arrival not found - removing the trace', tr)
                        epicentral_stream.remove(tr)
                elif normalize == None:
                    pass
                else:
                    print('This type of normalize not implemented yet', normalize)
                    return

                if observations:
                    ref_time = tr.stats.impact_time
                else:
                    ref_time = UTCDateTime(0)

                # plot the envelope


                time = tr.times(reftime=ref_time)

#               # old way (was originally not squared)
#                 data = tr.data*scale + tr.stats.distance_in_degree
#                 data1 = tr.data*scale*-1.0 + tr.stats.distance_in_degree
#                 fill in area between the two lines
#                 plt.fill_betweenx(y=time,x1=data,x2=data1,label='{} - {} Hz'.format(pre_filt[1],pre_filt[2]), color=filter_colors[ii],alpha=0.5, zorder=9)

                if scale_list is not None:
                    data = tr.data*scale_list[i]*scale + tr.stats.distance_in_degree
                else:
                    data = tr.data*scale + tr.stats.distance_in_degree

                if fill_envelope:
                    data1 = tr.stats.distance_in_degree
                    # fill in area between the two lines
                    if seismograms_vertical:
                        plt.fill_betweenx(y=time,x1=data,x2=data1,label='Envelope ({}-{} Hz)'.format(pre_filt[1],pre_filt[2]), color=filter_colors[ii],alpha=0.5, zorder=9)
                    else:
                        plt.fill_between(x=time,y1=data,y2=data1,label='Envelope ({}-{} Hz)'.format(pre_filt[1],pre_filt[2]), color=filter_colors[ii],alpha=0.5, zorder=9)
                else:
                    if seismograms_vertical:
                        plt.plot(data,time,label='Envelope ({}-{} Hz)'.format(pre_filt[1],pre_filt[2]),color=filter_colors[ii],alpha=0.7,zorder=9)
                    else:
                        plt.plot(time,data,label='Envelope ({}-{} Hz)'.format(pre_filt[1],pre_filt[2]),color=filter_colors[ii],alpha=0.7,zorder=9)

    if catalogs is not None:
        station_latitude, station_longitude, station_elevation = get_station_details(inv)
        for catalog in catalogs:
            cat = read_events(catalog)
            for event in cat.events:
                if event.event_type == 'crash':

                    for st1 in ['S12', 'S14','S15','S16']:
                        begin = None
                        end = None
                        for pick in event.picks:
                            station_code = pick.waveform_id.station_code
                            if st1 != station_code:
                                continue

                            pick_public_id = pick.resource_id.id
    #                         if pick.phase_hint == 'P' and 'instrument' in pick_public_id:
            #                 and pick.waveform_id.station_code==station_code:

                            pick_mark = pick.time - event.origins[0].time
                            if pick_mark == 0.0:
                                print('valid pick not found')
                                continue

                            event_latitude = event.origins[0].latitude
                            event_longitude = event.origins[0].longitude



                            distance, azimuth_A_B, azimuth_B_A =  gps2dist_azimuth(
                              event_latitude, event_longitude, station_latitude[station_code], station_longitude[station_code], MOON_RADIUS_IN_KM*1000., MOON_FLATTENING)
                            distance_in_km = distance/1000.
                            distance_in_degree = kilometers2degrees(distance_in_km,radius=MOON_RADIUS_IN_KM)

                            zorder=8
                            if 'begin' in pick_public_id:
                                color = 'r'
                                label1 = 'Nunn - P - begin'
                                zorder=9
                                begin = pick_mark
                            elif 'end' in pick_public_id:
                                color = '#FFA500'
                                label1 = 'Nunn - P - end'
                                zorder=10
                                end = pick_mark
                            elif 'Logn' in catalog:
                                color = phase_list_colors[pick.phase_hint]
                                label1 = 'Logn. - {}'.format(pick.phase_hint)
                            else:
                                color='k'
                                label1 = 'Nunn - P - unknown'


                            print('pick mark ', pick_mark)
                            # if seismograms_vertical:
                            #     plt.hlines(y=pick_mark,xmin=distance_in_degree-taup_height, xmax=distance_in_degree+taup_height, color=color, label=label1, zorder=zorder)
                            # else:
                            #     plt.vlines(x=pick_mark,ymin=distance_in_degree-taup_height, ymax=distance_in_degree+taup_height, color=color, label=label1, zorder=zorder)
                    #
        #                     print(event.event_descriptions[0].text,distance_in_degree,'deg',pick_mark, 's', station_code)
        #     #                     plt.gca().axvline(x=pick_markP,
        # #                                           color=taup_color, linewidth=2)

                        if (begin is not None) and (end is not None):
    #                         plt.fill_betweenx(x=distance_in_degree,y1=,y2=data1,color='r',alpha=0.5, zorder=10
                            if seismograms_vertical:
                                plt.gca().add_patch( Rectangle(xy=(distance_in_degree-(taup_height), begin),
                                    width=taup_height*2, height=(end-begin),
                                    fc ='r', alpha=0.5, linewidth=None,label='P onset (est.)',zorder=1))
                            else:
                                # in case it's really small and we can't see it, add a line first
                                plt.vlines(begin, distance_in_degree-(taup_height), distance_in_degree+(taup_height),color='#EE8683',linewidth=2,zorder=1)
                                plt.gca().add_patch( Rectangle(xy=(begin,distance_in_degree-(taup_height)),
                                    width=(end-begin), height=taup_height*2,
                                    fc ='#EE8683', linewidth=None, label='P onset (est.)'))

    if taup_show:

        phase_color_dict = get_phase_colors(phase_list)

        # correct TauModel for source depth
        depth_corrected_model = model_taup.model.depth_correct(source_depth)

        # phase_names = parse_phase_list(phase_list)
        for i, phase in enumerate(phase_list):

            if 'kmps' in phase:
                arrivals = model_taup.get_travel_times(source_depth_in_km=source_depth, distance_in_degree=degree_max, phase_list=[phase])
                x = [0, degree_max ]
                y = [0, arrivals[0].time]
                phase_name = phase.replace('kmps', ' km/s')
                if seismograms_vertical:
                    ax.plot(x, y, label=phase_name, color=phase_color_dict[phase], linestyle='dashed', zorder=1)
                else:
                    ax.plot(y, x, label=phase_name, color=phase_color_dict[phase], linestyle='dashed', zorder=1)
            else:

                ph = SeismicPhase(phase, depth_corrected_model)
                # don't join lines across shadow zones
                for s in ph._shadow_zone_splits():




                    dist_deg = (180.0/np.pi)*ph.dist[s]
    #                 time_min = ph.time[s]/60
                    time_sec = ph.time[s]

                    if len(dist_deg) > 0:
                        phase_name = phase_name_dict.get(phase, phase)
                        linestyle = 'dashed'
                        zorder = 1
                        if just_taup and phase_name in ["P","S"]:
                            linestyle = 'solid'
                            zorder = 0.5

                        if seismograms_vertical:
                            ax.plot(dist_deg, time_sec, label=phase_name, color=phase_color_dict[phase],linestyle=linestyle, zorder=zorder)
                        else:
                            ax.plot(time_sec, dist_deg, label=phase_name, color=phase_color_dict[phase],linestyle=linestyle, zorder=zorder)

    x = [0,20]
    y = [0,280]
#     ax.plot(x, y, label='Middle of peaks', color='r',linestyle='solid')


    if raw:
        annotation = 'Unfiltered trace'
    else:
        annotation = 'Bandpass filtered between {:.2f} and {:.2f} Hz'.format(freqmin,freqmax)


    if normalize == 'relative':
        annotation = 'Normalized' + '\n' + annotation
    elif normalize == 'P_arrival':
        annotation = 'Normalized to P-arrival amplitude' + '\n' + annotation
    elif normalize == 'S_arrival':
        annotation = 'Normalized to S-arrival amplitude' + '\n' + annotation
    elif normalize == 'PP_arrival':
        annotation = 'Normalized to amplitude of first arrivals' + '\n' + annotation
    elif normalize in ['first_arrivals','first_arrivals_longer']:
        annotation = 'Normalized to first arrivals' + '\n' + annotation
    else:
        annotation = 'Traces not normalized' + '\n' + annotation

    # if model_taup is not None:
    #     annotation = taup_label + '\n' + annotation

    print(annotation)

    if seismograms_vertical:
        # ax.annotate(text=annotation, xy=(0.99,0.01), xycoords='axes fraction',
        #             horizontalalignment='right', verticalalignment='bottom',fontsize=10, color='k')
        plt.xlabel('Epicentral Distance [degrees]')
        plt.ylabel('Time [s]')
    else:
        # ax.annotate(text=annotation, xy=(0.99,0.01), xycoords='axes fraction',
        #             horizontalalignment='right', verticalalignment='bottom',fontsize=10, color='k')
        plt.xlabel('Time [s]')
        plt.ylabel('Epicentral Distance [degrees]')



    # plot legend, avoiding duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if show_legend:
        ax_legend = ax.legend(by_label.values(), by_label.keys(), loc=legend_loc, numpoints=1, ncol=ncol,columnspacing=0.8,handletextpad=0.1)
        ax_legend.set_zorder(100)
    if seismograms_vertical:
        plt.xlim(degree_min,degree_max)
        plt.ylim(startsecond,endsecond)
    else:
        plt.xlim(startsecond,endsecond)
        plt.ylim(degree_min,degree_max)

#     print(freq)
    if title is not None:
        plt.title(title)
#     plt.title('Propagation of an Explosive Source\nLowpass filtered below {:.1f} Hz'.format(freq))

    if save_fig:
        plt.tight_layout()
        fig_name = 'fig_{}.png'.format(UTCDateTime.now())
        fig_name = os.path.join('temp/',fig_name)
        plt.savefig(fig_name)
        print(fig_name)

    plt.show()
    plt.close()

    return epicentral_stream

### Other Code
channel_dict = {
    'U1' : 'R',
    'U2' : 'T',
    'U3' : 'Z',
}

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

# get_station_details()
# print(station_latitude)
# print(station_longitude)
# print(station_elevation)

# def get_phase_colors(phase_list, name='muted'):
#
#     phase_list_extra = []
#     for phase in phase_list:
#         if phase not in ['P','S','pS']:
#             phase_list_extra.append(phase)
#
#     palette = seaborn.color_palette(name, len(phase_list_extra)+3).as_hex()
#     phase_list_dict =  {
#     'P': palette[0],
#     'S': palette[1],
#
#     'pS': palette[2],
#     }
#
#     for i, phase in enumerate(phase_list_extra):
#         if phase not in ['P','S','pS']:
#             phase_list_dict[phase] = palette[i+3]
#
#     return (phase_list_dict)

def get_phase_colors(phase_list):

    # give some basic phases the same colors everytime

    phase_list_extra = []
    for phase in phase_list:
        if phase not in ['P','p','S', 's', 'PS']:
            phase_list_extra.append(phase)

    phase_list_dict =  {
    'P': '#4363d8', # blue
    'p': '#42d4f4', # cyan
    'S': '#911eb4', # purple
    's': '#dcbeff', # lavender
    'PS':'#f032e6', # magenta

    }

    remaining_list = [
        '#000075', # navy
        '#e6194B', # red
        '#f58231', # orange
        '#ffe119', # yellow
        '#bfef45', # lime
        '#3cb44b', # green
        '#469990', # teal
        '#800000', # maroon
        '#9A6324', # brown
        '#808000', # olive
        '#fabed4', # pink
        '#aaffc3', # mint
    ]

    for i, phase in enumerate(phase_list_extra):
        if phase not in ['P','p','S', 's', 'PS']:
            phase_list_dict[phase] = remaining_list[i]

    return (phase_list_dict)

# plot taup and the simulations or observations - the plot shows the seismograms horizontally
def plot_envelope_taup(original_stream=None,original_stream_dict=None,run_list=[],
    observation_stream=None,
    distance_in_degree=60,inv=None,
    model_taup=None,model_taup_label=None,taup_show=True,title=None,
    freqmin=None,freqmax=None,channel='Z',obs_start=-1800,startsecond=0,endsecond=1800,normalize='relative',
    scale_list=None,
    taper_len=10,phase_list=["P", "PP", "PcP", "Pdiff", "PvmP"],source_depth=0.001,
    smooth_periods=10,plot_seismogram=False,plot_envelope=False,plot_envelope_one_color=False,plot_derivative=False,annotate_relative=False,save_fig=False,figsize=(11, 17)):

    annotate_endsecond = endsecond*.99
    annotate_endsecond2 = endsecond*1.1

    if model_taup_label is not None:
        print('Calculated using model : ', model_taup_label)
        taup_label1 = model_taup_label.replace('_taup','')
        taup_label = "Seismic phases: {}".format(taup_label1)
    else:
        taup_label1 = ''
        taup_label = ''

    # cmap_phaselist = get_cmap_phaselist(len(phase_list))
    # seaborn_colors = seaborn.color_palette()
    seaborn_colors = seaborn.color_palette(n_colors=14)



    color_match = {
    '120_VPREMOON_atten_explosion_2':'#1E4538',
    '128b_VPREMOON_atten_surface_2': '#BCBD45',
    '125_VPREMOON_atten_Moho_2': 'g',
    '127_VPREMOON_atten_linear20_2': '#73A848',
    '126_VPREMOON_atten_linear50_2': '#2E6C4B',
    '124_VPREMOON_atten_linear80_2': '#3C867F',
    '152_VPREMOON_atten_combi_50_2': '#8CE1B8',

    '141_ISSI_atten_explosion_2': '#05007B',
    '148_ISSI_atten_surface_2': '#2C2B4F',
    '147_ISSI_atten_Moho_2': '#52B2F9',
    '146_ISSI_atten_linear20_2': '#536F86',
    '145_ISSI_atten_linear50_2': '#377DF6',
    '140_ISSI_atten_linear80_2': '#74FAFD',
    '150_ISSI_atten_combi_50_2': '#60B4CA',
    }



    # seaborn_colors = seaborn.palplot(seaborn.color_palette('Spectral',20))

#     # find a good smoothing size
#     mid_period  = 1/((freqmin + freqmax) / 2)
#     smooth_kernel_size = mid_period * smooth_periods / epicentral_stream_env[0].stats.delta
#     smooth_kernel_size = int(smooth_kernel_size) + 1

#     print('smooth_periods={}, smooth_kernel_size={} smoothing length={:.1f} s'.format(smooth_periods,smooth_kernel_size, smooth_kernel_size*tr.stats.delta))


    # start a figure
    fig, ax = plt.subplots(figsize=figsize,dpi=300)
    # twin the x axis
    ax2 = ax.twinx()

    # default limit, if there's no observations
    y_upper_lim = -1

# select the streams to work on
# working on the epicentral stream, remove the instrument response, remove the mean if necessary
# and taper if necessary
# plot the epicentral stream
    # First process the observations, if present
    observation_stream_local = Stream()
    if observation_stream is not None:

        if channel == 'Z':
            channel1 = 'MHZ'
        else:
            print('needs more work to rotate')
            return
        for tr in observation_stream.select(channel=channel1):
            if math.isclose(tr.stats.distance_in_degree,distance_in_degree,abs_tol=0.01):
                tr1 = tr.copy()
                observation_stream_local.append(tr1)

        # make some basic preparations for the traces (remove the -1 values, remove the single sample gaps, taper, remove the mean)
        remove_negative_ones(observation_stream_local)
        for tr in observation_stream_local:
            # interpolate across the gaps of one sample
            linear_interpolation(tr,interpolation_limit=1)

        observation_stream_local.taper(max_percentage=None,max_length=taper_len,side='both')

        # find a good smoothing size
        mid_period  = 1/((freqmin + freqmax) / 2)
        smooth_kernel_size = mid_period * smooth_periods / observation_stream_local[0].stats.delta
        smooth_kernel_size = int(smooth_kernel_size) + 1

        print('Observations: smooth_periods={}, smooth_kernel_size={} smoothing length={:.1f} s'.format(smooth_periods,smooth_kernel_size, smooth_kernel_size*observation_stream_local[0].stats.delta))

        for i, tr in enumerate(observation_stream_local):
            if scale_list is not None:
                scale = scale_list[0]
            else:
                scale = 1

            y_upper_lim = scale
            y_lower_lim = -1*scale
            # add linear interpolation but keep the original mask
            original_mask = linear_interpolation(tr,interpolation_limit=None)

            if tr.stats.channel in ['MH1', 'MH2', 'MHZ','SHZ']:
                # remove the instrument response
                pre_filt = [freqmin/1.5,freqmin,freqmax,freqmax*1.5]
                tr.remove_response(inventory=inv, pre_filt=pre_filt, zero_mean=True, taper=True, output="DISP",
                           water_level=None, plot=False)

                if plot_envelope or plot_derivative or plot_envelope_one_color:
                    tr_envelope = tr.copy()
                    # find the envelope
                    tr_envelope.data=envelope(tr_envelope.data)

                    # apply the mask back to the trace
                    tr_envelope.data = np.ma.masked_array(tr_envelope, mask=original_mask)
                    tr_envelope.trim(starttime=tr.stats.impact_time + startsecond, endtime=tr_envelope.stats.impact_time + endsecond)

                    # normalize BEFORE smoothing
                    if normalize == 'relative':
                        # find the normalization factor
                        max_data = abs(tr_envelope.data).max()
                        # normalize the trace
                        tr_envelope.data = tr_envelope.data/max_data
    #                     print('trace should be normalized')

                    if smooth_periods > 1:
                        kernel = np.ones(smooth_kernel_size) / smooth_kernel_size
                        data_convolved = np.convolve(tr_envelope.data, kernel, mode='same')

                        tr_envelope.data = data_convolved

#                     # y_avg = np.zeros((1, len(tr)))
#                     # avg_mask = np.ones(smooth_length)/smooth_length
#                     # y_avg = np.convolve(tr.data, avg_mask, 'same')

#                     # tr.data = y_avg





                if plot_seismogram:
                    # apply the mask back to the trace
                    tr.data = np.ma.masked_array(tr, mask=original_mask)

                    tr.trim(starttime=tr.stats.impact_time + startsecond, endtime=tr.stats.impact_time + endsecond)

                    if normalize == 'relative':
                        # find the normalization factor
                        max_data = abs(tr.data).max()
                        # normalize the trace
                        tr.data = tr.data/max_data
    # #                     print('trace should be normalized')
    #                     if annotate_relative:
    #                         # scale_factor = scale/max_data
    #                         ax.annotate(text='x{:.1e}'.format(max_data), xy=(annotate_endsecond2,-1), xycoords='data',
    #                                     horizontalalignment='left', verticalalignment='bottom',fontsize=10, color='k',annotation_clip=False,
    #                                     )
    #

                if plot_envelope:
                    ref_time = tr.stats.impact_time
                    time = tr_envelope.times(reftime=ref_time)

                    # plot the envelope
                    ax.plot(time,tr_envelope.data*scale,label='Observation - Envelope',color='k',alpha=0.5,zorder=8)

                if plot_envelope_one_color:
                    ref_time = tr.stats.impact_time
                    time = tr_envelope.times(reftime=ref_time)

                    # plot the envelope
                    ax.plot(time,tr_envelope.data*scale,color='#710193',alpha=0.5,zorder=9)

                if plot_derivative:

                    tr_diff = tr.copy()
                    tr_diff.differentiate()
                    if smooth_periods > 1:
                        kernel = np.ones(smooth_kernel_size) / smooth_kernel_size
                        data_convolved = np.convolve(tr_diff.data, kernel, mode='same')
                        tr_diff.data = data_convolved
                    ref_time = tr.stats.impact_time
                    time = tr.times(reftime=ref_time)
                    ax2.plot(time,tr_diff.data*scale*100,label='Observation - Derivative',color='r',alpha=0.5,zorder=9)

                if plot_seismogram:
                    # plot the envelope
                    ref_time = tr.stats.impact_time
                    time = tr.times(reftime=ref_time)
                    ax.plot(time,tr.data*scale,color='k',alpha=0.5,zorder=9)
                    ax.annotate(text=tr.stats.title, xy=(annotate_endsecond,0.5), xycoords='data',
                                horizontalalignment='right', verticalalignment='bottom',fontsize=10, color='k',annotation_clip=False,
                                path_effects=[patheffects.withStroke(linewidth=4, foreground="white")]
                                )

    ################################################
    # Plot the simulations
    epicentral_stream = Stream()

    if original_stream is not None:
        for tr in original_stream.select(channel=channel):
            if math.isclose(tr.stats.distance_in_degree,distance_in_degree,abs_tol=0.01):
                tr1 = tr.copy()
                epicentral_stream.append(tr1)

    else:
        for run in run_list:
            for tr in original_stream_dict[run].select(channel=channel):
                if math.isclose(tr.stats.distance_in_degree,distance_in_degree,abs_tol=0.01):
                    tr1 = tr.copy()
                    epicentral_stream.append(tr1)


#     if pre_filt_env is not None:
#         epicentral_stream_env = epicentral_stream.copy()



#     annotate_time = (endsecond - 1.07*(endsecond-startsecond))

#     fig = plt.figure(dpi=150)


#     nasa_red = '#FC3D21'
#     nasa_blue = '#1B3C8C'
#     bright_blue = '#1B3CF8'
#     blue_grey =  '#454756'


    epicentral_stream.trim(starttime=UTCDateTime(-200.0),pad=True,fill_value=0.0)
    epicentral_stream.taper(max_percentage=None,max_length=taper_len,side='right')

    if observation_stream is None:
        scale_len = 0
    else:
        scale_len = 1

    for i, tr in enumerate(epicentral_stream):
        if scale_list is not None:
            scale = scale_list[i+scale_len]
        else:
            scale = 1

        color_i = run_list.index(tr.stats.run)

        color_match_col = color_match.get(tr.stats.run, 'r')

        # process the traces for the envelopes

        # bandpass the simulated trace
        tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=False)

        if plot_envelope or plot_derivative or plot_envelope_one_color:
            tr_envelope = tr.copy()
            # find the envelope
            tr_envelope.data=envelope(tr_envelope.data)

        if plot_seismogram:
            tr.trim(starttime=UTCDateTime(startsecond), endtime=UTCDateTime(endsecond))

            if normalize == 'relative':
                # find the normalization factor
                max_data = abs(tr.data).max()
                # normalize the trace
                tr.data = tr.data/max_data

                max_data_norm = max_data / scale
                print(i, max_data_norm)

                if i == 0:

                    max_data_norm0 = max_data_norm

                # the first in the list has a amplification factor of 1 - the
                # rest are relative to it - it relates to how they are displayed
                # so if it's strongly magnified (for example to better see the
                # coda), the multiplificaiton factor would be greater
                ampl_factor = max_data_norm0/max_data_norm

                if annotate_relative:
                    # scale_factor = scale/max_data
                    ax.annotate(text='x{:.1f}'.format(ampl_factor), xy=(annotate_endsecond2,-(i+1)*2+0.5), xycoords='data',
                                horizontalalignment='right', verticalalignment='bottom',fontsize=10, color=color_match_col,annotation_clip=False,
                                )

    #                     print('trace should be normalized')

        if plot_envelope or plot_derivative or plot_envelope_one_color:

            # normalize BEFORE smooothing
            if normalize == 'relative':
                # find the normalization factor
                max_data = abs(tr_envelope.data).max()
                # normalize the trace
                tr_envelope.data = tr_envelope.data/max_data
    #                     print('trace should be normalized')
                print('plot envelope ', abs(tr_envelope.data).max())

            tr_envelope.trim(starttime=UTCDateTime(startsecond), endtime=UTCDateTime(endsecond))

            # find a good smoothing size
            mid_period  = 1/((freqmin + freqmax) / 2)
            smooth_kernel_size = mid_period * smooth_periods / tr_envelope.stats.delta
            smooth_kernel_size = int(smooth_kernel_size) + 1

            print('Simulation: {}, smooth_periods={}, smooth_kernel_size={} smoothing length={:.1f} s'.format(tr_envelope.stats.run,smooth_periods,smooth_kernel_size, smooth_kernel_size*epicentral_stream[0].stats.delta))

            if smooth_periods > 1:
                kernel = np.ones(smooth_kernel_size) / smooth_kernel_size
                data_convolved = np.convolve(tr_envelope.data, kernel, mode='same')
                tr_envelope.data = data_convolved

        ref_time = UTCDateTime(0)

        if plot_seismogram:
            # plot the seismogram
            time = tr.times(reftime=ref_time)
            # XXXX
            ax.plot(time,tr.data*scale-(i+1)*2,color=color_match_col,alpha=0.5,zorder=7)
            ax.annotate(text=tr.stats.short_title, xy=(annotate_endsecond,-(i+1)*2+0.5), xycoords='data',
                        horizontalalignment='right', verticalalignment='bottom',fontsize=10, color=color_match_col,annotation_clip=False,
                        )


        if plot_envelope_one_color:
            # plot the envelope
            time = tr_envelope.times(reftime=ref_time)
            ax.plot(time,tr_envelope.data*scale-(i+1)*2,color='#710193',alpha=0.5,zorder=9)

        if plot_envelope:
            # plot the envelope
            time = tr_envelope.times(reftime=ref_time)
            ax.plot(time,tr_envelope.data*scale-(i+1)*2,label=tr_envelope.stats.run,color=color_match_col,alpha=0.5,zorder=7)

        if plot_derivative:
            tr_diff = tr_envelope.copy()
#             for ix in range(0,len(tr_diff.data)):
#                 tr_diff.data[ix] = 2*ix**2
#             for ix in range(0,200):
#                 tr_diff.data[ix] = tr_diff.data[199]
            tr_diff.differentiate()
            if smooth_periods > 1:
                kernel = np.ones(smooth_kernel_size) / smooth_kernel_size
                kernel = np.ones(smooth_kernel_size) / smooth_kernel_size
                data_convolved = np.convolve(tr_diff.data, kernel, mode='same')
                tr_diff.data = data_convolved

            time = tr_diff.times(reftime=ref_time)
            ax.plot(time,tr_diff.data*scale*100-(i+1)*2,label='{} Derivative'.format(tr_diff.stats.run),color='g',alpha=0.5,zorder=9)

        y_lower_lim = -(i+1)*2-scale

    if taup_show:
        phase_list_dict = get_phase_colors(phase_list)
        print('Source Depth={} km, Distance in degrees {}'.format(source_depth, distance_in_degree))
#         print(type(model_taup))
#         depth_corrected_model = model_taup.model.depth_correct(source_depth)
#         print(type(depth_corrected_model))
# #         phase_names = sorted(parse_phase_list(phase_list))
#         print(source_depth, distance_in_degree)
        arrivals = model_taup.get_travel_times(source_depth_in_km=source_depth,
                                  distance_in_degree=distance_in_degree,phase_list=phase_list)
        phases_found = []
        for arrival in arrivals:
            if arrival.phase.name not in phases_found:
                print('Suppressing secondary arrivals')
                ax.axvline(arrival.time, label=arrival.phase.name, color=phase_list_dict[arrival.phase.name], zorder=2, linestyle='dashed')
                phases_found.append(arrival.phase.name)

    # Here!!!

    if normalize == 'relative':
        annotation = 'Normalized'
    else:
        annotation = 'Traces not normalized'
    annotation = '{}\n{:.2f}-{:.2f} Hz'.format(annotation,freqmin,freqmax)


    # if model_taup is not None:
    #     annotation = taup_label + '\n' + annotation

    ax.annotate(text=annotation, xy=(1,-0.06), xycoords='axes fraction',
                horizontalalignment='right', verticalalignment='top',fontsize=10, color='k')

    if plot_derivative:
        plt.ylabel('Smoothed Derivative*100', color='g')


    # plot legend, avoiding duplicate labels
    # handles, labels = ax.get_legend_handles_labels()
    # by_label = dict(zip(labels, handles))


    # ax_legend = ax.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, numpoints=1, framealpha=0.5)
    # ax_legend.set_zorder(100)

    ax_legend = ax.legend(loc='upper left', bbox_to_anchor=(0, -0.03), numpoints=1, ncol=2, frameon=False, fontsize=10)
    ax_legend.set_title(model_taup_label,prop={'size':'small'})

    ax.get_yaxis().set_ticks([])
    ax2.get_yaxis().set_ticks([])

    if title is not None:
        plt.title(title)
    else:
        plt.title('Observations and Simulations: {}$\degree$'.format(round(distance_in_degree,2)))
    plt.xlim(startsecond,endsecond)

    ax.set_ylim(y_lower_lim,y_upper_lim)

    ax.set_xlabel('Time [s]')

    if save_fig:
        plt.tight_layout()
        print('now using tight_layout - check!!!')
        fig_name = 'fig_{}.png'.format(UTCDateTime.now())
        fig_name = os.path.join('./temp/',fig_name)
        plt.savefig(fig_name)
        print(fig_name)

    plt.show()
    plt.close()

    # XXXX

    return epicentral_stream


# plot both taup models
def plot_taup_epicentral(
    calculated_taup_dict=None,taup_show=True,title=None,
    freqmin=None,freqmax=None,channel='Z',obs_start=-1800,startsecond=0,endsecond=1800,normalize='relative',scale=20,taup_height=100,
    observations=True,taper_len=10,phase_list=["P", "PP", "PcP", "Pdiff", "PvmP"],source_depth=0.01,degree_min=0,degree_max=180):


    # Pretty paired colors. Reorder to have saturated colors first and remove
    # some colors at the end.
    cmap = get_cmap('Paired', lut=12)
    COLORS = ['#%02x%02x%02x' % tuple(int(col * 255) for col in cmap(i)[:3])
              for i in range(12)]
    COLORS = COLORS[1:][::2][:-1] + COLORS[::2][:-1]

    for i in range(len(phase_list)):
        c = COLORS[i % len(COLORS)]
        # print(c)

    prefix = ['VPREMOON', 'ISSI']
    styles = ['solid', 'dashed']

    fig, ax = plt.subplots(figsize=(15, 15))
    for ii, m in enumerate(['VPREMOON_taup','ISSI_MOON_M1_taup']):
        print('Calculated using model : ', m)
        taup_calculated_model = calculated_taup_dict[m]

        style = styles[ii]



    #     fig = plt.figure(dpi=150)


        if taup_show:
            depth_corrected_model = taup_calculated_model.model.depth_correct(source_depth)

            phase_names = sorted(parse_phase_list(phase_list))
            for i, phase in enumerate(phase_names):
                ph = SeismicPhase(phase, depth_corrected_model)
                # don't join lines across shadow zones
                for s in ph._shadow_zone_splits():
                    dist_deg = (180.0/np.pi)*ph.dist[s]
    #                 time_min = ph.time[s]/60
                    time_sec = ph.time[s]
                    c = COLORS[i % len(COLORS)]
                    if len(dist_deg) > 0:
                        plot_all = True
                        if plot_all:
                            # wrap-around plotting
                            while dist_deg[0] > 360.0:
                                dist_deg = dist_deg - 360.0
                            ax.plot(dist_deg, time_sec, label='{} {}'.format(prefix[ii],phase), color=c, linestyle=style)
                            ax.plot(dist_deg - 360.0, time_sec,
                                    label=None, color=c)
                            ax.plot(360.0 - dist_deg, time_sec,
                                    label=None, color=c)
                        else:
                            ax.plot(dist_deg, time_sec, label='{} {}'.format(prefix[ii],phase), color=c, linestyle=style)


    plt.xlabel('Epicentral Distance [degrees]')
    plt.ylabel('Time [s]')
    # plot legend, avoiding duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right', numpoints=1)
    plt.xlim(degree_min,degree_max)
    plt.ylim(startsecond,endsecond)
#     print(freq)
    if title is not None:
        plt.title(title)
#     plt.title('Propagation of an Explosive Source\nLowpass filtered below {:.1f} Hz'.format(freq))
    plt.show()
    plt.close()
