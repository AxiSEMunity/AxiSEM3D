#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Util code for the jupyter notebooks - observations 
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

import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
from matplotlib.dates import date2num 
from matplotlib.patches import Rectangle

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
from obspy.clients.fdsn.client import Client
from obspy.clients.fdsn.header import FDSNNoDataException

from irfpy.moon import moon_map

from pdart.util import linear_interpolation, remove_negative_ones, remove_negative_ones_trace

import seaborn


MOON_RADIUS_IN_KM = 1737.1
# MOON_FLATTENING = 1/8255
MOON_FLATTENING = 0.0 # flattening for the spherical model

def get_station_details(inv):

    station_latitude = {}
    station_longitude = {}
    station_elevation = {}
    
    for network in inv:
        for station in network: 
            station_latitude[station.code] = station.latitude
            station_longitude[station.code] = station.longitude
            station_elevation[station.code] = station.elevation 

    return (station_latitude, station_longitude, station_elevation)


## Section 3 - Code to get a catalog and connect to IRIS to get a stream from the catalog
def get_observations_any(channel='MHZ',model_taup=None,catalog=None,inv=None,):

    if inv is None: 
        client = Client("IRIS")
        # get the response file (wildcards allowed)
        inv = client.get_stations(starttime=UTCDateTime('1969-01-01'),endtime=UTCDateTime('1977-09-30T23:59:59'),
            network='XA', sta='*', loc='*', channel='*',
            level="response")

    station_latitude, station_longitude, station_elevation = get_station_details(inv)

    cat = read_events(catalog)    
    
    # print(inv)
    degrees = []
    original_observation_stream = Stream()
    for event in cat:

        if len(event.origins) > 1:
            print('More than one origin for this event. Using the first one.')

        origin = event.origins[0]

        event_latitude = origin.latitude
        event_longitude = origin.longitude
        event_depth = origin.depth

        origin_time = UTCDateTime(origin.time)
        # here 

        print(origin_time, event_latitude, event_longitude, event_depth)

        for network in inv:
            for station in network: 
                if ((origin_time >= station.start_date) and (origin_time <= station.end_date)):
    #                 print(station)
    #                 print(impact, event_latitude, event_longitude, event_elevation)
    #                 print(station.code, station.latitude, station.longitude, station.elevation, station.start_date, station.end_date)
                    
                    distance, azimuth_A_B, azimuth_B_A =  gps2dist_azimuth(
                          event_latitude, event_longitude, station.latitude, station.longitude, MOON_RADIUS_IN_KM*1000., MOON_FLATTENING)
                    distance_in_km = distance/1000.
                    distance_in_degree = kilometers2degrees(distance_in_km,radius=MOON_RADIUS_IN_KM)
                
                    
                    # if event_depth is not None: 
                    #     elevation_diff = station.elevation - event_depth
                    # else:
                    #     elevation_diff = ''
                    print(station.code, distance_in_km, distance_in_degree, event_depth, station.elevation)
                    degrees.append(distance_in_degree)

                    st = None
#                     if distance_in_km < 1000:
#                     TODO
                    try:  
                        client = Client("IRIS")
                        st = client.get_waveforms('XA', station.code, '*', channel, origin_time-1800, origin_time+5400) 
                    except FDSNNoDataException: 
                        print('Data not found for ', station.code, str(origin_time))
                        
                    if st is not None: 
                        for tr in st:
                            tr.stats.distance_in_km = distance_in_km
                            tr.stats.distance_in_degree = distance_in_degree
                            tr.stats.impact_time = origin_time
                            tr.stats.title = str(origin_time) + '\n' + station.code.strip()
                            tr.stats.run = 'Observations'
                            tr.stats.depth = event_depth
                        
                        if channel == 'MH*':
                            stream_MHZ = st.select(channel='MHZ')
                            stream_MH1 = st.select(channel='MH1')
                            if len(stream_MHZ) > 0:
                                original_observation_stream.extend(stream_MHZ)
                            else:
                                original_observation_stream.extend(stream_MH1)

                # temporarily just getting the first one
                # break 

        # temporarily just getting the first one
        # break 


#     degrees.sort()
#     # print(degrees)
#     for d in degrees:
#         print("{:.4f}".format(d))
        
    original_observation_stream.sort(keys=['distance_in_degree'])
        
    return original_observation_stream

## Section 3 - Code to get list of Impact Parameters and connect to IRIS to get a stream of the artificial impacts
def get_observations(channel='MHZ'):

    df = pd.read_csv('./input_files/ImpactParameters.csv')

    print(df.to_string())
#     return
    
    #  view the impact parameters 
    print(df[["Impact","Latitude","Longitude","Elevation (m)","UST (time recorded on Earth)"]].to_string())
    
    client = Client("IRIS")
    # get the response file (wildcards allowed)
    inv = client.get_stations(starttime=UTCDateTime('1969-01-01'),endtime=UTCDateTime('1977-09-30T23:59:59'),
        network='XA', sta='*', loc='*', channel='*',
        level="response")

    station_latitude, station_longitude, station_elevation = get_station_details(inv)    
    
    # print(inv)
    degrees = []
    original_observation_stream = Stream()
    for index, row in df.iterrows():

    #     print(row[["Impact","Latitude","Longitude","Elevation (m)","UST (time recorded on Earth)"]])
        impact = row["Impact"]
        print(impact)
    #     if impact == 'A11/LM':
    #         print(row)
    #     continue

        event_latitude1 = row["Latitude"]
        event_longitude1 = row["Longitude"]
        event_elevation1 = row["Elevation (m)"]
      
        try:
            event_latitude = float(event_latitude1)
        except ValueError:
            event_latitude = None
            
        try:
            event_longitude = float(event_longitude1)
        except ValueError:
            event_longitude = None
        
        try:
            event_elevation = float(event_elevation1)
        except ValueError:
            event_elevation = None


        try: 
            impact_time = UTCDateTime(row["UST (time recorded on Earth)"]) 
        except TypeError:
            print('Impact not found for {}'.format(impact))
            continue
            
        print(impact, impact_time, event_latitude, event_longitude, event_elevation)


        for network in inv:
            for station in network: 
                if ((impact_time >= station.start_date) and (impact_time <= station.end_date)
                    and 'ASCENT' not in impact.upper()):
    #                 print(station)
    #                 print(impact, event_latitude, event_longitude, event_elevation)
    #                 print(station.code, station.latitude, station.longitude, station.elevation, station.start_date, station.end_date)
                    
                    distance, azimuth_A_B, azimuth_B_A =  gps2dist_azimuth(
                          event_latitude, event_longitude, station.latitude, station.longitude, MOON_RADIUS_IN_KM*1000., MOON_FLATTENING)
                    distance_in_km = distance/1000.
                    distance_in_degree = kilometers2degrees(distance_in_km,radius=MOON_RADIUS_IN_KM)
                
                    
                    if event_elevation is not None: 
                        elevation_diff = station.elevation - event_elevation
                    else:
                        elevation_diff = ''
                    print(station.code, distance_in_km, distance_in_degree, elevation_diff, event_elevation, station.elevation)
                    degrees.append(distance_in_degree)

                    st = None
#                     if distance_in_km < 1000:
#                     TODO
                    try:  
                        st = client.get_waveforms('XA', station.code, '*', channel, impact_time-1800, impact_time+5400) 
                    except FDSNNoDataException: 
                        print('Data not found for ', station.code, impact)
                        
                    if st is not None: 
                        for tr in st:
                            tr.stats.distance_in_km = distance_in_km
                            tr.stats.distance_in_degree = distance_in_degree
                            tr.stats.impact_time = impact_time
                            tr.stats.title = impact.strip() + '\n' + station.code.strip()
                            tr.stats.run = 'Observations'
                        
                        if channel == 'MH*':
                            stream_MHZ = st.select(channel='MHZ')
                            stream_MH1 = st.select(channel='MH1')
                            if len(stream_MHZ) > 0:
                                original_observation_stream.extend(stream_MHZ)
                            else:
                                original_observation_stream.extend(stream_MH1)


#     degrees.sort()
#     # print(degrees)
#     for d in degrees:
#         print("{:.4f}".format(d))
        
    original_observation_stream.sort(keys=['distance_in_degree'])
        
    return original_observation_stream

# X  Station S11 (ALSEP 11, Mare Tranquillitatis, Moon)
# 	Station Code: S11
# 	Channel Count: 0/5 (Selected/Total)
# 	1969-07-21T03:00:00.000000Z - 1969-08-26T04:00:00.000000Z
# 	Access: open 
# 	Latitude: 0.67, Longitude: 23.47, Elevation: -1929.0 m
# 	Available Channels:   

# get the observations from the local drive and add the origin times and distances 

def get_observations_local(channel='MH*',model_taup=None,catalog=None,inv=None,
    local_data_dir='input_files/local_MSEED/'
):
    
    cat = read_events(catalog)

    overall_stream = read(os.path.join(local_data_dir,'*.MSEED'))



    station_latitude, station_longitude, station_elevation = get_station_details(inv)
    
    # include tmax from Gillet et al., 2017
    # df_model = pd.read_csv('./files/tmaxtd_impacts.csv',names=['Station','Event','Comp','Year','Month','Day','Hour','Min','Sec','Lat','Lon','Dep','EDist','Edisterr','Tmax','Tmaxerr','Tcoda','Taud','Tauderr','CC'],sep=None,skiprows=1,skipinitialspace=0,comment='#')
        
    observation_stream = Stream()
            
    for event in cat.events:
        if event.event_type == 'crash':
            impact_time = event.origins[0].time
            impact = event.event_descriptions[0].text

            for sta in ['S12','S14','S15','S16']:
            
                event_stream = Stream()
                for tr1 in overall_stream.select(station=sta,channel=channel):
                    if tr1.stats.starttime < impact_time and tr1.stats.endtime > impact_time: 
                        event_stream.append(tr1)
                        
                # check the special case of 'MH*' and keep only one (MHZ if available, MH1 if not)  
                if channel == 'MH*':
                    mhz_found = False
                    if len(event_stream.select(channel='MHZ')) == 1:
                        mhz_found = True
                    for tr in event_stream:
                        if mhz_found: 
                            if tr.stats.channel in ('MH1', 'MH2'):
                                event_stream.remove(tr)
                        else:
                            if tr.stats.channel in ('MH2'):
                                event_stream.remove(tr)
                                

                if len(event_stream) == 0:
                    continue
                    
                for tr in event_stream:
                    
                    event_latitude = event.origins[0].latitude 
                    event_longitude = event.origins[0].longitude 
                    event_depth = event.origins[0].depth 
                
                    distance, azimuth_A_B, azimuth_B_A =  gps2dist_azimuth(
                        event_latitude, event_longitude, station_latitude[sta], station_longitude[sta], MOON_RADIUS_IN_KM*1000., MOON_FLATTENING)

                    distance_in_km = distance/1000.

                    distance_in_degree = kilometers2degrees(distance_in_km,radius=MOON_RADIUS_IN_KM)

                    tr.stats.distance_in_km = distance_in_km
                    tr.stats.distance_in_degree = distance_in_degree
                    tr.stats.event_depth = event_depth
                    tr.stats.title = event.event_descriptions[0].text
                    print(impact, tr.stats.distance_in_degree)
                    
                    tr.stats.impact_time = impact_time
                    
                    tr.stats.impact = impact
                    
                    tr.stats.run = 'Observations (Local)'
            
                    if model_taup is not None: 
                        arrivals = model_taup.get_travel_times(source_depth_in_km=0, distance_in_degree=distance_in_degree, phase_list=(['P']))
                        if len(arrivals) > 0: 
                            tr.stats.P_arrival = impact_time + arrivals[0].time
                        arrivals = model_taup.get_travel_times(source_depth_in_km=0, distance_in_degree=distance_in_degree, phase_list=(['S']))
                        if len(arrivals) > 0: 
                            tr.stats.S_arrival = impact_time +  arrivals[0].time
                        arrivals = model_taup.get_travel_times(source_depth_in_km=0, distance_in_degree=distance_in_degree, phase_list=(['PP']))
                        if len(arrivals) > 0: 
                            tr.stats.PP_arrival = impact_time + arrivals[0].time
                            

                    # tmax from Gillet et al., 2017
                    # if tr.stats.channel == 'MHZ':
                    #     index1 = df_model.loc[(df_model.Event == impact) & (df_model.Station == sta) & (df_model.Comp == 'BRB1Z')].index.tolist()
                    #     if len(index1) > 0:
                    #         i = index1[0]
                    #         tr.stats.tmax = df_model.Tmax.iloc[i]
                    # 
                    # if tr.stats.channel in ['MH1','MH2']:
                    #     index1 = df_model.loc[(df_model.Event == impact) & (df_model.Station == sta) & (df_model.Comp == 'BRB1XY')].index.tolist()
                    #     if len(index1) > 0:
                    #         i = index1[0]
                    #         tr.stats.tmax = df_model.Tmax.iloc[i]
                    # 
                    # if tr.stats.channel in ['SHZ']:            
                    #     index1 = df_model.loc[(df_model.Event == impact) & (df_model.Station == sta) & (df_model.Comp == 'SP')].index.tolist()
                    #     if len(index1) > 0:
                    #         i = index1[0]
                    #         tr.stats.tmax = df_model.Tmax.iloc[i]
                        

                    observation_stream.append(tr)
                    
    observation_stream.sort(keys=['distance_in_degree'])
                    
    return observation_stream


def plot_spectrogram_obs(original_trace1,title=None,startsecond=0,endsecond=1800,clip=[0.0, 1.0],remove_response=False,freqmin=0.3,freqmax=0.9,inv=None):

    original_trace1_copy = original_trace1.copy()

    # make some basic preparations for the traces (remove the -1 values, remove the single sample gaps, taper, remove the mean)
    remove_negative_ones_trace(original_trace1_copy)
    # interpolate across the gaps of one sample 
    linear_interpolation(original_trace1_copy,interpolation_limit=1)

    if remove_response:
        pre_filt = [freqmin/1.5,freqmin,freqmax,freqmax*1.5]
        original_trace1_copy.remove_response(inventory=inv, pre_filt=pre_filt, zero_mean=True, taper=True, output="DISP",
                       water_level=None, plot=False)

    if endsecond is not None:
        original_trace1_copy.trim(endtime=original_trace1_copy.stats.impact_time + endsecond)
    if startsecond is not None:
        original_trace1_copy.trim(starttime=original_trace1_copy.stats.impact_time + startsecond)

    # HEREX

    # # taper
    # if observations: 
    #     epicentral_stream.taper(max_percentage=None,max_length=taper_len,side='both')
    # else: 
    #     epicentral_stream.trim(starttime=UTCDateTime(-200.0),pad=True,fill_value=0.0)
    #     epicentral_stream.taper(max_percentage=None,max_length=taper_len,side='right')
    # 
    # 
    # if pre_filt_env is not None: 
    #     epicentral_stream_env = epicentral_stream.copy()



#     if endsecond is not None:
#         original_trace1_copy = original_trace1_copy.slice(endtime=UTCDateTime(endsecond))
#     if startsecond is not None:
#         original_trace1_copy = original_trace1_copy.slice(starttime=UTCDateTime(startsecond))

    # print('{:} Distance from Source {:.1f} km'.format(station, original_trace1_copy.stats.distance_in_km))
    # original_trace1_copy.resample(lowpass)



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
#     plt.ylim(top=lowpass)
    plt.show()
                



