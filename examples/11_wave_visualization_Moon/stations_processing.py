#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Process stations of AxiSEM3D model and make png files for surface projection
Author: Matheo Fouchet, JPL
"""

import os
import pyvista as pv
from pyvista import examples
# import pyvistaqt 
import vtk
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as nc4
from netCDF4 import chartostring
import matplotlib
from datetime import datetime
import sys
import open3d as o3d
import warnings
import gc
from obspy import read
from obspy.signal.filter import bandpass,lowpass

# pyvista.start_xvfb(wait=0.1)

warnings.filterwarnings("ignore", category=DeprecationWarning) 

MOON_RADIUS_IN_KM = 1737.1

# specify a run name
# run = '157w_ISSI_atten_linear50_10' # to adapt to the simulation
run = '161pre_ISSI_linear50_full_2'

# run_title = 'M1 with ±50% heterogeneity (linear to 50 km), surface explosion'
# run_title = 'Very Simple Moon, surface explosion'
# run_title = 'Very Simple Moon, deep explosion'
run_title = 'Lunar Model M1 with heterogeneity '
# run_title = 'Test Model'

# model for TauP
model_taup='homogeneous_Moon_taup' # it has no boundaries

# top level dir 
top_dir = '/Users/mfouchet/Documents/Simulations/' # to adapt with user's directory
# top_dir = '/scratch/planetseismology/mfouchet/'
folder='simu3D'

# output.txt file 
output_file = os.path.join(top_dir,run,'output.txt')
target_line = '============================ Sources ==========================='
lines_to_skip = 3 # lines to skip after the target line, do not change unless the output.txt file changes

#######
# specify station_group
station_group = 'stations_array'
# specify station file
station_file = 'stations_array.txt'
# specify a station key (network.name)
station_key = 'A.A0'

# source name 
source_name='explosion_south_east_quadrant' # this is an incorrect name, depth=0kms
element_name='stations_array'
clim={'R': [-1e-15, 1e-15], 'T': [-1e-15, 1e-15], 'Z': [-1e-6, 1e-6]}
model_type = 'MOON'

# channels to calculate
# include_channels = ['R', 'T', 'Z']
# include_channels = ['R','Z']
# include_channels = ['R']
# include_channels = ['T']
include_channels = ['Z']

# Filtering parameters
freqmin = 1/100  # Minimum frequency in Hz
freqmax = 1/2  # Maximum frequency in Hz
corners = 6  # Number of corners
zerophase = False  # Apply filter in both directions
dt = 0.0967807
fs = 1/dt

# Function for the calculation of the mesh surface
def triangle_area(v0, v1, v2):
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

# Main processing function
def processing(top_dir=None,run=None,element_name=None,station_group=None,station_file=None, run_title=None):
    
    # map channel names to something human readable
    channel_dict = {'U1' : 'R', 'U2' : 'T', 'U3' : 'Z'}
    channel_dict2 = {'U1' : 'Radial', 'U2' : 'Transverse', 'U3' : 'Vertical'}
    
    ###############################
    # read in the netcdf file     
    data_dir = os.path.join(top_dir,run,folder,'output','stations',element_name)
    print('about to read')
    file1 = 'axisem3d_synthetics.nc.rank_all.nc'
    file = os.path.join(top_dir,run,folder,'output','stations',element_name, file1)
    print(file)

    ds = nc4.Dataset(file)

    for ch_inc in include_channels:
        
        ###############################
        # make the mesh from the netcdf file 
        element_lat = ds.variables['list_lat'][:]
        element_lon = ds.variables['list_lon'][:]

        wave1 = ds.variables['data_wave'][:, :, :]
        dim_data_wave = np.shape(wave1)
        
        # read time
        data_time = ds.variables['data_time'][:]
        
        # time steps
        ntime = len(data_time)

        # read the list of stations
        list_station = ds.variables['list_station'][:,:]
 
        # create mesh folder
        mesh_dir = data_dir + '/mesh'
        os.makedirs(mesh_dir, exist_ok=True)
         
        #############################
        # make the mesh          

        #  Element coordinates in cartesian
        R = MOON_RADIUS_IN_KM
        element_coords = np.ndarray((dim_data_wave[0], dim_data_wave[1]))
        element_lat = element_lat*np.pi/180
        element_lon = element_lon*np.pi/180
        element_coords[:,0] = R*np.cos(element_lat)*np.cos(element_lon)
        element_coords[:,1] = R*np.cos(element_lat)*np.sin(element_lon)
        element_coords[:,2] = R*np.sin(element_lat)

        # read in the wave simulation
        for wave_dim, channel1 in enumerate(ds['list_channel']):

            # wave dimension to animate
            channel = str(chartostring(channel1))
            wave_channel = channel_dict[channel]  

            if wave_channel != ch_inc:
                continue

            ch = 'U' + wave_channel

            # setting up the data
            wave1 = ds.variables['data_wave'][:, wave_dim, :]
            
            # Apply the bandpass filter 
            print('filtering...')
            wave2 = bandpass(wave1, freqmin, freqmax, fs, corners=corners, zerophase=zerophase)
        
            print('setup complete')

            # Initializing the surface thresholds
            S_mesh = 0
            S_threshold = 13000000
            
            # loop over time
            for itime in range(1586,ntime):          
            # for itime in [2726]: 

                # plotter = pyvistaqt.BackgroundPlotter()

                # Analyzing the non zero data
                points = element_coords[:,:]
                wave3 = wave2[:,itime]

                mesh = pv.PolyData(points)
                mesh[ch] = wave3[:]
                
                if S_mesh <= S_threshold :
                    threshold = 0.4e-7
                else:
                    threshold = 0
                nonzeroind = np.where(np.abs(wave3[:]) >= threshold)[0]
                points = points[nonzeroind]
                
                if len(points) >= 6:
                    mesh = pv.PolyData(points)   
                    mesh[ch] = wave3[nonzeroind] 
                    
                    # open3d mesh creation
                    point_cloud = points
                    pcd = o3d.geometry.PointCloud()
                    pcd.points = o3d.utility.Vector3dVector(point_cloud[:,:3])
                    # print('points ok')
                    pcd.estimate_normals()
                    # print('normals ok')
                    pcd.orient_normals_consistent_tangent_plane(5)
                    # print('orientation ok')
                    distances = pcd.compute_nearest_neighbor_distance()
                    # print('ok mesh open3d')
                    avg_dist = np.mean(distances)
                    radius = 3 * avg_dist
                    bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector([radius, radius * 2]))
                    vertices = np.asarray(bpa_mesh.vertices)
                    triangles = np.asarray(bpa_mesh.triangles)
                    S_mesh = sum(triangle_area(vertices[t[0]], vertices[t[1]], vertices[t[2]]) for t in triangles)
                    filename = "bpa_mesh.ply"
                    o3d.io.write_triangle_mesh(filename, bpa_mesh)
                    bpa_mesh2 = pv.read(filename)
                    interpolated = bpa_mesh2.interpolate(mesh, radius=10,sharpness=10)
                    
                else :
                    interpolated = pv.PolyData()

                print('done meshing')
                
                mesh_name = 'mesh_{:04d}.vtk'.format(itime)
                mesh_file = os.path.join(mesh_dir,mesh_name)
                interpolated.save(mesh_file)

                print('mesh_' + str(itime) + ' saved')
               

        print('\nDone %s' % wave_channel)        
        print(data_time[itime])
    return mesh, data_time[itime]

if __name__ == '__main__':
    mesh, sample_time = processing(top_dir=top_dir,
            run=run,element_name=element_name,run_title=run_title)