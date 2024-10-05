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

####### PARAMETERS #######

# specify a run name
run1 = '161pre_ISSI_linear50_full_2'
run2 = '160_ISSI_2'

# model for TauP
model_taup='homogeneous_Moon_taup' # it has no boundaries

# top level dir 
# top_dir = '/Users/mfouchet/Documents/Simulations/' # to adapt with user's directory
top_dir = '/scratch/planetseismology/mfouchet/'
folder='simu3D'

# channels to calculate
include_channels = ['Z']

# Filtering parameters
freqmin = 1/100  # Minimum frequency in Hz
freqmax = 1/2 # Maximum frequency in Hz
corners = 6  # Number of corners
zerophase = False  # Apply filter in both directions
dt = 0.24675 # sampling period found in the temporal section of the output.txt
fs = 1/dt

# Moon surface threshold
S_threshold = 13000000 # in km2, less than half of the moon's hemisphere

##########################

# output.txt file 
output_file = os.path.join(top_dir,run1,'output.txt')
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
element_name='stations_array'
clim={'R': [-1e-15, 1e-15], 'T': [-1e-15, 1e-15], 'Z': [-1e-6, 1e-6]}
model_type = 'MOON'


# Function for the calculation of the mesh surface
def triangle_area(v0, v1, v2):
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

# Main processing function
def processing(top_dir=None,run1=None,run2=None,element_name=None,station_group=None,station_file=None):
    
    # map channel names to something human readable
    channel_dict = {'U1' : 'R', 'U2' : 'T', 'U3' : 'Z'}
    channel_dict2 = {'U1' : 'Radial', 'U2' : 'Transverse', 'U3' : 'Vertical'}
    
    ###############################
    # read in the netcdf file     
    data_dir = os.path.join(top_dir,run1,folder,'output','stations',element_name)
    print('about to read')
    file11 = 'axisem3d_synthetics.nc.rank_all.nc'
    file1 = os.path.join(top_dir,run1,folder,'output','stations',element_name, file11)
    print(file1)
    ds = nc4.Dataset(file1)
    
    data_dir2 = os.path.join(top_dir,run2,folder,'output','stations',element_name)
    print('about to read')
    file12 = 'axisem3d_synthetics.nc.rank_all.nc'
    file2 = os.path.join(top_dir,run2,folder,'output','stations',element_name, file12)
    print(file2)
    ds2 = nc4.Dataset(file2)
    
    for ch_inc in include_channels:
        
        ###############################
        # make the mesh from the netcdf file 
        element_lat = ds.variables['list_lat'][:]
        element_lon = ds.variables['list_lon'][:]

        element_lat2 = ds2.variables['list_lat'][:]
        element_lon2 = ds2.variables['list_lon'][:]

        # print("Matrices lat are equal:", np.array_equal(element_lat,element_lat2))
        # print("Matrices lon are equal:", np.array_equal(element_lon,element_lon2))
        
        wave11 = ds.variables['data_wave'][:, :, :]
        # print(np.shape(wave11))
        # print('wave11 read')
        wave12 = ds2.variables['data_wave'][:, :, :]
        # print(np.shape(wave12))
        # print('wave12 read')
        # wave1 = wave11 - wave12
        # print('wave1 read')
        dim_data_wave = np.shape(wave11)
        
        # read time
        data_time = ds.variables['data_time'][:]
        
        # time steps
        ntime = len(data_time)

        # read the list of stations
        list_station = np.array(ds.variables['list_station'][:])
        decoded_list = [''.join(b.decode('utf-8') for b in i).strip() for i in list_station]
        
        list_station2 = np.array(ds2.variables['list_station'][:])
        decoded_list2 = [''.join(b.decode('utf-8') for b in i).strip() for i in list_station2]

        # create mesh folder
        mesh_dir = data_dir + '/mesh'
        os.makedirs(mesh_dir, exist_ok=True)
         
        #############################
        # make the mesh          

        #  Element coordinates in cartesian
        R = MOON_RADIUS_IN_KM
        element_coords = np.ndarray((dim_data_wave[0], dim_data_wave[1]))
        element_lat2 = element_lat2*np.pi/180
        element_lon2 = element_lon2*np.pi/180
        element_coords[:,0] = R*np.cos(element_lat2)*np.cos(element_lon2)
        element_coords[:,1] = R*np.cos(element_lat2)*np.sin(element_lon2)
        element_coords[:,2] = R*np.sin(element_lat2)

        # read in the wave simulation
        for wave_dim, channel1 in enumerate(ds['list_channel']):

            # wave dimension to animate
            channel = str(chartostring(channel1))
            wave_channel = channel_dict[channel]  

            if wave_channel != ch_inc:
                continue

            ch = 'U' + wave_channel

            # setting up the data
            wave11 = ds.variables['data_wave'][:, :, :]
            print('wave11 read')
            # print('filtering...')
            # wave12 = bandpass(wave11, freqmin, freqmax, fs, corners=corners, zerophase=zerophase)

            
            wave21 = ds2.variables['data_wave'][:, :, :]
            print('wave21 read')
            # print('filtering...')
            # wave22 = bandpass(wave21, freqmin, freqmax, fs, corners=corners, zerophase=zerophase)

            # print("Matrices are equal:", np.array_equal(wave11,wave12))

            # Step 1: Create a dictionary for quick lookup
            index_map = {station: idx for idx, station in enumerate(decoded_list)}
            
            # Step 2: Find indices of elements of A in B
            indices = [index_map.get(station) for station in decoded_list2]
            
            output_matrix = wave11.copy()

            print('joining the matrix')
            for i, idx in enumerate(indices):
                output_matrix[i] = wave11[idx]

            print('joining done')
            wave2 = np.linalg.norm(output_matrix - wave21,axis=1)
            print(np.shape(wave2))
            # print('filtering...')
            # wave2 = bandpass(wave22, freqmin, freqmax, fs, corners=corners, zerophase=zerophase)
            
            # Apply the bandpass filter 
            # print('filtering...')
            # wave2 = bandpass(wave1, freqmin, freqmax, fs, corners=corners, zerophase=zerophase)
        
            print('setup complete')

            # Initializing the surface thresholds
            S_mesh = 0
            
            # loop over time
            for itime in range(1342,ntime):          
            # for itime in [2726]: 

                # plotter = pyvistaqt.BackgroundPlotter()

                # Analyzing the non zero data
                points = element_coords[:,:]
                wave3 = wave2[:,itime]

                # points2 = element_coords2[:,:]
                # wave23 = wave22[:,itime]

                # mesh1 = pv.PolyData(points)
                # mesh1[ch] = wave13[:]

                # mesh2 = pv.PolyData(points2)
                # mesh2[ch] = wave23[:]

                # wave = wave13[:]-wave23[:]
                
                # mesh = pv.PolyData(points)
                # mesh[ch] = mesh1[ch] - mesh2[ch]
                
                # if S_mesh <= S_threshold :
                #     threshold = 0.4e-7
                # else:
                #     threshold = 0
                threshold = 0.4e-7
                nonzeroind = np.where(np.abs(wave3) >= threshold)[0]
                points = points[nonzeroind]

                # nonzeroind2 = np.where(np.abs(wave23[:]) >= threshold)[0]
                # points2 = points2[nonzeroind2]
                
                if len(points) >= 6:
                    mesh = pv.PolyData(points)   
                    mesh[ch] = wave3[nonzeroind] 

                    # mesh2 = pv.PolyData(points2)   
                    # mesh2[ch] = wave23[nonzeroind2] 

                    # if mesh1.n_points == mesh2.n_points:
                    #     # Subtract the values
                    #     subtracted_values = mesh1[ch] - mesh2[ch]
                    
                    #     # Create a result mesh with the same geometry as mesh1
                    #     result_mesh = mesh1.copy()  # Copy mesh1 to retain geometry
                    #     result_mesh[ch] = subtracted_values
                    # else:
                    #     print('error in the number of points')
                    
                    # open3d mesh creation
                    point_cloud = points
                    # point_cloud2 = points2
                    
                    pcd = o3d.geometry.PointCloud()
                    # pcd2 = o3d.geometry.PointCloud()
                    
                    pcd.points = o3d.utility.Vector3dVector(point_cloud[:,:3])
                    # pcd2.points = o3d.utility.Vector3dVector(point_cloud2[:,:3])
                    
                    # print('points ok')
                    pcd.estimate_normals()
                    # pcd2.estimate_normals()
                    
                    # print('normals ok')
                    pcd.orient_normals_consistent_tangent_plane(5)
                    # pcd2.orient_normals_consistent_tangent_plane(5)
                    
                    # print('orientation ok')
                    distances = pcd.compute_nearest_neighbor_distance()
                    # distances2 = pcd2.compute_nearest_neighbor_distance()
                    
                    # print('ok mesh open3d')
                    avg_dist = np.mean(distances)
                    # avg_dist2 = np.mean(distances2)
                    
                    # radius = 3 * avg_dist
                    # bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector([radius, radius * 2]))
                    radii = [avg_dist/2,avg_dist,3*avg_dist,8*avg_dist,25*avg_dist]
                    # radii2 = [avg_dist2/2,avg_dist2,3*avg_dist2,8*avg_dist2,25*avg_dist2]
                    
                    bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector(radii))
                    # bpa_mesh2 = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd2,o3d.utility.DoubleVector(radii2))
                    
                    vertices = np.asarray(bpa_mesh.vertices)
                    # vertices2 = np.asarray(bpa_mesh2.vertices)
                    
                    triangles = np.asarray(bpa_mesh.triangles)
                    # triangles2 = np.asarray(bpa_mesh2.triangles)
                    
                    S_mesh = sum(triangle_area(vertices[t[0]], vertices[t[1]], vertices[t[2]]) for t in triangles)
                    # S_mesh2 = sum(triangle_area(vertices2[t[0]], vertices2[t[1]], vertices2[t[2]]) for t in triangles2)
                    
                    filename = "bpa_mesh.ply"
                    o3d.io.write_triangle_mesh(filename, bpa_mesh)
                    bpa_mesh12 = pv.read(filename)
                    interpolated = bpa_mesh12.interpolate(mesh, radius=10,sharpness=10)
                    
                    # filename2 = "bpa_mesh2.ply"
                    # o3d.io.write_triangle_mesh(filename2, bpa_mesh2)
                    # bpa_mesh22 = pv.read(filename2)
                    # interpolated2 = bpa_mesh22.interpolate(mesh2, radius=10,sharpness=10)

                    # result_mesh = mesh1.copy()
                    # result_mesh[ch] = interpolated[ch] - interpolated2[ch]
                    
                else :
                    interpolated = pv.PolyData()

                print('done meshing')

                
                
                mesh_name = 'mesh_{:04d}.vtk'.format(itime)
                mesh_file = os.path.join(mesh_dir,mesh_name)
                interpolated.save(mesh_file)
                # result_mesh.save(mesh_file)

                print('mesh_' + str(itime) + ' saved')
               

        print('\nDone %s' % wave_channel)        
        print(data_time[itime])
    return mesh, data_time[itime]

if __name__ == '__main__':
    mesh, sample_time = processing(top_dir=top_dir,
            run1=run1,run2=run2,element_name=element_name)