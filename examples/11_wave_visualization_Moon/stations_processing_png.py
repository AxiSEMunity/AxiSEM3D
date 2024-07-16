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

warnings.filterwarnings("ignore", category=DeprecationWarning) 

MOON_RADIUS_IN_KM = 1737.1

# specify a run name
run = '163_very_simple_moonShroedinger_hemsiphere_10' # to adapt to the simulation

# run_title = 'M1 with ±50% heterogeneity (linear to 50 km), surface explosion'
# run_title = 'Very Simple Moon, surface explosion'
# run_title = 'Very Simple Moon, deep explosion'
# run_title = 'Lunar Model M1 - no heterogeneity '
run_title = 'Test Model'

# model for TauP
model_taup='homogeneous_Moon_taup' # it has no boundaries

# top level dir 
top_dir = '/Users/mfouchet/Documents/Simulations/' # to adapt with user's directory
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

def find_lines_number(filename, target_line, lines_to_skip):
    with open(filename, 'r') as file:
        for line_number, line in enumerate(file, 1):
            if target_line in line:
                for line_number_bis, line_bis in enumerate(file, line_number):
                    
                    if line_number_bis-line_number<= lines_to_skip:
                        continue
                    line1 = line_bis
                    line_number1 = line_number_bis
                    line_number2 = line_number_bis + 1
                    break
                return line_number1,line_number2
    return None, None

def find_coordinates(filename, target_line, lines_to_skip):
    with open(filename, 'r') as f:
        line_number1, line_number2 = find_lines_number(filename, target_line, lines_to_skip)
        lines=f.readlines()
        line1 = lines[line_number1]
        line2 = lines[line_number2]
        _, _, value1 = line1.partition('=')
        _, _, value2 = line2.partition('=')
        latitude = float(value1.strip())
        longitude = float(value2.strip())
    return latitude,longitude

def triangle_area(v0, v1, v2):
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

def animate_pyvista_png(top_dir=None,run=None,element_name=None,station_group=None,station_file=None, run_title=None):

    # make a colormap
    pink1 = '#EF8683'
    blue0 = '#D0EFFF'
    blue1 = '#A9C7E9'
    blue2 = '#021076'
    black = '#000000'
    yellow = '#FFD966'
    red1 = '#E50000'
    orange1 = '#FFA500'
    lightgray = '#D3D3D3'
    transparent = '#03D3D3D3'
    white = '#FFFFFF'
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [blue1,lightgray,pink1])
    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [blue0,blue1,blue2,black,red1,orange1,yellow])

    # create dictionary of parameters to control scalar bar
    sargs = dict(title='Displacement in meters',position_x=0.21,color='white',title_font_size=40,label_font_size=50)
    
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
 
        # create images folder
        images_dir = data_dir + '/images'
        os.makedirs(images_dir, exist_ok=True)
         
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
            
            # Define the bandpass filter parameters
            freqmin = 1/100  # Minimum frequency in Hz
            freqmax = 1/2  # Maximum frequency in Hz
            corners = 6  # Number of corners
            zerophase = False  # Apply filter in both directions
            dt = 0.0967807
            fs = 1/dt
            
            # Apply the bandpass filter 
            print('filtering...')
            wave2 = bandpass(wave1, freqmin, freqmax, fs, corners=corners, zerophase=zerophase)
            
            # Initializing the plotter
            plotter = pv.Plotter(off_screen=True)
            print('setup complete')

            # Initializing the surface thresholds
            S_mesh = 0
            S_threshold = 13000000
            
            # loop over time
            for itime in range(0,ntime):          
            # for itime in [2726]: 

                # plotter = pyvistaqt.BackgroundPlotter()

                # Analyzing the non zero data
                points = element_coords[:,:]
                wave3 = wave2[:,itime]

                mesh = pv.PolyData(points)
                mesh[ch] = wave3[:]
                
                # if itime <= np.floor(ntime/3) :
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
                    point_cloud= points
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
                    plotter.add_mesh(interpolated,clim=clim['Z'],cmap=cmap,scalar_bar_args=sargs,opacity=0.6)
                else :
                    mesh = 0

                print('done meshing')
                # Setting the back ground for the image
                # image_path = examples.planets.download_stars_sky_background(load=True)
                # plotter.add_background_image(image_path)
                plotter.set_background('black')

                # Setting the text for the image
                plotter.add_text(run_title, position=(0.01,0.90), color='white',font_size=40, viewport=True,font='arial')
                plotter.add_text(channel_dict2[channel], position='lower_right', color='white',font_size=40, font='arial')
                plotter.add_text('{:.02f} s'.format(data_time[itime]), position='lower_left', color='white',font_size=40, font='arial')

                # Setting the Moon texture
                R = MOON_RADIUS_IN_KM-2
                sphere = pv.Sphere(radius=R, theta_resolution=300, phi_resolution=300,start_theta=270.001, end_theta=270)
                sphere.active_texture_coordinates = np.zeros((sphere.points.shape[0], 2))
                sphere.active_texture_coordinates[:, 0] = 0.5 + np.arctan2(-sphere.points[:, 0], sphere.points[:, 1])/(2 * np.pi)
                sphere.active_texture_coordinates[:, 1] = 0.5 + np.arcsin(sphere.points[:, 2]/R) / np.pi
                moon = pv.Texture('lroc_color_poles_16k.png')
                plotter.add_mesh(sphere, texture=moon, smooth_shading=False)
                 
                # camera position
                event_lat, event_lon = find_coordinates(output_file, target_line, lines_to_skip)
                event_lat = event_lat*np.pi/180
                event_lon = event_lon*np.pi/180
                event_pos = (5*R*np.cos(event_lat)*np.cos(event_lon),5*R*np.cos(event_lat)*np.sin(event_lon),5*R*np.sin(event_lat))
                plotter.camera_position = [event_pos,(0, 0, 88), (0, 1, 0)]
                plotter.camera.zoom(1)

                # Labeling the Schrodinger Bassin
                fss_lat = -71.379*np.pi/180
                fss_lon = 138.248*np.pi/180
                fss__pos = [R*np.cos(fss_lat)*np.cos(fss_lon),R*np.cos(fss_lat)*np.sin(fss_lon),R*np.sin(fss_lat)]
                sch = pv.PolyData(fss__pos)
                sch["My Labels"] = ['Schrodinger Basin']
                plotter.add_point_labels(sch, "My Labels", point_size=25,
                italic=True,font_size=65,text_color='white',point_color='black',shape_opacity=0,render_points_as_spheres=True,always_visible=True)

                # Extracting the video
                png = 'simulation_{}_{:04d}.png'.format(wave_channel,itime)
                png_file = os.path.join(images_dir,png)
                plotter.window_size = [3000, 3000]
                plotter.screenshot(png_file)
                print(png_file)
                # plotter.deep_clean()
                plotter.clear()
                
            plotter.close()
            gc.collect()

        print('\nDone %s' % wave_channel)        
        print(png_file)
        print(data_time[itime])
    return mesh, data_time[itime]

if __name__ == '__main__':
    mesh, sample_time = animate_pyvista_png(top_dir=top_dir,
            run=run,element_name=element_name,run_title=run_title)