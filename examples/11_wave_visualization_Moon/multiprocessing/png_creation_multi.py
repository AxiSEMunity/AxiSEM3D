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
import re
import multiprocessing as mp

print('done importing')

warnings.filterwarnings("ignore", category=DeprecationWarning) 

MOON_RADIUS_IN_KM = 1737.1

# specify a run name
run = '161pre_ISSI_linear50_full_2' # to adapt to the simulation

# Specify an associated run title
# run_title = 'M1 with ±50% heterogeneity (linear to 50 km), surface explosion'
# run_title = 'Very Simple Moon, surface explosion'
# run_title = 'Very Simple Moon, deep explosion'
run_title = 'Lunar Model M1 with heterogeneity '
# run_title = 'Test Model'

# model for TauP
model_taup='homogeneous_Moon_taup' # it has no boundaries

# top level dir 
# top_dir = '/Users/mfouchet/Documents/Simulations/' # to adapt with user's directory
top_dir = '/scratch/planetseismology/mfouchet/'
folder='simu3D'

# output.txt file 
output_name = 'output.txt'
output_file = os.path.join(top_dir,run,'output.txt')
target_line = '============================ Sources ==========================='
lines_to_skip = 3 # lines to skip after the target line, do not change unless the output.txt file changes

# Testing the output.txt file
if not os.path.isfile(output_file):
    raise FileNotFoundError(f"The file '{output_name}'do not exists in the directory '{os.path.join(top_dir,run)}'.")

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

pv.global_theme.allow_empty_mesh = True

def find_lines_number(filename, target_line, lines_to_skip):
    with open(filename, 'r', encoding='utf-8') as file:
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
    with open(filename, 'r', encoding='utf-8') as f:
        line_number1, line_number2 = find_lines_number(filename, target_line, lines_to_skip)
        lines=f.readlines()
        line1 = lines[line_number1]
        line2 = lines[line_number2]
        _, _, value1 = line1.partition('=')
        _, _, value2 = line2.partition('=')
        latitude = float(value1.strip())
        longitude = float(value2.strip())
    return latitude,longitude

def process_mesh(mesh_file, mesh_dir, data_dir, output_file, target_line, lines_to_skip, run_title, clim, MOON_RADIUS_IN_KM, include_channels):

    ###############################
    
    # make a colormap
    pink1 = '#EF8683'
    blue1 = '#A9C7E9'
    lightgray = '#D3D3D3'
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [blue1,lightgray,pink1])
    
    # create dictionary of parameters to control scalar bar
    sargs = dict(title='Displacement in meters',position_x=0.21,color='white',title_font_size=40,label_font_size=50)

    # Directory    
    data_dir = os.path.join(top_dir,run,folder,'output','stations',element_name)

    # Radius of the texture (km)
    R = MOON_RADIUS_IN_KM-2
    
    # Position of the event 
    event_lat, event_lon = find_coordinates(output_file, target_line, lines_to_skip)
    event_lat = event_lat*np.pi/180
    event_lon = event_lon*np.pi/180
    event_pos = (5*R*np.cos(event_lat)*np.cos(event_lon),5*R*np.cos(event_lat)*np.sin(event_lon),5*R*np.sin(event_lat))
    
    ###############################
    
    # create images folder
    images_dir = data_dir + '/images'
    os.makedirs(images_dir, exist_ok=True)

    file_path = os.path.join(mesh_dir, mesh_file)

    # number of the mesh
    match = re.search(r'\d+', mesh_file)

    if match:
        number = int(match.group())
        
    # Initializing the plotter
    plotter = pv.Plotter(off_screen=True)

    # Loading the mesh
    interpolated_mesh = pv.read(file_path)

    #Adding the mesh to the plotter
    plotter.add_mesh(interpolated_mesh,clim=clim['Z'],cmap=cmap,scalar_bar_args=sargs,opacity=0.6)
    
    # Setting the back ground for the image
    plotter.set_background('black')

    # Setting the text for the image
    plotter.add_text(run_title, position=(0.01,0.90), color='white',font_size=40, viewport=True,font='arial')
    plotter.add_text('Vertical', position='lower_right', color='white',font_size=40, font='arial')
    # plotter.add_text('{:.02f} s'.format(data_time[itime]), position='lower_left', color='white',font_size=40, font='arial')

    # Setting the Moon texture
    sphere = pv.Sphere(radius=R, theta_resolution=300, phi_resolution=300,start_theta=270.001, end_theta=270)
    sphere.active_texture_coordinates = np.zeros((sphere.points.shape[0], 2))
    sphere.active_texture_coordinates[:, 0] = 0.5 + np.arctan2(-sphere.points[:, 0], sphere.points[:, 1])/(2 * np.pi)
    sphere.active_texture_coordinates[:, 1] = 0.5 + np.arcsin(sphere.points[:, 2]/R) / np.pi
    moon = pv.Texture('lroc_color_poles_16k.png')
    plotter.add_mesh(sphere, texture=moon, smooth_shading=False)
     
    # camera position
    plotter.camera_position = [event_pos,(0, 0, 88), (0, 1, 0)]
    plotter.camera.zoom(1)

    # Labeling the Schrodinger Bassin
    fss_lat = -71.379*np.pi/180
    fss_lon = 138.248*np.pi/180
    fss__pos = [R*np.cos(fss_lat)*np.cos(fss_lon),R*np.cos(fss_lat)*np.sin(fss_lon),R*np.sin(fss_lat)]
    sch = pv.PolyData(fss__pos)
    sch['My Labels'] = ['Schrodinger Basin']
    plotter.add_point_labels(sch, "My Labels", point_size=25,
    italic=True,font_size=60,text_color='white',point_color='black',shape_opacity=0,render_points_as_spheres=True,always_visible=True)

    # Extracting the video
    png = 'simulation_{}_{:04d}.png'.format(include_channels[0],number)
    png_file = os.path.join(images_dir,png)
    plotter.window_size = [3000, 3000]
    plotter.screenshot(png_file)
    print(png_file)
    plotter.clear()
        
    plotter.close()
    gc.collect()

    print(count)
    return count

if __name__ == '__main__':
    data_dir = os.path.join(top_dir, run, folder, 'output', 'stations', element_name)
    mesh_dir = data_dir + '/mesh'
    all_files = os.listdir(mesh_dir)
    mesh_files = [f for f in all_files if f.endswith(('.vtk', '.stl', '.obj', '.ply'))]
    mesh_files.sort()

    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.starmap(process_mesh, [(mesh_file, mesh_dir, data_dir, output_file, target_line, lines_to_skip, run_title, clim, MOON_RADIUS_IN_KM, include_channels) for mesh_file in mesh_files])
