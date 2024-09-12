#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Process stations of AxiSEM3D model and make png files for surface projection
Author: Matheo Fouchet, JPL
"""

import os
import pyvista as pv
from pyvista import examples
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
from obspy.taup import TauPyModel
import time
import cv2
import matplotlib.colors as mcolors

calculated_VPREMOON = TauPyModel(model='VPREMOON_atten_taup')
calculated_ISSI_M1 = TauPyModel(model='ISSI_MOON_M1_atten_taup')
calculated_VPREMOON_atten_no_LVZ_taup = TauPyModel(model='VPREMOON_atten_no_LVZ_taup')

calculated_taup_dict =	{
  "VPREMOON_taup": calculated_VPREMOON,
  "ISSI_MOON_M1_taup": calculated_ISSI_M1,
    'VPREMOON_atten_no_LVZ_taup':calculated_VPREMOON_atten_no_LVZ_taup
}

from obspy.core import Stream, UTCDateTime
from obspy import read_inventory

from postprocessing_util import get_all_streams_from_netcdf
from postprocessing_util import plot_epicentral_distance_taup
from postprocessing_util import plot_envelope_taup

warnings.filterwarnings("ignore", category=DeprecationWarning) 

MOON_RADIUS_IN_KM = 1737.1

################# PARAMETERS #################

#### specify a run name ###
# run = '157_ISSI_atten_linear50_slice_10' # to adapt to the simulation
run = '161pre_ISSI_linear50_full_2'
# run = '158_ISSI_atten_slice_10'
# run = '160_ISSI_2'

#### specify a run title ####
# run_title = 'Lunar Model M1 without heterogeneity, surface explosion'
# run_title = 'Very Simple Moon, surface explosion'
# run_title = 'Very Simple Moon, deep explosion'
# run_title = 'Lunar Model M1 with heterogeneity min period 10.49'
run_title = '3D-1D Model M1 with heterogeneity'

#### specify a short title and model TauP ####
short_title = 'Model M1'
model_taup='ISSI_MOON_M1_taup' # it has no boundaries

#### specify top level dir and folder ####
top_dir = '/Users/mfouchet/Documents/Simulations/' # to adapt with user's directory
# top_dir = '/scratch/planetseismology/mfouchet/'
folder='simu3D'

#### specify camera parameter ####
pos_cam = 'tilted' # position of the camera can be either straight or tilted

#### specify channels to calculate ####
# include_channels = ['R', 'T', 'Z']
# include_channels = ['R','Z']
# include_channels = ['R']
# include_channels = ['T']
include_channels = ['Z']

#### specify borders for the colourscale ####
# clim={'R': [-1e-15, 1e-15], 'T': [-1e-15, 1e-15], 'Z': [-1e-6, 1e-6]}
clim={'R': [-1e-15, 1e-15], 'T': [-1e-15, 1e-15], 'Z': [0, 1e-6]}

################# END PARAMETERS #################

# output.txt file 
output_name = 'output.txt'
output_file = os.path.join(top_dir,run,'output.txt')
target_line = '============================ Sources ==========================='
lines_to_skip = 3 # lines to skip after the target line, do not change unless the output.txt file changes

# Testing the output.txt file
if not os.path.isfile(output_file):
    raise FileNotFoundError(f"The file '{output_name}'do not exists in the directory '{os.path.join(top_dir,run)}'.")

# specify station_group
station_group = 'stations_array'
new_station_group = 'seismograms_stations'
# specify station file
station_file = 'stations_array.txt'
# specify a station key (network.name)
station_key = 'A.A0'

# source name 
source_name='explosion_south_east_quadrant' # this is an incorrect name, depth=0kms
element_name='stations_array'
model_type = 'MOON'

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

def sismometer_observation(top_dir,run,folder,short_title,new_station_group,station_group):
    # create new station folder
    data_dir = os.path.join(top_dir, run, folder, 'output', 'stations')
    seismo_dir = data_dir + '/seismograms_stations'
    os.makedirs(seismo_dir, exist_ok=True)
    fn_file = os.path.join(top_dir, run, folder, 'output', 'stations', station_group, 'axisem3d_synthetics.nc.rank_all.nc')
    
    ds = nc4.Dataset(fn_file)
    
    # Open the original NetCDF4 file
    dst_file = os.path.join(top_dir, run, folder, 'output', 'stations', new_station_group, 'axisem3d_synthetics.nc.rank_all.nc')
    
    # Create the destination file
    dst = nc4.Dataset(dst_file, 'w', format='NETCDF4')
    
    dst.createDimension('new_dim_station', 1)
    dst.createDimension('new_dim_channel', ds.dimensions['dim_channel'].size)
    dst.createDimension('new_dim_time', ds.dimensions['dim_time'].size)
    dst.createDimension('new_dim_station_str_length', ds.dimensions['dim_station_str_length'].size)
    dst.createDimension('new_dim_channel_str_length', ds.dimensions['dim_channel_str_length'].size)
    
    # Target station to find
    target_station = 'A.A0'
    
    # Read the list_station variable and decode it
    list_station = np.array(ds.variables['list_station'][:])
    decoded_list = [''.join(b.decode('utf-8') for b in i).strip() for i in list_station]
    
    # Find the index of the target station
    target_station = target_station.strip()
    index = np.where(np.array(decoded_list) == target_station)[0]
    
    # print(f"Index of the target station '{target_station}':", index)
    
    for var_name in ds.variables:
        variable = ds.variables[var_name]
        # print(variable.dimensions)
        # print(var_name)
    
        if var_name == 'data_time' or var_name == 'list_channel':
            
            # Check if dimensions exist in the destination file, if not, create them
            for dim_name in variable.dimensions:
                if dim_name not in dst.dimensions:
                    dst.createDimension(dim_name, len(ds.dimensions[dim_name]) if not ds.dimensions[dim_name].isunlimited() else None)
            
            # Create the variable in the destination file with the same name, type, and dimensions
            new_var = dst.createVariable(var_name, variable.datatype, variable.dimensions)
        
            # Copy variable attributes
            for attr_name in variable.ncattrs():
                new_var.setncattr(attr_name, variable.getncattr(attr_name))
                
                # Copy variable data
            new_var[:] = variable[:]
    
        else:
            # Create the variable in the destination file with the new dimension
            if var_name == 'list_lat' or var_name == 'list_lon':    
                new_var = dst.createVariable(var_name, variable.datatype,('new_dim_station'))
                
                # Copy variable attributes
                for attr_name in variable.ncattrs():
                    new_var.setncattr(attr_name, variable.getncattr(attr_name))
                new_var[:] = variable[index]
            
            elif var_name == 'list_station':
                new_var = dst.createVariable(var_name, variable.datatype,('new_dim_station','new_dim_station_str_length'))
                
                # Copy variable attributes
                for attr_name in variable.ncattrs():
                    new_var.setncattr(attr_name, variable.getncattr(attr_name))
                new_var[:] = variable[index,:]
            else:
                new_var = dst.createVariable(var_name, variable.datatype,('new_dim_station','new_dim_channel','new_dim_time'))
                
                # Copy variable attributes
                for attr_name in variable.ncattrs():
                    new_var.setncattr(attr_name, variable.getncattr(attr_name))
                new_var[:] = variable[index,:,:]
                
    # get all the streams from the netcdf file
    original_stream_dict = {}
    original_stream_run, st_rank_all = get_all_streams_from_netcdf(top_dir=top_dir,
    run=run,short_title=short_title,station_group=new_station_group,station_file=station_file,
    source_name=source_name,folder=folder,model_taup=calculated_taup_dict[model_taup])
    original_stream_dict[run] =  original_stream_run

    # Filtering
    freqmin=0.1
    freqmax=0.5
    taper_len = 10
    
    epicentral_stream = Stream()
    for tr in original_stream_dict[run]:
        if tr.stats.station == 'A0':
            max_data = np.max(abs(tr.data))
            print(tr.stats)
            tr.data = tr.data/max_data
            epicentral_stream.append(tr)
            # epicentral_stream.trim(starttime=UTCDateTime(-200.0),pad=True,fill_value=0.0)
            epicentral_stream.taper(max_percentage=None,max_length=taper_len,side='right')
            # epicentral_stream.filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=False)

    time = epicentral_stream[2].times()
    data = epicentral_stream[2].data

    # removing the file
    if os.path.exists(dst_file):
        os.remove(dst_file)
        
    return time, data
        
def animate_pyvista_png(top_dir=None,run=None,element_name=None,short_title=None,station_group=None,station_file=None, run_title=None,new_station_group=None):

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
    colors = ['#EF8683', '#D0EFFF', '#A9C7E9', '#021076', '#000000', 
          '#FFD966', '#E50000', '#FFA500', '#D3D3D3', '#FFFFFF']
    # Define the start and end colors
    start_color = '#E50000'  # Red
    end_color = '#021076'    # Blue

    # # Generate 10 colors between red and blue
    # colors = [mcolors.to_hex(c) for c in plt.cm.RdBu(np.linspace(0, 1, 15))]

    # # Create a colormap from the list of colors
    # cmap = mcolors.LinearSegmentedColormap.from_list('red_to_blue_cmap', colors, N=15)

    # # Generate 15 colors between pale purple to bright purple using the 'Purples' colormap
    # colors = [mcolors.to_hex(c) for c in plt.cm.Reds_r(np.linspace(0, 1, 15))]
    
    # # Create a sequential colormap from the list of colors
    # cmap = mcolors.LinearSegmentedColormap.from_list('pale_to_bright_cmap', colors, N=15)
    

    # Define the start and end colors
    start_color = '#cc0000'   # Larger, pale red
    end_color = '#ffe6e6' #'#ff9999'   # Smaller, bright red
    
    # Generate a continuous gradient of colors between start_color and end_color
    colors = [start_color, end_color]
    cmap = mcolors.LinearSegmentedColormap.from_list('pale_to_bright_red', colors, N=256)


    
    # Create a colormap from the list of colors
    
    # cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=10)

    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [blue1,lightgray,pink1])
    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [blue0,blue1,blue2,black,red1,orange1,yellow])

    # Sismograms data
    time, data = sismometer_observation(top_dir,run,folder,short_title,new_station_group,station_group)
    
    # create dictionary of parameters to control scalar bar
    sargs = dict(title='Displacement in meters',position_x=0.21,color='white',title_font_size=40,label_font_size=50)
    
    # map channel names to something human readable
    channel_dict = {'U1' : 'R', 'U2' : 'T', 'U3' : 'Z'}
    channel_dict2 = {'U1' : 'Radial', 'U2' : 'Transverse', 'U3' : 'Vertical'}
    
    ###############################
    # Directory    
    data_dir = os.path.join(top_dir,run,folder,'output','stations',element_name)
    print('about to read')
    file1 = 'axisem3d_synthetics.nc.rank_all.nc'
    file = os.path.join(top_dir,run,folder,'output','stations',element_name, file1)
    print(file)
    ds = nc4.Dataset(file)
    data_time = ds.variables['data_time'][:]
    
    # Radius of the texture (km)
    R = MOON_RADIUS_IN_KM-2
    
    # Position of the event 
    event_lat, event_lon = find_coordinates(output_file, target_line, lines_to_skip)
    event_lat = event_lat*np.pi/180
    event_lon = event_lon*np.pi/180
    event_pos = (5*R*np.cos(event_lat)*np.cos(event_lon),5*R*np.cos(event_lat)*np.sin(event_lon),5*R*np.sin(event_lat))
    
    # create images folder
    images_dir = data_dir + '/images'
    os.makedirs(images_dir, exist_ok=True)

    # create images folder
    result_dir = data_dir + '/result'
    # result_dir = data_dir + '/coda'
    os.makedirs(result_dir, exist_ok=True)
    
    # create sismo folder
    sismo_dir = data_dir + '/sismo'
    os.makedirs(sismo_dir, exist_ok=True)
    
    # Initializing the plotter
    plotter = pv.Plotter(off_screen=True)

    #Adding background
    image_path = examples.planets.download_stars_sky_background(load=False)
    # plotter.add_background_image(image_path)

    if pos_cam == 'straight':
        plotter.camera.zoom(1)
    elif pos_cam == 'tilted':
        plotter.camera.zoom(2)
    else:
        print('ERROR: camera postion has to be straight or tilted')
        sys.exit()
    print('setup complete')

    # List all files in the directory
    mesh_dir = data_dir + '/mesh'
    all_files = os.listdir(mesh_dir)

    # Filter the list to include only supported mesh files
    mesh_files = [f for f in all_files if f.endswith(('.vtk', '.stl', '.obj', '.ply'))]
    mesh_files.sort()
    # List to hold the loaded meshes
    meshes = []

    # Load each mesh file
    for mesh_file in mesh_files:

        plotter.clear()
        
        file_path = os.path.join(mesh_dir, mesh_file)

        # number of the mesh
        match = re.search(r'\d+', mesh_file)

        if match:
            number = int(match.group())
        
        # Loading the mesh
        interpolated_mesh = pv.read(file_path)
        
        #Adding the mesh to the plotter
        plotter.add_mesh(interpolated_mesh,clim=clim['Z'],cmap=cmap,scalar_bar_args=sargs,opacity=0.6)

        plotter.set_background('black')
        
        # Setting the text for the image
        plotter.add_text(run_title, position=(0.01,0.90), color='white',font_size=40, viewport=True,font='arial')
        plotter.add_text('Vertical', position='lower_right', color='white',font_size=40, font='arial')
        text_actor = plotter.add_text('{:.02f} s'.format(data_time[number]), position='lower_left', color='white',font_size=40, font='arial')

        # Setting the Moon texture
        sphere = pv.Sphere(radius=R, theta_resolution=300, phi_resolution=300,start_theta=270.001, end_theta=270)
        sphere.active_texture_coordinates = np.zeros((sphere.points.shape[0], 2))
        sphere.active_texture_coordinates[:, 0] = 0.5 + np.arctan2(-sphere.points[:, 0], sphere.points[:, 1])/(2 * np.pi)
        sphere.active_texture_coordinates[:, 1] = 0.5 + np.arcsin(sphere.points[:, 2]/R) / np.pi
        moon = pv.Texture('lroc_color_poles_16k.png')
        # moon = pv.Texture('Moon.png')
        plotter.add_mesh(sphere, texture=moon, smooth_shading=False)

        # camera position
        if pos_cam == 'straight':
            plotter.camera_position = [event_pos,(0, 0, 88), (0, 1, 0)]
        elif pos_cam == 'tilted':
            plotter.camera_position = [[event_pos[0],event_pos[1]-7000,event_pos[2]-7000],(0, 0, 88), (0, 1, 0)]
        
        # Add a scene light (fixed in the scene)
        light = plotter.add_light(pv.Light(light_type='scenelight', position=event_pos, color='white', intensity=1.0))
        
        # Labeling the Schrodinger Bassin
        fss_lat = -71.379*np.pi/180
        fss_lon = 138.248*np.pi/180
        fss__pos = [R*np.cos(fss_lat)*np.cos(fss_lon),R*np.cos(fss_lat)*np.sin(fss_lon),R*np.sin(fss_lat)]
        sch = pv.PolyData(fss__pos)
        sch['My Labels'] = ['Schrodinger Basin']
        plotter.add_point_labels(sch, "My Labels", point_size=25,
        italic=True,font_size=60,text_color='white',point_color='black',shape_opacity=0,render_points_as_spheres=True,always_visible=True)

        # Plotting the sismograms
        plt.figure()
        plt.plot(time,data)
        plt.plot(time[0:number],data[0:number],color='red')
        plt.axis('off')
        plot_file = os.path.join(sismo_dir,"sismo_plot.png")
        plt.savefig(plot_file,transparent=True)
        plt.close()
        
        # Extracting the image for the video
        png = 'simulation_{}_{:04d}.png'.format(include_channels[0],number)
        png_file = os.path.join(images_dir,png)
        plotter.window_size = [3000, 3000]
        plotter.screenshot(png_file)
        print(png_file)
        

        # Merging sismograms image and plotter
        background = cv2.imread(png_file)
        foreground = cv2.imread(plot_file, -1)  # Load with alpha channel

        if background is None:
            print(f"Error: Cannot load the background image from {png_file}")
        if foreground is None:
            print(f"Error: Cannot load the foreground image from {plot_file}")
        
        if background is not None and foreground is not None:
            # Resize the foreground image if needed (for example, 200x200 pixels)
            foreground = cv2.resize(foreground, (1100, 500))
        
            # Get dimensions of the foreground
            h, w = foreground.shape[:2]
        
            # Define the position for overlay (top-right corner)
            y, x = 0, background.shape[1] - w
        
            # Ensure the foreground image's alpha channel is used for blending
            for c in range(0, 3):
                background[y:y+h, x:x+w, c] = (foreground[:, :, c] * (foreground[:, :, 3] / 255.0) + 
                                               background[y:y+h, x:x+w, c] * (1.0 - foreground[:, :, 3] / 255.0))
        
            # Save the result image with the same size as the background
            result_image = 'result_image_{}_{:04d}.png'.format(include_channels[0],number)
            result_path = os.path.join(result_dir, result_image)
            cv2.imwrite(result_path, background)
        
    # plotter.close()
    gc.collect()

    print(number)
    return number

if __name__ == '__main__':
    number = animate_pyvista_png(top_dir=top_dir,
run=run,element_name=element_name,short_title=short_title,station_group=station_group,run_title=run_title,new_station_group=new_station_group)