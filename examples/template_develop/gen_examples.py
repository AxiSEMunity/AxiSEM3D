#
#  vicinity.py
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 6/20/20.
#  Copyright © 2020 Kuangdai Leng. All rights reserved.
#

#  regenerate the examples when any input parameters change

import os
import numpy as np


###################### tools ######################
# read a file
def read(fname):
    fin = open(fname, 'r')
    return fin.read()
    
# replace in a string
def replace_in_string(contents, from_strs, to_strs):
    for istr, from_str in enumerate(from_strs):
        contents = contents.replace(from_str, to_strs[istr])
    return contents

# replace in a file
def replace_in_file(fname, from_strs, to_strs):
    contents = read(fname)
    contents = replace_in_string(contents, from_strs, to_strs)
    fout = open(fname, 'w')
    fout.write(contents)
    fout.close()


###################### item templates ######################
item_source_VIR = read('input/item_source_VIR.yaml')[:-1]
item_stations_GSN = read('input/item_stations_GSN.yaml')[:-1]
item_elements_mantle = read('input/item_elements_mantle.yaml')[:-1]


################ 01_spherical_Earth_PREM_100s ################
# copy
ex_name = '01_spherical_Earth_PREM_50s'
input_dir = '../%s/input' % ex_name
os.system('cp input/inparam.*.yaml %s' % input_dir)

# source
replace_in_file(input_dir + '/inparam.source.yaml',
                ['list_of_sources: []'],
                ['list_of_sources:\n' + item_source_VIR])

# output
item_stations_US = replace_in_string(item_stations_GSN,
['global_seismic_network_GSN', 'GSN.txt', 'channels: [U]',
 'format: ASCII_STATION', 'sampling_period: DT'],
['USArray_transportable', 'US_TA.txt', 'channels: [U3, E_I1, R3]',
 'format: ASCII_CHANNEL', 'sampling_period: 5.'])
replace_in_file(input_dir + '/inparam.output.yaml',
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' +
                 item_stations_GSN + '\n' + item_stations_US])


################ 03_Cartesian_SEG_EAGE_salt_5Hz ################
# copy
ex_name = '03_Cartesian_SEG_EAGE_salt_5Hz'
input_dir = '../%s/input' % ex_name
os.system('cp input/inparam.*.yaml %s' % input_dir)

# model
replace_in_file(input_dir + '/inparam.model.yaml',
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'attenuation: CG4'],
                ['exodus_mesh: local_mesh__SEG_salt__5Hz.e',
                 'attenuation: NONE'])
                
# nr
replace_in_file(input_dir + '/inparam.nr.yaml',
                ['constant: 5'],
                ['constant: 1'])
                
# source
item_source_mono = replace_in_string(item_source_VIR,
['- VIRGINIA_201108231751A:',
 'latitude_longitude: [37.91, -77.93]', 'depth: 12e3',
 'ellipticity: true', 'depth_below_solid_surface: true',
 'undulated_geometry: true',
 'data: [4.71e24, 3.81e22, -4.74e24, 3.99e23, -8.05e23, -1.23e24]',
 'half_duraction: 50.', 'use_derivative_integral: ERF'],
['- the_only_source:',
 'latitude_longitude: [90., 0.]', 'depth: 1000.',
 'ellipticity: false', 'depth_below_solid_surface: false',
 'undulated_geometry: false',
 'data: [1e20, 0., 0., 0., 0., 0.]',
 'half_duraction: 0.2', 'use_derivative_integral: GAUSSIAN'])
replace_in_file(input_dir + '/inparam.source.yaml',
                ['list_of_sources: []'],
                ['list_of_sources:\n' + item_source_mono])
replace_in_file(input_dir + '/inparam.source.yaml',
                ['record_length: 1800.', 'Courant_number: 0.6'],
                ['record_length: 10.', 'Courant_number: 0.3'])
                
# output
item_elements_inplane = replace_in_string(item_elements_mantle,
['- Fourier_coefficients_spherical_Earth_whole_mantle:',
 'horizontal_range: [0, 3.15]', 'vertical_range: [3480e3, 6371e3]',
 'GLL_points_one_edge: FULL',
 'phi_list: []',
 'sampling_period: DT'],
['- orthogonal_azimuthal_slices:',
 'horizontal_range: [0, 1e10]', 'vertical_range: [0, 1e10]',
 'GLL_points_one_edge: [0, 2, 4]',
 'phi_list: [0, 1.57079632679, 3.14159265359, 4.71238898038]',
 'sampling_period: 0.05'])

item_elements_ocean_floor = replace_in_string(item_elements_mantle,
['- Fourier_coefficients_spherical_Earth_whole_mantle:',
 'horizontal_range: [0, 3.15]', 'vertical_range: [3480e3, 6371e3]',
 'edge_dimension: BOTH', 'edge_position: 6371e3',
 'sampling_period: DT'],
['- Fourier_coefficients_ocean_floor:',
 'horizontal_range: [0, 1e10]', 'vertical_range: [0, 1e10]',
 'edge_dimension: VERTICAL', 'edge_position: 6370.8e3',
 'sampling_period: 0.01'])

replace_in_file(input_dir + '/inparam.output.yaml',
                ['list_of_element_groups: []'],
                ['list_of_element_groups:\n' +
                 item_elements_inplane + '\n' + item_elements_ocean_floor])

# advanced
replace_in_file(input_dir + '/inparam.advanced.yaml',
                ['nproc_per_group: 1'],
                ['nproc_per_group: 24'])

# copy 3D
input3D_dir = input_dir + '3D'
os.system('cp %s/inparam.*.yaml %s' % (input_dir, input3D_dir))

# model
item_3D_model_list = read(input3D_dir + '/list_of_3D_models.yaml')[:-1]
replace_in_file(input3D_dir + '/inparam.model.yaml',
                ['list_of_3D_models: []'],
                ['list_of_3D_models:\n' + item_3D_model_list])

# nr
replace_in_file(input3D_dir + '/inparam.nr.yaml',
                ['constant: 1'],
                ['constant: 50'])
                
# output
replace_in_file(input3D_dir + '/inparam.output.yaml',
                ['buffer_size: 1000'],
                ['buffer_size: 100'])

# advanced
replace_in_file(input3D_dir + '/inparam.advanced.yaml',
                ['loop_info_interval: 1000'],
                ['loop_info_interval: 100'])
