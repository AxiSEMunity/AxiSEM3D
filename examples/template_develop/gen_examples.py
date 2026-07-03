#!/usr/bin/env python3
#
#  gen_examples.py
#  AxiSEM3D
#
#  Originally created by Kuangdai Leng on 6/20/20.
#  Copyright © 2020 Kuangdai Leng. All rights reserved.
#

#  Heavily edited by Jonathan Wolf on 2026-04-13
#
#  Regenerate ALL examples from base templates + snippet files.
#  Safe to run from any directory.

import os
import glob
import shutil

# Always operate relative to this script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

EXAMPLES_DIR = os.path.dirname(SCRIPT_DIR)


###################### tools ######################
def read(fname):
    with open(fname, 'r') as f:
        return f.read()


def snippet(fname):
    """Read an input snippet, stripping the trailing newline."""
    return read('input/' + fname).rstrip('\n')


def replace_in_string(contents, from_strs, to_strs):
    if len(from_strs) != len(to_strs):
        raise ValueError('replace_in_string: length mismatch (%d vs %d)'
                         % (len(from_strs), len(to_strs)))
    for i, token in enumerate(from_strs):
        if token not in contents:
            raise ValueError('replace_in_string: token not found:\n  %s'
                             % token)
        contents = contents.replace(token, to_strs[i])
    return contents


def replace_in_file(fname, from_strs, to_strs):
    contents = read(fname)
    contents = replace_in_string(contents, from_strs, to_strs)
    with open(fname, 'w') as f:
        f.write(contents)


def copy_templates(dst_dir, which=None):
    """Copy base inparam template YAML files to dst_dir.
    If *which* is given, copy only those filenames (e.g. ['model'])."""
    os.makedirs(dst_dir, exist_ok=True)
    if which is None:
        for src in sorted(glob.glob('input/inparam.*.yaml')):
            shutil.copy2(src, dst_dir)
    else:
        for name in which:
            shutil.copy2('input/inparam.%s.yaml' % name, dst_dir)


def copy_input_set(src_dir, dst_dir):
    os.makedirs(dst_dir, exist_ok=True)
    for src in sorted(glob.glob(os.path.join(src_dir, 'inparam.*.yaml'))):
        shutil.copy2(src, dst_dir)


def ex_path(*parts):
    return os.path.join(EXAMPLES_DIR, *parts)


def jp(base, fname):
    return os.path.join(base, fname)


###################### item templates ######################
item_source_VIR = snippet('item_source_VIR.yaml')
item_stations_GSN = snippet('item_stations_GSN.yaml')
item_elements_mantle = snippet('item_elements_mantle.yaml')

# Derived: on-axis monopole source (used by ex03)
item_source_onaxis = replace_in_string(item_source_VIR,
    ['- VIRGINIA_201108231751A:',
     'latitude_longitude: [37.91, -77.93]', 'depth: 12e3',
     'ellipticity: true', 'depth_below_solid_surface: true',
     'undulated_geometry: true',
     'data: [4.71e24, 3.81e22, -4.74e24, 3.99e23, -8.05e23, -1.23e24]',
     'half_duration: 50.', 'use_derivative_integral: ERF'],
    ['- the_only_source:',
     'latitude_longitude: [90., 0.]', 'depth: 8000.',
     'ellipticity: false', 'depth_below_solid_surface: false',
     'undulated_geometry: false',
     'data: [1e20, 0., 0., 0., 0., 0.]',
     'half_duration: 0.1', 'use_derivative_integral: GAUSSIAN'])

# Derived: ON_AXIS source variant used by ex04
item_source_onaxis_ex04 = replace_in_string(
    item_source_onaxis,
    ['latitude_longitude: [90., 0.]'],
    ['latitude_longitude: ON_AXIS'])

# Derived: USArray station group (used by ex00, ex01)
item_stations_US = replace_in_string(item_stations_GSN,
    ['global_seismic_network_GSN', 'GSN.txt', 'channels: [U]',
     'format: ASCII_STATION', 'sampling_period: DT'],
    ['USArray_transportable', 'US_TA.txt', 'channels: [U3, E_I1, R3]',
     'format: ASCII_CHANNEL', 'sampling_period: 5.'])

# Derived: USArray station group for ex02, which keeps scalar ASCII stations
item_stations_US_ex02 = replace_in_string(item_stations_GSN,
    ['global_seismic_network_GSN', 'GSN.txt',
     'format: ASCII_STATION', 'sampling_period: DT'],
    ['USArray_transportable', 'US_TA.txt',
     'format: ASCII_CHANNEL', 'sampling_period: 5.'])

# Derived: inplane-slice element group (used by ex03, ex04a, ex04c, ex06)
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


################ 00_global_1D ################
print('Generating 00_global_1D ...')
d = ex_path('00_global_1D', 'input')
copy_templates(d)

replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []',
                 'Courant_number: 0.6'],
                ['list_of_sources:\n' + item_source_VIR,
                 'Courant_number: 1.0'])

replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' +
                 item_stations_GSN + '\n' + item_stations_US])


################ 01_S362ANI_EMC_global ################
print('Generating 01_S362ANI_EMC_global ...')
d = ex_path('01_S362ANI_EMC_global', 'input')
copy_templates(d)

replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []',
                 'Courant_number: 0.6'],
                ['list_of_sources:\n' + item_source_VIR,
                 'Courant_number: 1.0'])

replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' +
                 item_stations_GSN + '\n' + item_stations_US])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['list_of_3D_models: []'],
                ['list_of_3D_models:\n' + snippet('ex01_list_of_3D_models.yaml')])

replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['constant: 5'],
                ['constant: 50'])


################ 02_3d_crust_S362ANI_regional ################
print('Generating 02_3d_crust_S362ANI_regional ...')

# input_share: source, nr, output, advanced (no model)
d_share = ex_path('02_3d_crust_S362ANI_regional', 'input_share')
copy_templates(d_share, which=['source', 'nr', 'output', 'advanced'])

replace_in_file(jp(d_share, 'inparam.source.yaml'),
                ['list_of_sources: []'],
                ['list_of_sources:\n' + item_source_VIR])

replace_in_file(jp(d_share, 'inparam.nr.yaml'),
                ['constant: 5'],
                ['constant: 100'])

replace_in_file(jp(d_share, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' + item_stations_US_ex02])

# input_with_1d_crust: model only (regional mesh + S362ANI)
d_1d = ex_path('02_3d_crust_S362ANI_regional', 'input_with_1d_crust')
copy_templates(d_1d, which=['model'])

replace_in_file(jp(d_1d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e'],
                ['exodus_mesh: regional_mesh__prem_ani__50s.e'])

replace_in_file(jp(d_1d, 'inparam.model.yaml'),
                ['list_of_3D_models: []'],
                ['list_of_3D_models:\n' + snippet('ex01_list_of_3D_models.yaml')])

# input_with_3d_crust: model only (regional mesh + S362ANI + crust1.0)
d_3d = ex_path('02_3d_crust_S362ANI_regional', 'input_with_3d_crust')
copy_templates(d_3d, which=['model'])

replace_in_file(jp(d_3d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e'],
                ['exodus_mesh: regional_mesh__prem_ani__50s.e'])

replace_in_file(jp(d_3d, 'inparam.model.yaml'),
                ['list_of_3D_models: []'],
                ['list_of_3D_models:\n' +
                 snippet('ex02_list_of_3D_models_3d_crust.yaml')])


################ 03_salt_body_SEG_local ################
print('Generating 03_salt_body_SEG_local ...')

# --- input1D ---
d1 = ex_path('03_salt_body_SEG_local', 'input1D')
copy_templates(d1)

replace_in_file(jp(d1, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'attenuation: CG4'],
                ['exodus_mesh: local_mesh__SEG_salt__5Hz.e',
                 'attenuation: NONE'])

replace_in_file(jp(d1, 'inparam.nr.yaml'),
                ['constant: 5'],
                ['constant: 1'])

item_source_mono = replace_in_string(item_source_onaxis,
    ['depth: 8000.', 'half_duration: 0.1'],
    ['depth: 1000.', 'half_duration: 0.2'])
replace_in_file(jp(d1, 'inparam.source.yaml'),
                ['list_of_sources: []'],
                ['list_of_sources:\n' + item_source_mono])
replace_in_file(jp(d1, 'inparam.source.yaml'),
                ['record_length: 1800.', 'Courant_number: 0.6'],
                ['record_length: 10.', 'Courant_number: 0.5'])

item_elements_ocean_floor = replace_in_string(item_elements_mantle,
    ['- Fourier_coefficients_spherical_Earth_whole_mantle:',
     'horizontal_range: [0, 3.15]', 'vertical_range: [3480e3, 6371e3]',
     'edge_dimension: BOTH', 'edge_position: 6371e3',
     'sampling_period: DT'],
    ['- Fourier_coefficients_ocean_floor:',
     'horizontal_range: [0, 1e10]', 'vertical_range: [0, 1e10]',
     'edge_dimension: VERTICAL', 'edge_position: 6370.8e3',
     'sampling_period: 0.01'])

replace_in_file(jp(d1, 'inparam.output.yaml'),
                ['list_of_element_groups: []'],
                ['list_of_element_groups:\n' +
                 item_elements_inplane + '\n' + item_elements_ocean_floor])

replace_in_file(jp(d1, 'inparam.advanced.yaml'),
                ['nproc_per_group: 1'],
                ['nproc_per_group: 24'])

# --- input3D: copy from 1D, then customize ---
d3 = ex_path('03_salt_body_SEG_local', 'input3D')
copy_input_set(d1, d3)

replace_in_file(jp(d3, 'inparam.model.yaml'),
                ['list_of_3D_models: []'],
                ['list_of_3D_models:\n' + snippet('ex03_list_of_3D_models.yaml')])

replace_in_file(jp(d3, 'inparam.nr.yaml'),
                ['constant: 1'],
                ['constant: 50'])

replace_in_file(jp(d3, 'inparam.output.yaml'),
                ['buffer_size: 1000'],
                ['buffer_size: 100'])

replace_in_file(jp(d3, 'inparam.advanced.yaml'),
                ['loop_info_interval: 1000'],
                ['loop_info_interval: 100'])


################ 04_simple_3d_shapes ################
print('Generating 04_simple_3d_shapes ...')

# --- 04a: example_input_cartesian ---
d = ex_path('04_simple_3d_shapes', 'example_input_cartesian')
copy_templates(d)

replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []', 'record_length: 1800.'],
                ['list_of_sources:\n' + item_source_onaxis_ex04,
                 'record_length: 15'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'lat_lon_north_pole_mesh: SOURCE',
                 'boundaries: [RIGHT, BOTTOM]',
                 'attenuation: CG4',
                 'list_of_3D_models: []'],
                ['exodus_mesh: homogenous_cartesian.e',
                 'lat_lon_north_pole_mesh: [90., 0.]',
                 'boundaries: []',
                 'attenuation: NONE',
                 'list_of_3D_models:\n' + snippet('ex04a_list_of_3D_models.yaml')])

item_stn_04a = replace_in_string(item_stations_GSN,
    ['global_seismic_network_GSN', 'GSN.txt',
     'horizontal_x1_x2: LATITUDE_LONGITUDE',
     'ellipticity: true', 'depth_below_solid_surface: true',
     'undulated_geometry: true', 'flush: true'],
    ['station_1', 'stn.txt',
     'horizontal_x1_x2: DISTANCE_AZIMUTH',
     'ellipticity: false', 'depth_below_solid_surface: false',
     'undulated_geometry: false', 'flush: false'])

item_elem_04a = replace_in_string(item_elements_inplane,
    ['sampling_period: 0.05', 'buffer_size: 1000'],
    ['sampling_period: 0.1', 'buffer_size: 100'])

replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []',
                 'list_of_element_groups: []'],
                ['list_of_station_groups:\n' + item_stn_04a,
                 'list_of_element_groups:\n' + item_elem_04a])

replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['constant: 5'],
                ['constant: 1'])

replace_in_file(jp(d, 'inparam.advanced.yaml'),
                ['loop_info_interval: 1000', 'nproc_per_group: 1'],
                ['loop_info_interval: 100', 'nproc_per_group: 24'])

# --- 04b: example_release_paper ---
d = ex_path('04_simple_3d_shapes', 'example_release_paper', 'input')
copy_templates(d)

item_source_04b = replace_in_string(item_source_onaxis_ex04,
    ['depth: 8000.',
     'data: [1e20, 0., 0., 0., 0., 0.]',
     'unit: 1e-7',
     'half_duration: 0.1'],
    ['depth: 7500.',
     'data: [1.73e15, -2.81e14, -1.45e15, 2.12e15, 4.55e15, -6.57e15]',
     'unit: 1',
     'half_duration: 0.5'])
replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []',
                 'record_length: 1800.',
                 'Courant_number: 0.6'],
                ['list_of_sources:\n' + item_source_04b,
                 'record_length: 50',
                 'Courant_number: 0.4'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'lat_lon_north_pole_mesh: SOURCE',
                 'attenuation: CG4',
                 'list_of_3D_models: []'],
                ['exodus_mesh: paper_example.e',
                 'lat_lon_north_pole_mesh: [90., 0.]',
                 'attenuation: FULL',
                 'list_of_3D_models:\n' + snippet('ex04b_list_of_3D_models.yaml')])

item_stn_04b = replace_in_string(item_stations_GSN,
    ['global_seismic_network_GSN', 'GSN.txt',
     'horizontal_x1_x2: LATITUDE_LONGITUDE',
     'ellipticity: true', 'depth_below_solid_surface: true',
     'undulated_geometry: true', 'flush: true'],
    ['station_1', 'paperstns.txt',
     'horizontal_x1_x2: DISTANCE_AZIMUTH',
     'ellipticity: false', 'depth_below_solid_surface: false',
     'undulated_geometry: false', 'flush: false'])
replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' + item_stn_04b])

replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['constant: 5'],
                ['constant: 360'])

replace_in_file(jp(d, 'inparam.advanced.yaml'),
                ['loop_info_interval: 1000', 'nproc_per_group: 1'],
                ['loop_info_interval: 100', 'nproc_per_group: 40'])

# --- 04c: example_single_plume ---
d = ex_path('04_simple_3d_shapes', 'example_single_plume')
copy_templates(d)

replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []',
                 'record_length: 1800.',
                 'Courant_number: 0.6'],
                ['list_of_sources:\n' + item_source_onaxis_ex04,
                 'record_length: 1200',
                 'Courant_number: 0.5'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'lat_lon_north_pole_mesh: SOURCE',
                 'boundaries: [RIGHT, BOTTOM]',
                 'attenuation: CG4',
                 'list_of_3D_models: []'],
                ['exodus_mesh: homogenous_global.e',
                 'lat_lon_north_pole_mesh: [90., 0.]',
                 'boundaries: []',
                 'attenuation: NONE',
                 'list_of_3D_models:\n' + snippet('ex04c_list_of_3D_models.yaml')])

item_elem_04c = replace_in_string(item_elements_inplane,
    ['sampling_period: 0.05', 'buffer_size: 1000'],
    ['sampling_period: 5.0', 'buffer_size: 100'])
replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_element_groups: []'],
                ['list_of_element_groups:\n' + item_elem_04c])

replace_in_file(jp(d, 'inparam.advanced.yaml'),
                ['loop_info_interval: 1000', 'nproc_per_group: 1'],
                ['loop_info_interval: 10', 'nproc_per_group: 24'])


################ 05_anisotropy_global ################
print('Generating 05_anisotropy_global ...')

# --- 05a: sim_US32_olivineE ---
d = ex_path('05_anisotropy_global',
            '2012-07-03_paper_example_50s', 'sim_US32_olivineE', 'input')
copy_templates(d)

item_src_05a = snippet('ex05a_list_of_sources.yaml')
replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []', 'record_length: 1800.'],
                ['list_of_sources:\n' + item_src_05a, 'record_length: 2000.'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'attenuation: CG4',
                 'list_of_3D_models: []'],
                ['exodus_mesh: AxiSEM_prem_iso_50.e',
                 'attenuation: FULL',
                 'list_of_3D_models:\n' + snippet('ex05a_list_of_3D_models.yaml')])
replace_in_file(jp(d, 'inparam.model.yaml'),
                ['enable_Clayton_Enquist: true'],
                ['enable_Clayton_Enquist: false'])

item_stn_05a = replace_in_string(item_stations_GSN,
    ['global_seismic_network_GSN', 'GSN.txt', 'time_window: FULL'],
    ['JW_stations', 'eq3.txt', 'time_window: [0, 2000]'])
replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' + item_stn_05a])

replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['type_Nr: CONSTANT', 'constant: 5',
                 'nc_data_file: pointwise.nc'],
                ['type_Nr: POINTWISE', 'constant: 300',
                 'nc_data_file: scanning_output_Nr.nc'])
replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['Nr_at_control_depths: [100, 100, 50, 50]',
                 'nc_data_file: structured.nc',
                 'enable_scanning: false'],
                ['Nr_at_control_depths: [100, 100, 100, 100]',
                 'nc_data_file: bla.nc',
                 'enable_scanning: true'])

# --- 05b: sim_US32_olivineE_fullNU (same as 05a but CONSTANT Nr) ---
d = ex_path('05_anisotropy_global',
            '2012-07-03_paper_example_50s', 'sim_US32_olivineE_fullNU', 'input')
copy_templates(d)

replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []', 'record_length: 1800.'],
                ['list_of_sources:\n' + item_src_05a, 'record_length: 2000.'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'attenuation: CG4',
                 'list_of_3D_models: []'],
                ['exodus_mesh: AxiSEM_prem_iso_50.e',
                 'attenuation: FULL',
                 'list_of_3D_models:\n' + snippet('ex05a_list_of_3D_models.yaml')])
replace_in_file(jp(d, 'inparam.model.yaml'),
                ['enable_Clayton_Enquist: true'],
                ['enable_Clayton_Enquist: false'])

replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' + item_stn_05a])

replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['constant: 5'],
                ['constant: 300'])
replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['Nr_at_control_depths: [100, 100, 50, 50]',
                 'enable_scanning: false'],
                ['Nr_at_control_depths: [100, 100, 100, 100]',
                 'enable_scanning: true'])

# --- 05c: deep_mantle_anisotropy ---
d = ex_path('05_anisotropy_global',
            'deep_mantle_anisotropy_full_Cij_50s',
            'sim_lowermost_mantle_ani', 'input')
copy_templates(d)

item_src_05c = snippet('ex05c_list_of_sources.yaml')
replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []', 'record_length: 1800.'],
                ['list_of_sources:\n' + item_src_05c, 'record_length: 2000.'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'flattening_on_surface: WGS84',
                 'attenuation: CG4',
                 'list_of_3D_models: []'],
                ['exodus_mesh: AxiSEM_prem_iso_50.e',
                 'flattening_on_surface: SPHERE',
                 'attenuation: FULL',
                 'list_of_3D_models:\n' + snippet('ex05c_list_of_3D_models.yaml')])
replace_in_file(jp(d, 'inparam.model.yaml'),
                ['enable_Clayton_Enquist: true'],
                ['enable_Clayton_Enquist: false'])

item_stn_05c = replace_in_string(item_stations_GSN,
    ['global_seismic_network_GSN', 'GSN.txt',
     'coordinate_frame: RTZ', 'channels: [U]', 'time_window: FULL'],
    ['JW_stations', 'stations.txt',
     'coordinate_frame: ENZ', 'channels: [U, G, E, S, R]',
     'time_window: [0, 3600]'])
replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' + item_stn_05c])

# --- 05d: sim1_ani_prem_mesh (no 3D models) ---
d = ex_path('05_anisotropy_global',
            'PREM_anisotropy_w_and_wo_full_Cij_50s',
            'sim1_ani_prem_mesh', 'input')
copy_templates(d)

replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []', 'record_length: 1800.'],
                ['list_of_sources:\n' + item_src_05c, 'record_length: 2000.'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'flattening_on_surface: WGS84',
                 'attenuation: CG4'],
                ['exodus_mesh: AxiSEM_prem_ani_50.e',
                 'flattening_on_surface: SPHERE',
                 'attenuation: FULL'])
replace_in_file(jp(d, 'inparam.model.yaml'),
                ['enable_Clayton_Enquist: true'],
                ['enable_Clayton_Enquist: false'])

replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' + item_stn_05c])

# --- 05e: sim2_iso_prem_mesh_plus_ani ---
d = ex_path('05_anisotropy_global',
            'PREM_anisotropy_w_and_wo_full_Cij_50s',
            'sim2_iso_prem_mesh_plus_ani', 'input')
copy_templates(d)

replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []', 'record_length: 1800.'],
                ['list_of_sources:\n' + item_src_05c, 'record_length: 2000.'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'flattening_on_surface: WGS84',
                 'attenuation: CG4',
                 'list_of_3D_models: []'],
                ['exodus_mesh: AxiSEM_prem_iso_50.e',
                 'flattening_on_surface: SPHERE',
                 'attenuation: FULL',
                 'list_of_3D_models:\n' + snippet('ex05e_list_of_3D_models.yaml')])
replace_in_file(jp(d, 'inparam.model.yaml'),
                ['enable_Clayton_Enquist: true'],
                ['enable_Clayton_Enquist: false'])

replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' + item_stn_05c])


################ 06_finite_source_local ################
print('Generating 06_finite_source_local ...')
d = ex_path('06_finite_source_local', 'input')
copy_templates(d)

# source: finite-fault sources are too large for inline generation;
# they come from input_setup.ipynb.  We insert the pre-extracted snippet.
item_src_06 = snippet('ex06_list_of_sources.yaml')
replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []',
                 '# note: the start time depends on the source-time functions\n'
                 '    record_length: 1800.',
                 '# note: 1) Δt increases with the Courant number; decrease it when\n'
                 '    #          numerical instability occurs\n'
                 '    #       2) [safe] 0.5 <===> 1.0 [aggressive]; 0.6~0.7 normally works\n'
                 '    #       3) if Courant_number < 0.3 but instability still occurs,\n'
                 '    #          it is likely to be an issue caused by an input 3D model\n'
                 '    #          (e.g., mislocation near a model boundary)\n'
                 '    Courant_number: 0.6'],
                ['list_of_sources:\n' + item_src_06,
                 '# note: the start time depends on the source-time functions\n'
                 '    # 30 s (was 40): the figure\'s snapshot is at 60% of the record and the\n'
                 '    # near-field SFBA shaking/PGV for this Mw~4.4 source peaks within ~20 s, so\n'
                 '    # 30 s is ample and trims ~25% off the (mesh-Δt-bound) ~26 h runtime.\n'
                 '    record_length: 30.',
                 '# note: Δt increases with this number\n'
                 '    #       0.6 was marginally unstable on the fine 0.5 s mesh (late-onset\n'
                 '    #       exponential blow-up ~t=8-25 s); 0.5 gives a safe stability margin.\n'
                 '    Courant_number: 0.5'])

replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'lat_lon_north_pole_mesh: SOURCE',
                 'list_of_3D_models: []'],
                ['exodus_mesh: AxiSEMCartesian_sfba_2018_2s.e',
                 'lat_lon_north_pole_mesh: [37.7, -122.1]',
                 'list_of_3D_models:\n' + snippet('ex06_list_of_3D_models.yaml')])
replace_in_file(jp(d, 'inparam.model.yaml'),
                ['relative_spans: [.05, .05]',
                 'gamma_expr_solid: 1.1 / T0 * (VS / VP)^2 * exp(-0.04 * SPAN / (VP * T0))',
                 'gamma_expr_fluid: 0.88 / T0 * exp(-0.04 * SPAN / (VP * T0))'],
                ['relative_spans: [.1, .1]',
                 'gamma_expr_solid: 4.4 / T0 * (VS / VP)^2 * exp(-0.08 * SPAN / (VP * T0))',
                 'gamma_expr_fluid: 1.76 / T0 * exp(-0.08 * SPAN / (VP * T0))'])

# output: complex station + element groups from snippets
item_stn_06 = snippet('ex06_list_of_station_groups.yaml')
item_elem_06 = snippet('ex06_list_of_element_groups.yaml')
replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []',
                 'list_of_element_groups: []'],
                ['list_of_station_groups:\n' + item_stn_06,
                 'list_of_element_groups:\n' + item_elem_06])

# nr: switch to POINTWISE
replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['type_Nr: CONSTANT', 'constant: 5',
                 'nc_data_file: pointwise.nc'],
                ['type_Nr: POINTWISE', 'constant: 3000',
                 'nc_data_file: point_source_approx.nc'])
replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['enable_scanning: false',
                 'output_file: scanning_output_Nr.nc',
                 'relative_amplitude_skipped: 0.1',
                 'absolute_amplitude_skipped: 1e-12',
                 'max_num_peaks: 10000',
                 'vertex_only: true'],
                ['enable_scanning: true',
                 'output_file: SFBA_finite_rupture_Nr_learned.nc',
                 'relative_amplitude_skipped: 0.',
                 'absolute_amplitude_skipped: 1e-14',
                 'max_num_peaks: 10',
                 'vertex_only: false'])

# advanced
replace_in_file(jp(d, 'inparam.advanced.yaml'),
                ['level: ESSENTIAL', 'loop_info_interval: 1000', 'nproc_per_group: 1'],
                ['level: DETAILED', 'loop_info_interval: 100', 'nproc_per_group: 8'])


################ 08_atmosphere_Mars_global ################
print('Generating 08_atmosphere_Mars_global ...')
d = ex_path('08_atmosphere_Mars_global', 'axisem3d_mars_atm', 'input')
copy_templates(d)

# model: Mars atmospheric mesh, absorbing TOP only, single sponge layer
replace_in_file(jp(d, 'inparam.model.yaml'),
                ['exodus_mesh: global_mesh__prem_ani__50s.e',
                 'boundaries: [RIGHT, BOTTOM]',
                 'relative_spans: [.05, .05]'],
                ['exodus_mesh: tayak_60km_10s.e',
                 'boundaries: [TOP]',
                 'relative_spans: [.01]'])

# source: FLUID_PRESSURE at Mars surface
item_source_08 = replace_in_string(item_source_VIR,
    ['- VIRGINIA_201108231751A:',
     'latitude_longitude: [37.91, -77.93]',
     'depth: 12e3',
     'ellipticity: true',
     'type: MOMENT_TENSOR',
     'data: [4.71e24, 3.81e22, -4.74e24, 3.99e23, -8.05e23, -1.23e24]',
     'half_duration: 50.',
     'use_derivative_integral: ERF'],
    ['- s1094b:',
     'latitude_longitude: [0,0]',
     'radius: 3394000',
     'ellipticity: false',
     'type: FLUID_PRESSURE',
     'data: [1e1]',
     'half_duration: 0.0',
     'use_derivative_integral: GAUSSIAN'])
replace_in_file(jp(d, 'inparam.source.yaml'),
                ['list_of_sources: []',
                 'Courant_number: 0.6'],
                ['list_of_sources:\n' + item_source_08,
                 'Courant_number: 0.5'])

# output: MARS station group
item_stn_08 = replace_in_string(item_stations_GSN,
    ['global_seismic_network_GSN', 'GSN.txt',
     'ellipticity: true'],
    ['MARS', 'mars.txt',
     'ellipticity: false'])
replace_in_file(jp(d, 'inparam.output.yaml'),
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' + item_stn_08])

# nr: 1 for FLUID_PRESSURE (monopole)
replace_in_file(jp(d, 'inparam.nr.yaml'),
                ['constant: 5'],
                ['constant: 1'])


################ done ################
print('')
print('Done — all examples regenerated (00, 01, 02, 03, 04, 05, 06, 08).')