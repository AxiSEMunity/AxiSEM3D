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
item_source_VIR = read('input/item_source.yaml')[:-1]
item_stations_GSN = read('input/item_stations.yaml')[:-1]


################ 01_spherical_Earth_PREM_100s ################
# copy
ex_name = '01_spherical_Earth_PREM_50s'
input_dir = '../%s/input' % ex_name
os.system('cp input/inparam.*.yaml %s' % input_dir)

# source
replace_in_file(input_dir + '/inparam.source.yaml',
                ['list_of_sources: []'],
                ['list_of_sources:\n' + item_source_VIR])

# stations
item_stations_US = replace_in_string(item_stations_GSN,
['global_seismic_network_GSN', 'GSN.txt', 'channels: [U]',
 'format: ASCII_STATION', 'sampling_period: DT'],
['USArray_transportable', 'US_TA.txt', 'channels: [U3, E_I1, R3]',
 'format: ASCII_CHANNEL', 'sampling_period: 5.'])
replace_in_file(input_dir + '/inparam.output.yaml',
                ['list_of_station_groups: []'],
                ['list_of_station_groups:\n' +
                 item_stations_GSN + '\n' + item_stations_US])
