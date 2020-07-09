#
#  vicinity.py
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 3/19/20.
#  Copyright © 2020 Kuangdai Leng. All rights reserved.
#

#  print c++ code for constants used in vicinity.hpp

import numpy as np

# node to ipnt
def node_ipnt(npol, node):
    if node == 0:
        ipol = 0
        jpol = 0
    elif node == 1:
        ipol = npol
        jpol = 0
    elif node == 2:
        ipol = npol
        jpol = npol
    elif node == 3:
        ipol = 0
        jpol = npol
    return ipol * (npol + 1) + jpol
    
# edge to ipnt
def edge_ipnt(npol, edge):
    ipnts = []
    if edge == 0:
        for i in np.arange(npol + 1):
            ipol = i
            jpol = 0
            ipnts.append(ipol * (npol + 1) + jpol)
    elif edge == 1:
        for i in np.arange(npol + 1):
            ipol = npol
            jpol = i
            ipnts.append(ipol * (npol + 1) + jpol)
    elif edge == 2:
        for i in np.arange(npol + 1):
            ipol = npol - i
            jpol = npol
            ipnts.append(ipol * (npol + 1) + jpol)
    elif edge == 3:
        for i in np.arange(npol + 1):
            ipol = 0
            jpol = npol - i
            ipnts.append(ipol * (npol + 1) + jpol)
    return np.array(ipnts)
    
# print std::vector<int>
def print_array4(key, array, after_prefix=''):
    fmt={'int': lambda x: '%d' % x}
    print("const std::vector<int> %s = %s;" % (key,
    np.array2string(array, separator=', ', max_line_width=64,
    formatter=fmt).replace('[', '{' + after_prefix).replace(']', '}')))
    
# print std::vector<std::vector<int>>
def print_array4_vector(key, array):
    fmt={'int': lambda x: '%d' % x}
    print("const std::vector<std::vector<int>> %s = %s;" % (key,
    np.array2string(array, separator=', ', max_line_width=64,
    formatter=fmt).replace('[', '{').replace(']', '}').replace('{{', '{\n{')))

# print npol
def print_vicinity_const(npol):
    # node
    node_ipnts = []
    for node in np.arange(4):
        node_ipnts.append(node_ipnt(npol, node))
    node_ipnts = np.array(node_ipnts)
        
    # edge
    edge_ipnts = []
    for edge in np.arange(4):
        edge_ipnts.append(edge_ipnt(npol, edge))
    edge_ipnts = np.array(edge_ipnts)
    
    # print
    if npol == 1:
        print("#if _NPOL == %d" % npol)
    else:
        print("#elif _NPOL == %d" % npol)
    print("// nPol = %d" % npol)
    print_array4('gNodeIPnt', node_ipnts)
    print_array4_vector('gEdgeIPnt', edge_ipnts)
    print_array4('gEdgeIPntAll', np.unique(edge_ipnts), after_prefix='\n')
    if npol == 8:
        print("#endif")
    
# print all npol from 1 to 8
print_vicinity_const(1)
print_vicinity_const(2)
print_vicinity_const(3)
print_vicinity_const(4)
print_vicinity_const(5)
print_vicinity_const(6)
print_vicinity_const(7)
print_vicinity_const(8)
