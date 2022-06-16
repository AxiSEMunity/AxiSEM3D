#!bin/bash
# to generate mesh
python -m salvus_mesh_lite.interface AxiSEMCartesian --basic.model homogenous_cart.bm --basic.period .2 --cartesian2Daxisem.x 25. --cartesian2Daxisem.min_z 0.0 --attenuation.frequencies 0.01 10. --output_filename homogenous_cartesian.e

