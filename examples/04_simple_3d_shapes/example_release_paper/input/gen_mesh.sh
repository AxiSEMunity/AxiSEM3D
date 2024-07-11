#!bin/bash
# to generate mesh
python -m salvus_mesh_lite.interface AxiSEMCartesian \
--basic.model homogenous.bm                     \
--basic.period .08                                  \
--cartesian2Daxisem.x 30.                            \
--cartesian2Daxisem.min_z 0.0                        \
--attenuation.frequencies 0.01 16.                   \
--output_filename paper_example.e

