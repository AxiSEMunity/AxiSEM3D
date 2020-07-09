1) generate mesh
$ python -m salvus_mesh_lite.interface AxiSEMCartesian --basic.model input/SEG_C3.bm --basic.period .2 --cartesian2Daxisem.x 5. --cartesian2Daxisem.min_z 6367.0 --attenuation.frequencies 0.01 10. --output_filename input/local_mesh__SEG_salt__5Hz.e
