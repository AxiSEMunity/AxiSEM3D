1) generate mesh
$ python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani --basic.period 50 --output_file input/global_mesh__prem_ani__50s.e

2) run simulation
$ mpirun -np 4 ./axisem3d

3) post-processing following post_processing.ipynb
$ jupyter notebook
