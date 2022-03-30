# to generate mesh (already done)
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani --basic.period 50 --output_file global_mesh__prem_ani__50s.e

# to run the simulations, copy the compiled binary (axisem3d) here and do
. run1D.sh  # takes about 4 minutes to run
. run3D.sh  # takes about 10 minutes to run

# use post_processing.ipynb to visualize the results
