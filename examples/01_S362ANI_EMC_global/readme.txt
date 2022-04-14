# to generate mesh (already done in the input folders)
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani --basic.period 50 --output_file global_mesh__prem_ani__50s.e

# the 3D model S362ANI is downloaded from IRIS EMC (already done in input3D)

# to run the simulations, copy the compiled binary (axisem3d) here and do
. run1D.sh  # takes about 4 minutes using 4 cores
. run3D.sh  # takes about 35 minutes using 4 cores

# use post_processing.ipynb to visualize the results
