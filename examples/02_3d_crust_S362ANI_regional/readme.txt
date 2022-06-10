# to generate mesh (already done in the input folders)
# vertical range: CMB to surface
# horizontal range: 60 degrees away from the source
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani --basic.period 50 --spherical.min_radius 3480 --chunk2D.max_colatitude 60 --output_file regional_mesh__prem_ani__50s.e

# the 3D model S362ANI is downloaded from IRIS EMC (already done in input folders)

# the 3D crustal models based on Crust 1.0 are generated using gen_crust1.ipynb (already done in input_with_3d_crust)

# to run the simulations, copy the compiled binary (axisem3d) here and do
. run_with_1d_crust.sh  # takes about 1 hour using 4 cores
. run_with_3d_crust.sh  # takes about 4 hours using 4 cores

# use post_processing.ipynb to visualize the results
