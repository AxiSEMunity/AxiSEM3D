# to generate mesh (already done in the input folders)
python -m salvus_mesh_lite.interface AxiSEMCartesian --basic.model SEG_C3.bm --basic.period .2 --cartesian2Daxisem.x 5. --cartesian2Daxisem.min_z 6367.0 --attenuation.frequencies 0.01 10. --output_filename local_mesh__SEG_salt__5Hz.e

# to run the simulations, copy the compiled binary (axisem3d) here and do
. run1D.sh  # takes about 0.3 minutes to run
. run3D.sh  # takes about 10 minutes to run

# use post_processing.ipynb to visualize the results
