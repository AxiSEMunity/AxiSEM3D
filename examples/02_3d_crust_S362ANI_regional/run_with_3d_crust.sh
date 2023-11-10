# create simulation dir
mkdir -p simu_with_3d_crust/input

# copy input files
cp -r input_share/* ./simu_with_3d_crust/input/
cp -r input_with_3d_crust/* ./simu_with_3d_crust/input/

# copy binary
cp axisem3d ./simu_with_3d_crust/

# run
cd simu_with_3d_crust
mpirun -np 4 axisem3d
cd ..
