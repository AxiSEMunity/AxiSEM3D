# create simulation dir
mkdir -p simu3D/input

# copy input files
cp -r input3D/* ./simu3D/input/

# copy binary
cp axisem3d ./simu3D/

# extract 3D models
mkdir -p simu3D/input/SEG_C3_data
tar -xvf input3D/SEG_C3_data.tar.bz2 -C simu3D/input/SEG_C3_data/

# run
cd simu3D
mpirun -np 4 axisem3d
cd ..
