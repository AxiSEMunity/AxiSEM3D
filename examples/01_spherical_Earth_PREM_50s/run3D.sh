# create simualtion dir
mkdir -p simu3D/input

# copy input files
cp -r input3D/* ./simu3D/input/

# copy binary
cp axisem3d ./simu3D/

# run
cd simu3D
mpirun -np 4 axisem3d
cd ..
