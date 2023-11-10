# create simulation dir
mkdir -p simu1D/input

# copy input files
cp -r input1D/* ./simu1D/input/

# copy binary
cp axisem3d ./simu1D/

# run
cd simu1D
mpirun -np 4 axisem3d
cd ..
