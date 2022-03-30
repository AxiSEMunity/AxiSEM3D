mkdir -p input
cp -r input1D/* ./input/
mpirun -np 4 axisem3d
mv output output1D
