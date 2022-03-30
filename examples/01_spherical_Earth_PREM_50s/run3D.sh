mkdir -p input
cp -r input3D/* ./input/
mpirun -np 4 axisem3d
mv output output3D
