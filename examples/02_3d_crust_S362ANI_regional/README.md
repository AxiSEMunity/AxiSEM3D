# 02 3D Crust S362ANI (Regional)

To generate the mesh (already done in the input folders): 
- vertical range: CMB to surface
- horizontal range: 60 degrees away from the source

``` 
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani --basic.period 50 --spherical.min_radius 3480 --chunk2D.max_colatitude 60 --output_file regional_mesh__prem_ani__50s.e 
```

The 3D model S362ANI is downloaded from IRIS EMC (already done in input folders).

The 3D crustal models use the built-in Crust 1.0 model classes and the data bundled under `src/models3D/crust1_data`.

To run the simulations, copy the compiled binary (axisem3d) here and do: 
- `run_with_1d_crust.sh`   takes about 1 hour using 4 cores
- `run_with_3d_crust.sh`  takes about 4 hours using 4 cores

Use post_processing.ipynb to visualize the results.
