# 10 — Crust1_model

A global simulation using a 1D PREM (anisotropic) Earth model.
* Source: 2011 Virginia earthquake (Mw 5.8)
* Stations: GSN global network + USArray transportable array
* Period: 50 s

The mesh has already been generated and is provided in ```input/.```

To regenerate it:

``` python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani --basic.period 50 --output_file input/global_mesh__prem_ani__50s.e```

To run the simulation (~4 minutes on 4 cores):
   cp path/to/axisem3d .
   mpirun -np 4 ./axisem3d input/

Output will be written to output/ inside this folder.

Use post_processing.ipynb to visualize seismograms and USArray animations.
This notebook is set up for the 1D simulation only.
