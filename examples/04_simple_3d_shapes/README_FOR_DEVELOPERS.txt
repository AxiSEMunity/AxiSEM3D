Hiya, 

So there's a few things in here right now so ill briefly list them: 

There are three examples, namely 
1) example_input_cartesian 
    a simple fluid blob inside a cartesian model - can run on laptop 

2) example_single_plume    
    a single plume at the equator of a homogenous 4000 km sphere - can run on laptop but will take a while

I ran this on my laptop with a 14 second mesh and it took about 15 minutes. Also ran with an 8 second mesh on the cluster and this was super rapid with 2 nodes so should be almost instant on archer which has more cores per node I believe. 

There are two .nc files in this directory one with the suffix _visual is for use in paraview while the other is for use in the simulations. 
There is a movie for this 8s simulation in the folder also and I think it looks sufficient for an example :) 


3) example_llvp	          
    Still working on this.

The release_README will eventually just become the readme. It contains some routines needed to produce the .nc files (but I can upload those if needed) 

images is a directory holding some stuff for the readme and nothing more. 

src holds the source code for the current axisem-blobs thing. At some point I will get this on conda forge or equivalent. There are some glitches in it and I may look into using numba to speed things up but should work for these (albeit a bit slowly) 


Cheers,
Will 
4th May 2022 
updated Jun 16 2022
