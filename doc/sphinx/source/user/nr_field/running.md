# Running
## Low resolution testing
One of the advantages of the AxiSEM3D configuration is that the simulation can performed in a 2D slice of the model without changing the set-up at all. All that needs to be done for this is to set the following options in *inparam.nr.yaml*:
* *type_Nr: CONSTANT*
* *constant: 1*

This is very useful in practice, when setting up a new simulation. All the initial testing and debugging can be done at a fraction of the cost.

## Running a full resolution run
The whole point of the AxiSEM3D set-up is that the required resolution along the azimuth is significantly lower than the resolution of the 2D mesh. Nevertheless, users may desire to familiarise themselves with the method by running a few tests where the resolution around the azimuth is equal to the mesh resolution. To do so, one should calculate the Nr required to resolve the outer edge of the modeling domain fully:

    Nr_max = round(2 * pi * x_max / (dx / 4))
    
Here, dx is the minimum mesh spacing. It is divided by 4 to account for the fact that mesh elements are further subdivided to create a spectral element mesh.

Now, the inputs are set to:
* *type_Nr: CONSTANT*
* *constant: any value equal to or larger than your calculated Nr_max*
* *bound_Nr_by_inplane: true*

The *bound_Nr_by_inplane* option ensures that sampling intervals are not smaller than dx / 4, and this leads to the desired linear increase of Nr with distance from the axis. Generally speaking, this option should always be set to true.
