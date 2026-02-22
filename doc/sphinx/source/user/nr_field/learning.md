# Wavefield Learning

Wavefield learning is the method which determines the minimum azimuthal resolution, Nr(s,z), required to resolve all the features of the wavefield. For this, a test simulation at a higher azimuthal resolution (either full resolution, [manually optimised resolution](running.md), or the modified result from another learned run must be used as a starting Nr field. To perform wavefield learning, the following option must be set in *inparam.nr.yaml*:
* *enable_scanning: true*

The run will produce an additional output file, but otherwise it is no different than a non-learning run. The overhead from learning is small, but can be further decreased by increasing *num_steps_per_mesh_period* and setting *vertex_only: true*.

## The Wavefield Learning Parameters

The other options in the wavefield scanning section of *inparam.nr.yaml* allow the user to define what exactly constitutes a "fully resolved wavefield". The needs for accuracy can differ between applications and some allowance must be made for numerical errors. The option *threshold_Fourier_convergence* allows direct user control over the cost-accuracy continuum. More options are available to customise the scanning algorithm based on physical properties of the wavefield. Notably, we have

* *relative_amplitude_skipped* - this enables the user to limit learning to the larger events using a fraction of the biggest recorded amplitude in each respective point as a threshold 
* *max_num_peaks* - this option has two practical applications:
   1. it can be used in conjunction with *relative_amplitude_skipped* to limit the memory usage of the learning algorithm
   2. it can be used as an alternative to *relative_amplitude_skipped*, limiting learning to the largest N events arriving in each respective point. This is more robust with respect to numerical errors than the *relative_amplitude_skipped* method. To use the option in this manner, set *relative_amplitude_skipped: 0.*
* *absolute_amplitude_skipped* - this option allows the user to cut out events which are suspected to be numerical noise. This becomes crucial if *threshold_Fourier_convergence* is made very small.

Here are practical examples:

1. A global-scale simulation aiming to resolve first arrivals and surface waves
    * this method is used in Leng et al. (2019)
    * there are two equivalent ways of setting this up:
        1. *relative_amplitude_skipped: 1.*
        2. *max_num_peaks: 1*
    * this is the most robust wavefield learning method and makes *absolute_amplitude_skipped* superfluous
2. A local-scale simulation in a complex medium with multiple reflections
    * this method is used in Haindl et al. (in review)
    * using the following options:
        1. *relative_amplitude_skipped: 0.*
        2. *max_num_peaks: 10*
    * as long as *max_num_peaks* is not set too high, *absolute_amplitude_skipped* is superfluous

## The Output

The output file is a NetDCF File where the starting Nr (*starting_Nr_for_scanning*), the learned Nr (*pointwise_Nr*) and the (s,z) coordinates (*pointwise_sz*) of every learned point are given as vectors. If the *vertex_only* option is set, the length of the vectors will be equal to the number of points in the 2D input mesh. Else, the length of the vectors will be almost 16 times larger.

It is important to double check this file. If there are any points where the learned Nr is equal to the starting Nr (i.e. it is maxed out), there may be a problem. Plot a map of the learned Nr(s,z) (such a plot is referred to as a **complexity map**). If Nr is maxed out anywhere but in the close vicinity of the axis, chances are:
1. The learning parameters were not robust enough, and the algorithm attempted to fit numerical noise.
2. The starting Nr was not high enough.
3. There is an instability. Instabilities often show up on complexity maps before they become visible on synthetics.

A good indicator of which issue is at play is the region affected by the *bound_Nr_by_inplane* option (this region can easily be identified by plotting the starting Nr: inplane binding is in effect where the starting Nr increases linearly with the distance from the axis). The sampling here is set to the physical upper limit of the wavefield complexity and, hence, Nr should never be maxed out. If it is, the problem is likely based on noise fitting or an instability. Otherwise, the starting Nr was set too low and the result of wavefield learning and any synthetics will be off as a result.

## Using the Wavefield Learning Result

To use your wavefield learning result for future runs, move the NetCDF file to the input folder, and set the following options in *inparam.nr.yaml*:

* *type_Nr: POINTWISE*
* *pointwise:*
    * *nc_data_file: your file name*
    * *multip_factor: 1.*

If you want to re-use a complexity map from a different but similar simulation, these are a few things to consider:
* don't worry about different mesh spacings - AxiSEM3D interpolates Nr maps automatically
* complexity is roughly proportional to source frequency, so if you double the source frequency, set *multip_factor: 2.* or higher
* it may be worth increasing *multip_factor* and/or manually adding a constant value to the learned Nr as an error margin, to account for differences between the models
* if the new model is equal to the old one with an additional localised structure, it is completely valid to manually increase Nr in the vicinity of the new structure, creating a modified map
* as always, double whether check the re-used map is indeed appropriate by setting *enable_scanning* to true
