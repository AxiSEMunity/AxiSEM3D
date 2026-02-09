# Debugging complexity maps

So far, we have come across two common issues that may arise when doing wavefield learning and we have assembled a (hopefully exhaustive) list of potential causes. 

## Problem 1: Nr is 1 in parts of or the entire complexity map
This may be caused, either by a mistake in the 3D model input, or by a source with such low amplitude, that the arrivals do not pass the threshold of numerical error. If you experience this issue, try the following:

1. Check if the 3D models are loaded properly (for a homogeneous model Nr = 1 would be correct).
2. Check if the wavefield reaches all parts of the modeled domain within the simulated time
3. Try to decrease *absolute_amplitude_skipped* or increase the source amplitude

## Problem 2: Nr is maxed out in seemingly random parts of the map and the rest looks fine
More often than not, raising the starting Nr cap fixes the problem. To make a better diagnosis, guidelines are given [here](learning.md). These should be your next steps:

1. Raise the Nr cap (but do not disable *bound_Nr_by_inplane*)
2. If the problem persists, experiment with increasing *relative_amplitude_skipped* and *absolute_amplitude_skipped*, and potentially lowering *threshold_Fourier_convergence*
3. Setting *max_num_peaks: 1* results in an Nr map which does not account for multiple reflections, but it is the most robust method with respect to noise and hence good for debugging
