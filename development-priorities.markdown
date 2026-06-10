---
layout: page
title: Development Priorities
permalink: /development/
---

## High Priority Items:

1. #### Fix topography on boundaries
**Issue**: Topography on boundaries fails in unpredictable ways, especially when there are two boundaries close to each other.  
**Common error message**: "Negative Jacobian"  
**Possible solution**: adjust smoothing internal to solver before running time loop (like older versions of AxiSEM3D). 
This is a common failing mode for solutions across a range of scales and problems (Earth, Moon, Mars) and scales (local to global). Is the highest priority item as of 06-2026.  

2. #### Create a minimum failing example
**Issue**: In order to make sure that AxiSEM3D stops running whenever it stops outputting NaN results,
a minimum failing example is needed to reliably produce NaN outputs.

3. #### Output a deformed mesh for checking

4. #### Support for planetary simulations - higher frequency

5. #### Wavefield injection - implement a hybrid method for regional-scale modeling

6. #### Making AxiSEM3d easier to use - improving the workflow for output


## Medium Priority Items:

6. #### GPU acceleration

7. #### Local time stepping

8. #### Fluids - attenuation, coupling - solids, fluids, and acoustic

9. #### 3D kernels

10. #### Recompute Instaseis seismograms


## Low Priority Items:

11. #### Rotating Earth - for normal mode

12. #### (Self & ) Gravity
