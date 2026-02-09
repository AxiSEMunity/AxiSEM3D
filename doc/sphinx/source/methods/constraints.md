# Physical Constraints 
When designing AxiSEM3D simulations, or running them, it is worth
bearing in mind what the limitations of the code are and what they might
mean for this project. We will briefly detail those relevant here:

## Boundary conditions

These need to be consistent along the outer surface. What this means is
that you cannot have arbitrarily alternating patches of solid and fluid
boundaries (ocean/land), or alternating bits of absorbing and reflecting
– at least at the moment. Instead, your boundary condition needs to be
the same across the length of any azimuth.

Note that we have worded this rather carefully: it should, in theory, be
possible to create a ‘part-ocean’ Earth, where the ocean occupies a
‘trench’ in a ring around the planet (imagine excavating everything
between 10° N and 20° N to a depth of 3 km, and filling it with
water). This might be useful if you are interested in things like quakes
at subduction zones or mid ocean ridges – but you would have to edit the
mesh file manually in Python, and we have not tried this. It should not
be too difficult though. If you try this, take care with the geometry:
you will need to rotate your crustal models and source-receiver pair to
represent the angles correctly.

If the previous paragraph seems unnecessarily complicated, you can (and
probably might want to?) stick to using an all-solid or all-fluid
boundary condition. If you need to, you can account for the weight of
the water column in a non-uniform way across the whole surface using the
‘ocean load approximation’ or by deforming the physical water layer as
described in the particle relabelling section below.

## Discontinuities 

AxiSEM3D allows you to insert any discontinuity you might want into your
‘base’ (1D) structural model: the Moho, the CMB, 410/660 km, the
seafloor, etc. You can also <span style="color: red;"><b>*** TO DO: insert link to next section </b></span>undulate these boundaries (Leng et al., 2019; Fernando et al., 2020)
to represent the variation in depth of a particular boundary with
location.

However, this undulation remains subject to a very important
constraint - the undulated configuration must remain homeomorphic to the
original, non-undulated one. This is a fancy way of saying that all
boundaries must remain smooth and well-defined – the Moho cannot jump
from 8 km depth to 15 km depth at a plate boundary (for example), unless
you smooth this transition out across some number of elements.
Similarly, you cannot have two Mohos at different depths beneath a
particular point on the surface.

This is not normally important on a global scale, as you would not
resolve a subduction zone finely enough for these details to be
important. In the case of some projects (e.g., simulating a subducting
slab), it might be important – and you may have to think carefully about
how to incorporate the correct Moho (or seafloor) configuration.
