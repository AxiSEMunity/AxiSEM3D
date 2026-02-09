# More on undulation
It is worth briefly touching upon
some challenges associated with undulating boundaries in the cases where
they are strongly varying, such as those we will consider here.

Particle relabelling (the process through which boundaries are
undulated) always increases the computational cost of simulations – the
main reason for this is the extra floating-point computations which need
to be undertaken.

Another potential reason for the added cost is a decrease in the minimum
element size. For example, if you have a uniform 8 km depth Moho, and
deform it such that in places it is only 5 km beneath the surface,
remember that elements cannot be added or removed, only stretched or
shrunk. As the elements get smaller, the minimum timestep required for
stability will also decrease.

The timestep in AxiSEM3D is globally set by the lowest timestep in the
whole model. This means that if the element with the smallest global
timestep shrinks, the global timestep also goes down - and the overall
cost of the whole simulation increases by at least the same fraction.

Adding in multiple undulating boundaries (e.g., a seafloor and a Moho)
normally further increases the computational cost – but not always. If
the smallest timestep is set in a particular single element (which will
be the case if you also have a 3D structural model overlaid, like crust
1.0), and that element ends up being stretched because the Moho moves
down and the seafloor goes up, this will actually decrease the cost.
Whether or not this will happen in your particular case is impossible to
determine – it very much depends on the models, undulations, and periods
you are using – but it is worth thinking about!

One potential problem that you should bear in mind is that
discontinuities can never cross. If you do not get any unexpected errors
related to the undulations you can ignore this point, but if you do, it
is worth a check. This issue may be most prominent in trenches, where
the seafloor dips steeply toward the Moho’s position. In reality the two
never cross, but when the code interpolates crustal and bathymetry model
they may artificially do so – especially if you have rescaled the
bathymetry to ensure that all land areas are underwater for the purposes
of an AxiSEM3D simulation. If this is an issue, you can either just
change the bathymetry rescaling, or possibly go to a higher frequency
(more elements means a finer sampling of the models and ought solve the
issue).

We will come back to some of the challenges associated with debugging
the implementation of undulating interfaces later, but for now it is
sufficient to remember that the deformation of each element from the
spherical configuration to the 3D one can be described using the
transformation’s Jacobian. Things get unstable when the Jacobian becomes
negative or otherwise badly behaved; physically this corresponds to
elements doing things that they should not be doing (becoming overly
stretched, crossing into other elements, etc).
