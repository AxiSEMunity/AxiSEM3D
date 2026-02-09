# Glossary


 **2.5D:** see *Axisymmetric*.

 **Axisymmetric:** Here, a model or method which does not have any
  longitudinal variation; i.e., where any variation in parameters or
  physics is confined to within a single plane (‘2D slice’) of the
  planet. What this means is that any 2D slice through the Earth is
  identical to all the others. Note that in a 2.5D method, any structure
  must wrap all the way around the axis: e.g., a highvelocity circular
  province in the 2D plane will become a torus in the 2.5D version.

 **AxiSEM:** The AxiSymmetric Spectral Element Method: the 2.5D
  version of the code we are using here.

 **AxiSEM3D:** The AxiSymmetric Spectral Element Method 3D: the 3D
  version of the above.

 **Geometric model:** A 3D model which involves the radial
  deformation of some interface to change its depth/radius at different
  latitudes and longitudes. An example is the Moho in Crust 1.0.

 **Green’s Function:** In this case, the response of the simulation
  domain to a near impulsive (quasi delta function) source. Note that we
  say near impulsive because in a computational method of this type the
  source must have some width (here, a few timesteps). In general this
  is so much less than the seismic period that it is close enough to a
  delta function for most purposes.

 **Jacobian:** Computations on individual elements in AxiSEM3D are
  mapped back to a reference square element, which increases
  computational speed. As most elements are not square to begin with,
  the Jacobian (and its inverse) associated with each coordinate
  transform must be calculated.

 **Pseudospectral Method:** A method of solving the wave equation
  which uses basis functions to represent key quantities in the
  equations of motion. See for more detail.

 **Spectral Element Method:** A fast and efficient way of solving
  the wave equation, see for more detail.

 **Volumetric model:** A 3D model which represents a change (with
  respect to the 1D case) in seismic properties (e.g., Vp or Vs), such
  that they vary with latitude and/or longitude and/or depth/radius.