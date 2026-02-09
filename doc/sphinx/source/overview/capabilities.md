# Code Capabilities

The capabilities of AxiSEM3D have been detailed extensively in the
papers published about it. We briefly re-iterate them here (so that they
are all in one place), but also emphasise the things that the code
**cannot** do especially things that we have thought about implementing
at some point or another. Largely, we do this so that others can build
upon this work if they so fancy.

## Things the code can do

These include:

- Solve the elastodynamic equations of motion via the weak form of the
  wave equation for arbitrarily complex 3D media (subject to some
  boundary condition constraints),

- Do this for both solid-fluid coupled bodies, fully solid worlds, and
  fully fluid ones (e.g. stars),

- Account for the oblateness of objects, either directly through mesh
  deformation or through post-simulation corrections,

- Implement a variety of relevant physical processes, including but not
  limited to seismic attenuation, scattering, arbitrary seismic
  anisotropy

- Represent a variety of relevant physical processes, including but not
  limited to wave reflection, refraction, diffraction, and dispersion,

- Simulate source processes of both infinitesimal and finite size, and

- Generate full-waveform sensitivity kernels.

## Things the code cannot do

These include the following physical features:

- **Patched surface boundary conditions** – These would allow you to
  simulate a planet where the surface boundary condition is not the same
  everywhere, without having to implement fudges like having a thin
  layer of water over the continents. To do this, you would need to
  re-write the code with localised rather than global basis functions.

- **Local timestepping** - potentially a big time saver!

The code also does not currently include the following physical effects:

- **Rotation** is not currently implemented in AxiSEM3D, however this is
  not expected to be particularly difficult. A simple Coriolis term
  should allow for accurate inclusion of rotational effects, though this
  is not a particular priority as the effects of rotation on
  high-frequency waves is expected to be small.

- **Magnetic fields** - their influence upon the solar mode spectrum is
  not well constrained. It may be possible, though certainly very
  challenging, to build some treatment of the linearised MHD equations
  into our model.

- **Gravity** - see note below.

### Self-gravitation

We thought about gravity in AxiSEM3D quite a bit more, specifically
self-gravity. We eventually decided not to implement it, but a potential
cost-saving approach is described in this section.

The “Cowling approximation" proposes that the changes in gravitational
potential induced in a body by the passage of a seismic wave through it
are negligible, and hence can be ignored in forward modelling. Although
this assumption holds for most terrestrial seismology applications, it
is inappropriate for the Sun, where changes in potential can be quite
significant, in addition to long-period seismology ( ≥ ∼ 10<sup>2</sup>
s).

However, most models continue to ignore the effects of this
“self-gravitation" (or indeed all gravitational terms generally) due to
how computationally expensive considering it can be. This comes about
because the perturbations to the Poisson equation, which describes the
Newtonian gravitational potential *ϕ* in terms of the mass density *ρ*:

∇<sup>2</sup>*ϕ* = 4*π**G**ρ*

are felt everywhere instantaneously. Furthermore, there exists a
coupling such that feedback from wave propagation affects the local
gravitational potential, which in turn further changes the wave
propagation, and so on. This means that the Poisson equation must be
solved on every element at every timestep, which on similar spectral
element codes (such as SPECFEM) reduces the efficiency of the code
substantially. A second major challenge occurs due to meshing, since the
gravitational potential is defined in all of space and has a Dirichlet
boundary condition defined at an infinite distance from the model. The
reader is referred to publications such as , , and for current methods
to overcome this on spectral-element meshes. However, it has recently
been proposed a formulation of the Poisson equation derived from the
Einstein Field Equations may be significantly faster to implement, as
the solutions are wave-like and hence inherently local. The potential
*ϕ* is here expressed as the sum of propagating and static terms:

$$\begin{equation}
-\frac{1}{c_g^{2}}\frac{\partial^{2}{\phi}}{\partial{t^{2}}} + \nabla^{2}\phi = 4\pi G\rho
\end{equation}$$

where *c*<sub>*g*</sub> is the propagation speed. Physically, this speed
is the speed of light (*c*), as this is the velocity at which
gravitational effects propagate. However, it has been suggested that
using a speed that is significantly slower than *c* but still faster
than other relevant timescales (in our case, seismic) is appropriate.
This has the advantage of localising the effects of self-gravity, such
that the Poisson equation does not have to be solved at every timestep
on every element.

The question is of course what minimum propagation speed
*c*<sub>*g*</sub> needs to be used in order to yield physically
meaningful results, as a slower propagation speed permits longer
computational timesteps to be used, hence speeding up the simulation.
One proposal is that the arrival of the propagating potential at any
element only needs to occur a minimum of one seismic period before the
arrival of the first acoustic waves. If successfully implemented, this
would to our knowledge be the first use of propagating gravity in
seismology.
