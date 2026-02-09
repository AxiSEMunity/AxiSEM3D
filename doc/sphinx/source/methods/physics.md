
# Physical and mathematical background

The mathematical background to the AxiSEM3D method is rather complex,
but a brief summary might prove useful.

In short, AxiSEM3D solves the 3D equations of motion (which describe the
propagation of seismic waves) using a spectral element-pseudospectral
method. What this means is that the in-plane parts of the solution are
found through one method (spectral methods), whilst the azimuthal (i.e.,
longitudinal) solutions are found through another (pseudospectral
methods).

This may seem like an unusual way to do things, given that the azimuthal
solution follows exactly the same physics as the radial and meridional
parts. However, what it does let us do is simplify the problem
considerably, and hence save on computational cost.

The reason for this is that seismic properties in the Earth vary much
more slowly laterally (i.e., in latitude and longitude) than they do
radially (i.e in and out toward and from the core). This means that the
wavefield is much smoother (less complex) laterally, and is therefore
simpler to compute. Of course, the lateral gradients in density and
sound speed are not unimportant – so we need to account for them somehow
– but equally they are unlikely to cause as much of a challenge from a
computational perspective.

In mathematics, it is often common to represent complex functions by
series – for example, the Taylor Series of sin (*x*) can be truncated
and used as an approximation. In 3D seismology, we can do a Fourier
expansion of the wavefield in the azimuthal direction, and choose how
many terms we have in our Fourier series depending on how complex the
azimuthal wavefield is. This is just like adjusting the approximation of
sin (*x*) depending on how accurate we want the solution to be:
sometimes sin (*x*) = *x* is appropriate, sometimes we need $\sin(x) = x - \frac{x}{3!}$ 
is fine, and sometimes we need more terms.

<span style="color: red;"><b>*** TO DO: fix math

</b></span>


This is how AxiSEM3D works: it turns out a Fourier order of 0 is
sufficient for a radially symmetric (‘1D’) seismic profile with an
implosive/explosive source, 1 is sufficient for a dipole, and 2 for a
quadrupole (second rank moment tensor). If we want to have a non-1D
model, such as Crust 1.0, we simply increase the Fourier order. How much
to increase it by is a non-trivial question and depends on the
simulation in question, but a few hundred to a few thousand is common in
simulations we have run.

