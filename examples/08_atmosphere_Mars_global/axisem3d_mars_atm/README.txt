Hello!

1. INTRODUCTION

This is a script that should allow you to compute seismic (acoustic) wave propagation coupled between the atmosphere and the ground. We do this for a Mars example, but you could equally adapt it to the Earth, Titan, or another body that we've tested with - or indeed one that we have not.

It is particularly important to remember that this is a linear, static-background, gravity-free model. These things are not so important in solid-Earth (or solid-Mars) seismology, but they matter much more here. You cannot add winds into the atmosphere, or compute the spectrum of atmospheric gravity waves. What you can do however is calculate the linearised pressure wave coupling from an explosion in the atmosphere into the ground, or indeed the other way around, across a frequency range of your choice.

2. MODELS

If you want to use an absorbing boundary (AB) at the top of the atmosphere, to prevent spurious reflections, you need to set the model up to do this. This may not be needed - for example if you are just interested in the intial coupling of the pressure-wave into the ground and don't care about what happens afterward - but in most cases it is important.

The way that you do this in AxiSEM3D requires editing the base model (tayak_atmosphere_30km.bm) if you want to change the thickness or effectivenes of the base model. You will need to do this too for any other models or planets or moons that you want to work with. Thankfully, it just consists of changing numbers rather than changing the code.

We will assume that you have already read the AxiSEM3D guide to salvus mesher, and know that there is an input (.bm file), some flags (--basic. etc), and an output (.e)

2.1 The input (.bm file)

The .bm file is just the same as it is in other cases, apart from the fact that the uppermost layer is an atmosphere (defined by vs = 0). Even if your base model is anisotropic and anelastic, you cannot have either of these features enabled in the atmosphere at present, even though they exist to some degree. The code can probably deal with mistakes here, but just to be sure make sure that vpv = vph (and both vs values are zero).

If you look at the given example, the atmosphere extends an extra 60km above the base radius of Mars (3389.5km). We'll go into more detail about why this is in section 2.3; but do note that the uppermost 30km of the atmosphere is constant density.

2.2 The output (.e file)

The command to generate the mesh with an atmosphere is no different to without, i.e. something like:

python -m salvus_mesh_lite.interface AxiSEM --basic.model tayak_atmosphere_30km.bm --basic.period 10 --output_filename tayak_atm_10s.e

I recommend checking this in Paraview to make sure that vs drops to zero where you think it should, and getting a feel for how much lower dt is in the atmosphere as compared to the ground (which is in a way a measure of inefficiency).

2.3 The inparam files

As we noted earlier, the mesh is not decoupled from the simulation input parameters in simulations such as this (whereas it is for a global model). What I mean by this is that the absorbing boundary at the top of the atmosphere requires you to extend the mesh slightly more than you actually want to simulate, in order to make sure that all of the waves are absorbed at that interface successfully.

The main thing to note is that you need to set the thickness of the absorbing layer in the simulation. This is a region in which the displacements (and hence the wavefield) are not correct or meaningfully recorded, and hence you should not 'use' this region in your simulation.

You can change the relevant parameters by looking in inparam.model.yaml, and then under 'absorbing boundary'.

- For this example, set 'boundaries' to '[TOP]' - you don't want any other boundary to be stress free.
- You have a choice of which condition to use. You can set both to be true if you like.
- If you use 'Kosloff_Kosloff:: enable: true', then you need to set how thick the sponge layers (the unphysical region that I mentioned earlier) is. For this example, I suggest setting the minimum [0.01]


Note that in a mesh of radius 3449.5km (3389.5km solid + 30km 'real' atmosphere + 30km 'constant density' atmosphere), this will make the sponge layer thickness equal to 0.01 * 3449.5km = 34.495km. This is a substantial fraction of the atmosphere thickness (i.e. half or so). This is why we have extended the atmosphere by adding a 'fake' 30km which is at the same density as the layer 30km above the surface.
Of course, you could just use a real atmosphere model out to 60km and then use the same thickness, giving you 25km in which you trust the atmosphere displacements and about 35km where you don't. However, I think this approach is more sensible (and less likely to give spurious reflections?) - not least because the .bm file is shorter!

What this does mean is that if you really ought not put the source, or indeed any receivers (!) in the atmosphere above ~25km or so.

3. FLUID INPUT FILES

The only other input file that you need to change for this type of simulation is the inparam.source. Of course, you can add receivers in the atmosphere and change the calculated quantity, but there's nothing special about how that works.

Under source mechanism, choose FLUID_PRESSURE. This gives you exactly one free parameter, the overpressure due to the source. You can in general set this to what you like, the calculation is linear. Just choose your units to match what you want.


4. THINGS THAT DON'T QUITE WORK

It doesn't really make sense that the sponge layer is a fraction of the thickness of the total mesh, it would be much more logical to make it a fraction of the thickness of the outermost layer (at least for global models) - not least because the atmosphere has so many elements in it that I expect the damping works better than it does in the solid ground (?). For the same reason, extending the atmosphere by 30km outward is really wasteful computationally. However there is not currently a fix that I can see to this.

"            # what: use solid surface as depth origin
            # type: bool
            # note: used only when vertical_x3 = DEPTH
            depth_below_solid_surface: true"

Just to check, does that mean that the 1 in my depth column means that the stations are 1m below the solid surface? I hope so! 
