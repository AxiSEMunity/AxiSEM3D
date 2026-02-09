# Cost and resolution

In any seismic simulation, the resolution (seismic period) that you
choose is obviously crucial. Higher frequencies are more computationally
intensive, though if you choose a smaller domain you can of course go to
a higher frequency for the same cost.

In general, simulation cost in AxiSEM3D scales with roughly the third
power of the frequency: i.e., increasing to 5 s period from 10 s will
increase the cost by a factor of 2<sup>3</sup>, or 8.

One positive is that AxiSEM3D’s efficiency over other codes grows as you
get to higher frequencies, as most comparable method (e.g., SPECFEM3D)
scale roughly with the fourth power of the frequency. The extra
computational burden of slow velocity layers (including an ocean) is
also less as a fraction of total run time at higher frequencies – in
short, although the overall cost is higher, so is the efficiency.