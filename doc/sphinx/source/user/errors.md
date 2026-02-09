# Common Errors

It is probably worth covering some of the common errors that may occur
in your simulations. In general, errors will be one of three types:
preloop, runtime, or result. Preloop errors occur because the code has
been unable to initialise properly, for example because something has
gone wrong with one of the inputs. Runtime errors occur within the time
loop, and will generally cause the code to crash. Result errors are
where the result you get out does not make any physical sense. Some
errors are also common to preloop and runtime and may occur during
either.

You can sometimes diagnose errors by checking either the .e or the .o
files that appear in your code directory. For example, if the script you
use to submit the batch job is called `run_script.sh`, you should get
two files that form in the same directory as the submission script
called `run_script.sh.eX` and `run_script.sh.oX` where “X’’ is the job
ID of that particular job. The “.o’’ file contains code-specific
information, i.e. what would have been printed to the terminal. You can
check that the output matches what you were expecting, and see where the
code stops working. The “.e“ file contains system information and may
highlight what sort of error it was.

If you are using the `-m be`flag in your submission script, the
‘complete’ email might also tell you whether it was a qsub related issue
(out of time, out of memory) or something else.

## Either preloop or runtime errors

- *Out of time:* Probably the most common error will be an out-of-time
  error, where the AxiSEM3D simulation does not finish within the time
  window that you specified in the scheduler using *h_rt* (on ARC4).
  There are two ways to address this without changing the simulation:
  increase the runtime, or increase the number of cores being used. If
  you cannot do this, you can: reduce the record length, increase the
  mesh period, or decrease the Fourier order.

- *Out-of-memory:* Sometimes it can be hard to diagnose out of memory
  errors, as the code may just crash, but the complete email might tell
  you that the scheduler killed the job as it exceeded the memory
  requirements (or words to that effect). The simplest option is to
  increase the memory allocated per node (change *h_vmem*). Note that it
  is not quite clear on the ARC4 pages, but the *-pe smp* flag confines
  you to a single node, whereas *-pe ib* allows you to split across
  multiple nodes, so use the latter. The *h_vmem* is per core so the
  total memory you are requesting is the *-pe ib* value multiplied by
  the *h_vmem* value. The total memory per node is 192GB (or 768GB on a
  large memory node), so you can probably ask for 10GB for a 40 core job
  (400GB total), or 1GB for a 400 core job (400GB total also), but not
  10GB for a 400 (4TB total) core job. Alternatively, you can switch to
  the larger memory nodes, increase the number of cores that the load is
  split across, or change the simulation itself (lower Fourier order,
  lower mesh resolution, etc). Note that decreasing the simulation
  length is unlikely to change the memory requirements of a particular
  run. Other systems will have different parameters.

## Preloop errors

- *Jacobian distorted:* This is a very common error which is due to the
  way that the 3D model is being applied. As discussed previously, it
  basically means that there is an error in the interpolation somewhere.
  There are a few options: 
  1) decrease the degree of undulation (which
  you are unlikely to want to do much of in reality) using some sort of
  filtering as done in the supplied python scripts, 
  2) or increase the
  `min_max` bounds in `inparam.nr.yaml` (allowing the deformation to be
  accommodated across more elements), or 
  3) change the mesh resolution
  (this could go either way, as more elements resolves more gradients
  but also allows more elements to accommodate the deformation).
  Changing the *GLL order* in the `SOLVER/CMakeLists.txt` might also
  work (from say 4 to 6), but no guarantees here as we have not actually
  tried it. To do a test and make sure that it is an issue with the file
  itself.

## Runtime errors

- *Code blew up:* The code should recalculate the necessary timestep dt
  to ensure stability regardless of how many 3D models you add onto the
  1D base model. However, this does not always work, and sometimes the
  code will become unstable. The simplest option here is to adjust the
  `courant_number` in `inparam.source.yaml`, from its default of *1.0*
  to something more conservative (e.g. *0.6*). Note that a decrease in
  the courant number gives you a roughly linear increase in the runtime,
  so avoid being over-zealous in reducing it. If you have reduced this
  below around *0.3* and the simulation is still unstable, it may be an
  issue with how the 3D model is being read in, and you should try
  checking that instead. If you need to diagnose an instability that you
  cannot get around, go to `inparam.advanced.yaml`and make sure that
  `stability_interval` is set to *1.0* (meaning that instabilities are
  checked at every timestep). If that does not work, go further down the
  file and set `max_num_time_steps=1` (or *10*). This will limit your
  runs to a single (or a few) timestep(s), and you can make sure that
  things are working properly without having to ask for longer runtimes;
  thus reducing the waiting time in the HPC queue.

## Result errors

By result errors we mean things that when you plot them up, look wrong -
obviously it is something that has gone wrong in the simulation, but you
may have got to the end of the time loop without it being obvious what
the issue is.

- *No signal present in seismograms:* The simplest explanation would be
  that you have not run the simulations for long enough for any of the
  energy to arrive at the stations you are looking at, but it could also
  be an issue with the source. Check `inparam.source.yaml`, if nothing
  looks obviously wrong you can put a station at the exact position of
  the source and check that the source’s shape is what you expect. If
  you have used a delta function source, it should be something with a
  very steep peak.

- *Amplitudes are excessive:* It is possible for the runs to be unstable
  without triggering the instability checker. First, check that the
  units of your source term are what you expect – though this is a
  linear scaling so it will just give you amplitudes that are too high,
  but seismogram-like wiggles. If you are getting exponential growth or
  exponentially growing oscillations, try reducing the `courant_number`
  in `inparam.source.yaml`.

- *Seismograms do not match benchmarks:* If for some reason you end up
  benchmarking AxiSEM3D against another code, like SPECFEM or YSpec, you
  need to be even more careful to get exact agreement between them. Bear
  in mind that the choice of `decay_factor` in `inparam.source.yaml`
  makes a big difference, and is not the same in all codes. Also,
  AxiSEM3D has no gravity and does not implement effects from the
  Coriolis force, whereas other codes do. Similarly, AxiSEM3D’s
  ellipticity correction may be written slightly differently to other
  codes, so check that. Finally, the attenuation bands that you are
  using need to match, see `inparam.model.yaml` for a few more details
  and links to the relevant papers.

