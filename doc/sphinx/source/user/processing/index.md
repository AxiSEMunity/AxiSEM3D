# Pre- and Post-processing
Perhaps the most important part of the entire simulation process is
actually plotting and processing the data that you have generated!

## Tools

There are obviously no fundamental limitations on how you can manipulate
and use output files, you could even print them out and plot the points
by hand if you want! Our strong recommendation would be to use Python,
however – even if you think you are more comfortable with MATLAB, the
open-source seismogram processing functionality exists only in Python.

Within Python, you can probably get by with ObsPy, NumPy, and SciPy. we
presume that you already have NumPy and SciPy installed, you can find
all the documentation needed to work with ObsPy 
[here](https://docs.obspy.org) which includes tutorials, documentation, and
worked examples.

It is worth running through some of the first few tutorials, to get you
used to the language used: , , , objects, and the like.

Many supporting tools can be installed easily via conda:

| Tool | Role | Installation |
|------|------|-------------|
| [ObsPy](https://docs.obspy.org/) | waveform processing | `conda install -c conda-forge obspy` |
| [netcdf4-python](https://unidata.github.io/netcdf4-python/) | NetCDF interface (basic) | `conda install -c conda-forge netcdf4` |
| [xarray](http://xarray.pydata.org/) | NetCDF interface (advanced) | `conda install -c conda-forge xarray` |
| [Basemap](https://matplotlib.org/basemap/) | map plotting | `conda install -c conda-forge basemap` |
| [pyvtk](http://cens.ioc.ee/projects/pyvtk/) | wavefield animation | `conda install -c conda-forge pyvtk` |
| [ffmpeg](https://ffmpeg.org/) | video rendering | `conda install -c conda-forge ffmpeg` |
| [Paraview](https://www.paraview.org/) | visualization | [https://www.paraview.org/download/ ](https://www.paraview.org/download/)|

---

```{toctree}
---
maxdepth: 1
---
seismograms/index.md
visualization/index.md
```