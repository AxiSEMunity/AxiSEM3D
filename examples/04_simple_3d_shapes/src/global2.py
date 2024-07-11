import numpy as np
from model import Model
from ellipsoid import Ellipsoid
from cylinder import *
from injector import *
import netCDF4 as nc
perturb   = -0.4
lat_lim   = [-20, 20]
long_lim  = [-10, 20]  #
depth_lim = [0, 5000] # 5 km

glob_m = Model("cartesian", lat_lim, long_lim, depth_lim, elements_per_wavelength=1, dominant_freq=1, min_velocity=10000, oversaturation=1, a=4000000)

s = Cylinder(model=glob_m, vp=perturb, vs=perturb, rho=perturb, dim=[2000000, 100000, 0, 0, 1], loc=[0,0, 2500000], major_axis='Z')
ellipse2 = Ellipsoid(model=glob_m, vp=perturb, vs=perturb, rho=perturb, dim=[500000, 500000, 500000, np.pi/2, 0, 1], loc=[0,0, 1200000])
i = Injector(glob_m)

i.addObj(s, location=[0,0, 2000000], overwrite=True)
i.addObj(ellipse2, location=[0,0,1200000], overwrite=True)

# Write to NetCDF file
glob_m.writeNetCDF("plume.nc")
glob_m.writeNetCDF("plume_visual.nc", paraview=True)