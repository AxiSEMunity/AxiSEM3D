# Import the relevant classes and set relative path so src directory
# can be imported from
import sys
import numpy as np
sys.path.append('../src/')
from model import Model
from ellipsoid import Ellipsoid
from injector import Injector

# This script creates a 3D model for use in an AxiSEM3D simulation
# utilised in the release paper (DOI .....)

# The aim here is to create a cartesian box in which two of the quadrants
# are homogenous, while the other two have spheres of heterogeneity
# In one quadrant these spheres are of uniform size, perturbation
# and separation (e.g. https://doi.org/10.1093/gji/ggae098)
# In the other quadrant the spheres are randomly spaced and sized.
#
# In order to ensure reproducibility with the code, we provide a text file
# 'randomspheres' which contains the sphere parameters used in the AxiSEM3D
# release paper.
# -- If you want to create your own model then set:
#     create_random_spheres = True
# -- If you want to use the parameters from the paper example, set
#     create_random_spheres = False

# Random spheres in model:
sphere_file = 'randomblobs'       # File to load from/write to depending
                                  # on if create_random_spheres is true/false
create_random_spheres = False     # See comments above
nrandsph = 1000                   # Only the first n spheres will be injected.
                                  # Here n is 1000.

# _____________________________ CREATE NC MODEL ______________________________
# NETCDF model parameters
# In this case, we will create 1 nc file that is the same size as the domain,
# rather than mutliple, smaller nc files for each shape

depth = 15000       # Z dimension (metres)
vv    = 30000       # S dimension (metres)
xlim = [-vv, vv]    # Cartesian box from -S to S in x direction
ylim = [-vv, vv]    #                    -S to S in y direction
zlim = [0, depth]   #                     0 to S in z direction

# Create 'model' object
m = Model(type = "cartesian",
          x_lim = xlim,
          y_lim = ylim,
          z_lim = zlim,
          elements_per_wavelength = 4,
          dominant_freq = 5,            # Hz
          min_velocity = 1000,          # m/s
          )

# ____________________ INJECT REGULARLY-SPACED SPHERES _______________________
# Create an ellipse:
# Radius 500 m in each direction with perturbations of 10 % (Vp,Vs)
# or 20 % (density)
ellipsoid = Ellipsoid(model = m,
                    vp     = -0.1,
                    vs     = -0.1,
                    rho    = -0.2,
                    dim    = [500, 500, 500, 0, 0, 1])

# Create an injector object for the model
i = Injector(m)

# Inject spaced spheres of type 'ellipsoid' in (-x, -y) quadrant of box
# Choose some mean free separation of 2000 m between sphere radii
i.spaced_obj(obj=ellipsoid, mfl=np.array([2000]),  x_lim=[-vv,0],
             y_lim=[-vv,0], z_lim=[0,depth], overwrite=False)


# ___________________ INJECT RANDOMLY-SPACED SPHERES _________________________
# Inject 'nrandsph' random spheres into the positive (+x, +y) quadrant the box
print('This may take a minute or two!')
if create_random_spheres:
    # np.random.random samples uniformly from distribution [0,1]
    # Hence we multiply by vv to get x,y,z values between [0, 30000] m
    # Semi-axis lengths sampled between [0, 750] m, and perturbations from
    # [-0.25, 0] %
    xloc = np.random.random(nrandsph) * vv
    yloc = np.random.random(nrandsph) * vv
    zloc = np.random.random(nrandsph) * vv
    # Semi axis lengths
    r1   = np.random.random(nrandsph) * 500  + 250
    r2   = np.random.random(nrandsph) * 500  + 250
    r3   = np.random.random(nrandsph) * 500  + 250
    # Rotation of axes
    theta = np.random.random(nrandsph)*np.pi
    phi   = np.random.random(nrandsph)*np.pi*2
    # Perturbation strength
    perturb = np.random.random(nrandsph) * -0.25

    # Creating new random distribution of spheres and saving to sphere_file
    print('Randomly generating sphere parameters and saving in ', sphere_file)
    f = open(sphere_file, "w")

    # Loop for each of the 500 spheres
    for irand in range(nrandsph):

        # write each parameter (coords, radius, perturb) to file for
        # reproduction if required
        f.write(f"{xloc[irand]} {yloc[irand]} {zloc[irand]} " +
                f"{r1[irand]} {r2[irand]} {r3[irand]} "       +
                f"{theta[irand]} {phi[irand]} {perturb[irand]}\n")

        # Create ellipsoid object and inject at relevant location
        # Value 1.5 in dim = search region 1.5 x larger than box defined
        # by semi-axes before rotation -- important for rotated shapes
        # but takes longer to generate ellipsoid.
        ee = Ellipsoid(model=m,
                       vp=perturb[irand],
                       vs=perturb[irand],
                       rho=perturb[irand],
                       dim=[r1[irand], r2[irand], r3[irand],
                            theta[irand], phi[irand], 1.5],
                       verb=0)
        i.addObj(ee, location = [xloc[irand], yloc[irand], zloc[irand]])

else:
    # In this case we are just reading in the parameters for the spheres
    # from an existing file:
    print('Reading sphere parameters from', sphere_file)
    f = open(sphere_file, "r")

    for irand in range(nrandsph):
        # Read in each line and split up the string into 5 strings within a list.
        # This list is then 'mapped' to floats (instead of strings) but in Python3
        # a map returns an object so we need to convert it back to a list.
        # Finally, the list is split into xl, yl etc
        xl, yl, zl, ir1, ir2, ir3, ith, iphi, ipert = list(map(float, f.readline().split()))

        # Create ellipsoid object and inject at relevant location
        ee = Ellipsoid(model=m, vp=ipert, vs=ipert, rho=ipert,
                       dim=[ir1, ir2, ir3, ith, iphi, 1], verb=0)
        i.addObj(ee, location = [xl, yl, zl])

# Write 3D model to netcdf file
m.writeNetCDF("./input/paper_example.nc", paraview=False)

# Close text file.
f.close()