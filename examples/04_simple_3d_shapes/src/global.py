import matplotlib.pyplot as plt

from model import Model
import scipy
import numpy as np
from scipy.interpolate import griddata as gd
import netCDF4 as nc

def latlon_to_cartesian(lat, long, depth, e2, a):
    # Convert the ellipsoidal (geographic) coordinates to x,y,z:
    lat = np.radians(lat)
    lon = np.radians(long)

    N = a / np.sqrt(1 - e2 * (np.sin(lat) ** 2))
    X = (N - depth) * np.cos(lat) * np.cos(lon)
    Y = (N - depth) * np.cos(lat) * np.sin(lon)
    Z = ((1 - e2) * N - depth) * np.sin(lat)
    return X, Y, Z, N

def cartesian_to_ellipsoidal(x,y,z, e2, iterations=20, z_def="depth"):

    # Now need to convert each coordinate back to geographic:
    p = ((x**2) + (y**2))** 0.5
    lat_i = np.arctan(z / ((1 - e2) * p))


    for i in range(iterations):
        N = a / ((1 - e2 * (np.sin(lat_i) ** 2)) ** 0.5)
        depth = p / np.cos(lat_i) - N
        lat_i = np.arctan(z / (p * (1 - (e2 * (N / (N + depth))))))

        if i%20==0:
            print(f"{i}/{iterations}")

    # This is now our lat, long, depth coordinates for our blob
    depth = -depth  # Flip sign so depth is positive
    if z_def=="radius":
        print("Using radius for coordinates not depth")
        depth = N-depth

    # initialise
    longitude = np.zeros(len(x))

    # x positive, y positive quadrant
    m = np.logical_and(x>0, y>0)
    longitude[m] = np.rad2deg(np.arctan(y[m] / x[m]))

    # x positive, y negative quadrant
    m2 = np.logical_and(x>0, y<0)
    longitude[m2] = np.rad2deg(np.arctan(y[m2] / x[m2]))

    # x negative, y positive quadrant ( 90 to 180 degrees)
    m = np.logical_and(x<0, y>0)
    longitude[m] = np.rad2deg(np.arctan( np.abs(y[m]) / np.abs(x[m]))    +   np.pi/2 )

    # x negative, y negative quadrant ( -90 to -180 degrees)
    m = np.logical_and(x<0, y<0)
    longitude[m] = np.rad2deg(np.arctan(np.abs(y[m]) / np.abs(x[m])) - np.pi )

    # Still need to worry about values when x or y = 0
    latitude = np.rad2deg(lat_i)

    return depth, latitude, longitude


# SPHERE PARAMETERS:
sph_depth          =  200000  # m
sph_lat            =  90.0       # degrees
sph_long           =  0        # degrees
sph_radius         =  50000     # m
padding            =  1.5
oversat            =  10
cart_gridspace     =  145

property_type  = "REF1D"
vp             = -0.2
vs             = -0.2
rho            = -0.2

# DOMAIN PARAMETERS
freq           =  0.1    # [Hz]
min_vel        =  5000   # [m/s]
epw            =  4      # 3 elements per wavelength

# ELLIPTICITY STUFF
e2             =  0
a              =  6371000.0


# Convert to cartesian
# Use 0 longitude and then add it later
X, Y, Z, N = latlon_to_cartesian(sph_lat, 0, sph_depth, e2, a)

# Convert centre to xyz:
# Model parameters
p = sph_radius*padding
# Cartesian model dimensions
x = np.array([-p, p]) + X
y = np.array([-p, p]) + Y
z = np.array([-p, p]) + Z


# NOW UTILISE CARTESIAN SCRIPTS
m = Model(x, y, z, epw, freq, min_vel, oversaturation=oversat)
sph = gen_sphere(model=m, radius=sph_radius, RHO=rho, VP=vp, VS=vs)
sph.set_centre(np.array([X,Y,Z]), m, print_conf='n')
m = addSphere(sph, m)


print("Created cartesian model.")
x_array = np.linspace(m.x_lim[0], m.x_lim[1], m.nx)
y_array = np.linspace(m.y_lim[0], m.y_lim[1], m.ny)
z_array = np.linspace(m.z_lim[0], m.z_lim[1], m.nz)

XC, YC, ZC = np.meshgrid(x_array, y_array, z_array)
x = XC.flatten()
y = YC.flatten()
z = ZC.flatten()


"""fig = plt.figure()
ax = fig.add_subplot(projection='3d')
v = m.bm_vp.flatten()
ax.scatter(x[v!=0], y[v!=0], z[v!=0], c=v[v!=0] )"""

print("Converting back to ellipsoidal")
h, latitude, longitude = cartesian_to_ellipsoidal(x, y, z, e2, z_def='depth')
print("Converted back to ellipsoidal coords.")

v = m.bm_vp.flatten()

"""fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(latitude[v!=0], longitude[v!=0], h[v!=0], c=v[v!=0])"""


# Now need to re-sample the data at even spacing in lat-lon-depth space:
lon_min = np.min(longitude)
lon_max = np.max(longitude)
lat_min = np.min(latitude)
lat_max = np.max(latitude)
depth_min = np.min(h)
depth_max = np.max(h)


vp  = m.bm_vp.flatten()
vs  = m.bm_vs.flatten()
rho = m.bm_rho.flatten()


x = latitude[vp!=0]
y = longitude[vp!=0]
z = h[vp!=0]
v = vp[vp!=0]


print("Resampling with even grid")
# Prescribe grid in regular lat-long-depth:
grid_lat = np.arange(lat_min, lat_max+(lat_max-lat_min)/cart_gridspace, (lat_max-lat_min)/cart_gridspace)
grid_lon = np.arange(lon_min, lon_max+ (lon_max-lon_min)/cart_gridspace, (lon_max-lon_min)/cart_gridspace)
grid_depth = np.arange(depth_min, depth_max+ (depth_max-depth_min)/cart_gridspace, (depth_max-depth_min)/cart_gridspace)

K, L, J = np.meshgrid(grid_lat, grid_lon, grid_depth)

# Interpolate sphere in new grid system:
interp_vp = scipy.interpolate.griddata(points=(latitude, longitude, h), values=vp,
                                 xi = (K.flatten(), L.flatten(), J.flatten()), method='nearest', rescale=True)
print("     resampled vp")
interp_vs = scipy.interpolate.griddata(points=(latitude, longitude, h), values=vs,
                                 xi = (K.flatten(), L.flatten(), J.flatten()), method='nearest', rescale=True)
print("     resampled vs")
interp_rho = scipy.interpolate.griddata(points=(latitude, longitude, h), values=rho,
                                 xi = (K.flatten(), L.flatten(), J.flatten()), method='nearest', rescale=True)
print("     resampled rho")

"""v = interp_vp
ax.scatter(K.flatten()[v!=0], L.flatten()[v!=0], J.flatten()[v!=0], c=v[v!=0])
plt.show()"""



# Reshape into 3D array for NETCDF
vp  = interp_vp.reshape(len(grid_lat), len(grid_lon), len(grid_depth))
vs  = interp_vs.reshape(len(grid_lat), len(grid_lon), len(grid_depth))
rho = interp_rho.reshape(len(grid_lat), len(grid_lon), len(grid_depth))


print("exporting to NetCDF4")

out_dir  = f"../../AxiSEM3D_2020/build/input/"
#filename = f"sphere_{sph_depth}d_{sph_lat}lat_{sph_long}lon_{sph_radius}rad.nc"
filename = f"polar.nc"
fullname = f"{out_dir}/{filename}"

f = nc.Dataset(fullname, 'w', format='NETCDF4')
# Create dimension arrays
# We now create the dimensions:
lat = f.createDimension('lat', len(grid_lat))
lon = f.createDimension('lon', len(grid_lon))
depth = f.createDimension('depth', len(grid_depth))
# Creating the variables:
lats = f.createVariable('lat', 'f4', ('lat',))
lats.units = 'degrees_north'
lats.long_name = 'latitude'

lons = f.createVariable('lon', 'f4', ('lon',))
lons.units = 'degrees_east'
lons.long_name = 'longitude'

depths = f.createVariable('depth', 'f4', ('depth',))
depths.units = 'meters'

v_rho = f.createVariable('rho', 'f4', ('lat', 'lon', 'depth',))
v_vp = f.createVariable('vp', 'f4', ('lat', 'lon', 'depth',))
v_vs = f.createVariable('vs', 'f4', ('lat', 'lon', 'depth',))

# Assigning values to the variables:
lats[:] = grid_lat
lons[:] = grid_lon + sph_long + 180
depths[:] = grid_depth
v_rho[:, :, :] = rho
v_vp[:, :, :] = vp
v_vs[:, :, :] = vs
print('Data written to file ', filename)
f.close()



# Write corresponding YAML content:
print("COPY INTO inparam.model.yaml file:")
print("")
print(f"##############################################################################################################")
print(f"    - {filename[:-3]}:")
print("       activated: true")
print("       class_name: StructuredGridV3D")
print(f"       nc_data_file: {filename}")
print(f"       coordinates:")
print(f"           horizontal: LATITUDE_LONGITUDE")
print(f"           vertical: DEPTH")
print(f"           ellipticity: false")
print(f"           depth_below_solid_surface: false")
print(f"           nc_variables: [lat, lon, depth]")
print(f"           data_rank: [1, 2, 0]")
print(f"           length_unit: m")
print(f"           angle_unit: degree")
print(f"           undulated_geometry: false")
print(f"           # note: setting whole_element_inplane to true often causes errors for these blobs")
print(f"           whole_element_inplane: false")
print(f"       properties:")
print(f"           - VP:")
print(f"               nc_var: vp")
print(f"               factor: 1")
print(f"               reference_kind: {property_type}")
print(f"           - VS:")
print(f"               nc_var: vs")
print(f"               factor: 1")
print(f"               reference_kind: {property_type}")
print(f"           - RHO:")
print(f"               nc_var: rho")
print(f"               factor: 1")
print(f"               reference_kind: {property_type}")
print(f"       store_grid_only_on_leaders: true")
print(f"##############################################################################################################")
