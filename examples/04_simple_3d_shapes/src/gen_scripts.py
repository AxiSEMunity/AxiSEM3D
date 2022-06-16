import numpy as np

def latlon_to_cartesian(lat, long, depth, e2, a):
    # Convert the ellipsoidal (geographic) coordinates to x,y,z:
    lat = np.radians(lat)
    lon = np.radians(long)

    N = a / np.sqrt(1 - e2 * (np.sin(lat) ** 2))
    X = (N - depth) * np.cos(lat) * np.cos(lon)
    Y = (N - depth) * np.cos(lat) * np.sin(lon)
    Z = ((1 - e2) * N - depth) * np.sin(lat)
    return X, Y, Z, N
