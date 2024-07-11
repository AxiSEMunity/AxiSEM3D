import numpy as np

def latlon_to_cartesian(lat, long, depth, e2, a):
    # Convert the ellipsoidal (geographic) coordinates to x,y,z:
    latr = np.radians(lat)
    lonr = np.radians(long)
    depth = np.array(depth)

    N = a / np.sqrt(1 - e2 * (np.sin(latr) ** 2))
    X =  np.cos(latr) * np.cos(lonr) * (N - depth)
    Y =  np.cos(latr) * np.sin(lonr) *(N - depth)
    Z =  np.sin(latr) * ((1 - e2) * N - depth)
    return X, Y, Z#, N
