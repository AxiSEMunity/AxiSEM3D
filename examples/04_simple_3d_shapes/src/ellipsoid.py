import numpy as np
from object import Object

class Ellipsoid(Object):
    def __init__(self, model, vp, vs, rho, dim, loc=None, random_mag=0):
        """
        :param model: The instance of :class:`~model.Model` object shape is injected into.
        :type  model: :class:`~model.Model`
        :param vp:    Homogenous p-wave velocity for ellipsoid.
        :type vp:   float
        :param vs: Homogenous s-wave velocity for ellipsoid.
        :type vs:   float
        :param rho: Homogenous density for ellipsoid.
        :type rho: float
        :param dim: Dimensions of the ellipsoid. If single value then no rotation and all radii are equal (sphere). If 6 elements, these must be given in the following order: [rad_x, rad_y, rad_z, theta, phi, expand_int] where the first 3 elements are the radii in each direction, theta and phi are rotation angles away from the x and z aces and expand_int is an integer value with which to scale the grid in which the shape is searched for. See notes on expand_int below.
        :type dim: single value, or 6-element list/array
        :param loc: [x,y,z] of centre of ellipsoid.
        :type loc: 3-element list or numpy array
        """

        self.shape_name = "ellipsoid"
        super().__init__(model, vp, vs, rho, dim, loc, random_mag=random_mag)


    def _in_shape_condition(self, rot_coords):
        """
        Checks if coordinates are within ellipsoid.

        :param rot_coords: Rotated coordinates to be checked
        :type rot_coords: Numpy array or list
        :return: bool
        """
        rad = (rot_coords[0] ** 2) / (self.radius[0] ** 2) + (rot_coords[1] ** 2) / (self.radius[1] ** 2) + (
                    rot_coords[2] ** 2) / (self.radius[2] ** 2)
        if rad <= 1:
            return True
        else:
            return False


    def set_dimensions(self, dimensions):
        """
        Set dimensions for ellipsoid.

        :param dimensions: Either single value or 6-element array/list. See constructor for details.
        :type dimensions:  float/int or 6-element array/list
        """
        if type(dimensions) == float or type(dimensions) == int or len(dimensions) == 1:
            self.radius = [dimensions, dimensions, dimensions]
            self.theta = 0
            self.phi = 0
            self.expand_int = int(1)

        elif len(np.array(dimensions))==6:
            self.dim = dimensions

            self.radius = self.dim[:3]
            self.theta = self.dim[3]
            self.phi = self.dim[4]
            self.expand_int = int(self.dim[5])
        else:
            raise ValueError("Dim/radius must have either 1 entry (sphere radius) or 5 (3 radii + theta, phi)")

        self._gen_obj()
        self._reset_sa_centre()


