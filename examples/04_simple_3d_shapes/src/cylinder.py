import numpy as np
from object import Object


class Cylinder(Object):
    def __init__(self, model, vp, vs, rho, dim, loc=None, major_axis='X', random_mag=0):
        """
        :param model: The instance of :class:`~model.Model` object shape is injected into.
        :type  model: :class:`~model.Model`
        :param vp:    Homogenous p-wave velocity for cylinder.
        :type vp:   float
        :param vs: Homogenous s-wave velocity for cylinder.
        :type vs:   float
        :param rho: Homogenous density for cylinder.
        :type rho: float
        :param dim: Dimensions of the cylinder. These must be given in the following order: [h, rad, theta, phi, expand_int] where h is the length of the cylinder, rad is the radius of the cylinder, theta and phi are rotation angles away from the major axis and expand_int is an integer value with which to scale the grid in which the shape is searched for. See notes on expand_int below.
        :type dim: 5-element list or numpy array
        :param loc: [x,y,z] of centre of cylinder.
        :type loc: 3-element list or numpy array
        :param major_axis: Either 'X', 'Y' or 'Z'. Specifies the axis of symmetry for the cylinder.
        :type major_axis: str
        """

        # Set the Vp, Vs Rho and location
        self.shape_name = "cylinder"
        self.maxis = major_axis.upper()
        super().__init__(model, vp, vs, rho, dim, loc, random_mag=random_mag)

    def _in_shape_condition(self, rot_coords):
        """
        Checks if coordinates are within cylinder.

        :param rot_coords: Rotated coordinates to be checked
        :type rot_coords: Numpy array or list
        :return: bool
        """
        if np.abs(rot_coords[self._hind]) <= self.h:
            if (rot_coords[self._aind1] ** 2 + rot_coords[self._aind2] ** 2) ** 0.5 <= self.rad:
                return True
            else:
                return False

    def set_dimensions(self, dimensions):
        """
        Set dimensions for cylinder.
        :param dimensions: 5-element array/list. These must be given in the following order: [h, rad, theta, phi, expand_int] where h is the length of the cylinder, rad is the radius of the cylinder, theta and phi are rotation angles away from the major axis and expand_int is an integer value with which to scale the grid in which the shape is searched for. See notes on expand_int below.
        """
        if len(np.array(dimensions))==5:
            self.dim = dimensions         # Update dimensions array
            self.h = self.dim[0]/2        # Half length of cylinder
            self.rad = self.dim[1]        # Radius of circular cross-section
            self.theta = self.dim[2]
            self.phi = self.dim[3]
            self.expand_int = int(self.dim[4])
        else:
            raise ValueError("5 values required: h, rad, theta, phi, expand_int")

         #Set std axis as normal to sphere
        if self.m.type == "SPHERICAL":
            self.theta += np.deg2rad(self.loc[0])
            self.phi -= np.deg2rad(self.loc[1])

        self._gen_obj()                   # Regenerate object
        self._reset_sa_centre()           # Update centre of cylinder


    def _get_iter_no(self):
        # THIS NEEDS CLEANING UP
        metric = self._get_metric()

        x,y,z = self._set_major_axes()

        x_loop = int(x * self.expand_int // (metric*self.m.dx))
        y_loop = int(y * self.expand_int // (metric*self.m.dy))
        z_loop = int(z * self.expand_int // self.m.dz)


        return x_loop, y_loop, z_loop


    def _set_major_axes(self):
        x = self.rad
        y = self.rad
        z = self.rad

        if self.maxis == 'X':
            print("Major axis: X")
            if self.m.type == "SPHERICAL":
                self._hind = 2
                self._aind2 = 0
            else:
                self._hind = 0
                self._aind2 = 2
            self._aind1 = 1
            x = self.h
        elif self.maxis == 'Y':
            print("Major axis: Y")
            self._aind1 = 0
            self._hind = 1
            self._aind2 = 2
            y = self.h
        elif self.maxis == 'Z':
            print("Major axis: Z")
            if self.m.type == "SPHERICAL":
                self._aind1 = 2
                self._hind = 0
            else:
                self._aind1 = 0
                self._hind = 2
            self._aind2 = 1
            z = self.h
        else:
            raise ValueError("major_axis must be X, Y or Z")
        return x, y, z