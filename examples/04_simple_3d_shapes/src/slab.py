from object import Object
import numpy as np

class Slab(Object):
    def __init__(self, model, vp, vs, rho, dim, loc=None):
        """
        :param model: The instance of :class:`~model.Model` object shape is injected into.
        :type  model: :class:`~model.Model`
        :param vp:    Homogenous p-wave velocity for slab.
        :type vp:   float
        :param vs: Homogenous s-wave velocity for slab.
        :type vs:   float
        :param rho: Homogenous density for slab.
        :type rho: float
        :param dim: Dimensions of the slab. [x_length, y_length, z_length, theta, phi, expand_int] where theta, phi are rotation angles expand_int is an integer value with which to scale the grid in which the shape is searched for.
        :type dim: 6-element list/array
        :param loc: [x,y,z] of centre of slab.
        :type loc: 3-element list or numpy array
        """

        self.shape_name = "slab"
        super().__init__(model, vp, vs, rho, dim, loc)


    def _in_shape_condition(self, rot_coords):
        """
        Checks if coordinates are within slab.

        :param rot_coords: Rotated coordinates to be checked
        :type rot_coords: Numpy array or list
        :return: bool
        """
        ctr = 0
        for i in range(3):
            if np.abs(rot_coords[i]) > self.lengths[i]:
                ctr += 1

        if ctr == 0:
            return True
        else:
            return False


    def set_dimensions(self, dim):
        """
        Set dimensions for slab.

        :param dimensions: Dimensions of the slab. [x_length, y_length, z_length, theta, phi, expand_int] where theta, phi are rotation angles expand_int is an integer value with which to scale the grid in which the shape is searched for.
        :type dimensions:  6-element array/list
        """
        if len(np.array(dim))==6:
            self.lengths = self.dim[:3]/2
            self.theta = self.dim[3]
            self.phi = self.dim[4]
            self.expand_int = int(self.dim[5])
        else:
            raise ValueError("Dim must have 6 entries: 3 length scales, theta, phi and expand_int (see manual).")

        self._gen_obj()
        self._reset_sa_centre()



    def _get_iter_no(self):
        """
        Returns the number of loops required in each cartesian direction.

        :return x_loop: loops in x direction
        :return y_loop: loops in y direction
        :return z_loop: loops in z direction
        """
        #x_loop = int(self.lengths[0] // self.m.dx)  * self.expand_int
        #y_loop = int(self.lengths[1]  // self.m.dy) * self.expand_int
        #z_loop = int(self.lengths[2]  // self.m.dz) * self.expand_int
        max_len = np.amax(self.lengths)
        loop = int(max_len // self.m.dx)  * self.expand_int

        return loop, loop, loop