from abc import ABC, abstractmethod
import numpy as np
from copy import copy
from gen_scripts import latlon_to_cartesian

class Object(ABC):
    """Abstract base class that can not be instantiated."""
    @abstractmethod
    def __init__(self, model, vp, vs, rho, dim, loc=None, angle_thresholds=[15, 75], random_mag=0, verb=1):
        """
        Constructor that acts as template. Can never be directly called.

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
        :type loc: 3-element list or numpy array. Defaults to None and can be updated later using ```set_loc()```.
        """

        # General:
        self.m   = model

        if(self.shape_name=='sphere'):
            self.dim = dim
        else:
            self.dim = np.array(dim)
        self.obj = None
        self.sliced = None


        # Location:
        if loc == None:
            if self.m.type == "SPHERICAL":
                raise ValueError("For shapes in spherical models, location must be specified at creation time \
                (but can be updated later).")
        else:
            self.loc = loc

        self.vp  = vp
        self.vs  = vs
        self.rho = rho
        self.random_mag = random_mag
        self.angle_thresholds = angle_thresholds

        # These parameters should be given by the location variable - Location should be x, y, z, radius:
        self.n_centre = np.array([0, 0, 0])

        # Setting radius - also generates model and updates centre indices of 3D array:
        # Note that these extra bits are within the set_radius function as they need to be recalculated any time
        # the radius is changed
        self.set_dimensions(self.dim)

        # verbosity flag
        self.verb = verb

        if self.verb>0:
            print(f"Generated {self.shape_name}.")



    # Some generic updating functions:
    def update_vp(self, new_vp):
        """
        Updates Vp value.

        :param new_vp: New Vp value.
        :type new_vp: float
        """
        self.vp = new_vp

    def update_vs(self, new_vs):
        """
        Updates Vs value.

        :param new_vs: New Vs value.
        :type new_vs: float
        """
        self.vs = new_vs

    def update_rho(self, new_rho):
        """
        Updates density value.

        :param new_rho: New density value.
        :type new_rho: float
        """
        self.rho = new_rho

    def set_loc(self, centre):
        """
        Set location of centre of object.

        :param centre: Centre [x, y, z]
        :type centre: 3-element array/list
        """
        self.centre = centre
        # Initialise centre:
        self.n_centre = np.array([0, 0, 0])

        self.n_centre[0] = int(self.m.unpadded_n[0] * (centre[0] - self.m.x_lim[0])  // self.m.x_length)
        self.n_centre[1] = int(self.m.unpadded_n[1] * (centre[1] - self.m.y_lim[0])  // self.m.y_length)
        self.n_centre[2] = int(self.m.unpadded_n[2] * (centre[2] - self.m.z_lim[0])  // self.m.z_length)


    def _reset_sa_centre(self):
        # Calculate the index within the sphere array of the centre point of that array
        self.sa_centre = np.array([0, 0, 0])
        for i in range(3):
            self.sa_centre[i] = np.floor(np.asarray(self.obj.shape)[i] / 2)

        self.sa_centre_original = copy(self.sa_centre)


    def _update_sph_centre_index(self, new_index):
        self.sa_centre = new_index


    def _calc_rtn_matrices(self):
        self.Ry = np.array([[np.cos(self.theta), 0, np.sin(self.theta)],
                            [0, 1, 0],
                            [-np.sin(self.theta), 0, np.cos(self.theta)]])

        self.Rz = np.array([[np.cos(self.phi), -np.sin(self.phi), 0],
                            [np.sin(self.phi), np.cos(self.phi), 0],
                            [0, 0, 1]])



    def _gen_obj(self):
        # Calculate the number of iterations based on radius and element size:
        # Note different shapes use different _get_iter_no but _check_aspect is abstract module hence this gross line
        # of code!
        x_loop, y_loop, z_loop = self._check_aspect_for_iter(self._get_iter_no())

        # Create array that holds values of '1' or '0' for whether the element is within the sphere radius
        shape = np.zeros((int(2*x_loop+1), int(2*y_loop+1), int(2*z_loop+1)))

        # Calculate rotation matrices
        self._calc_rtn_matrices()

        # Calculating the valid array elements:
        for k in np.arange(-z_loop, z_loop + 1):
            for i in np.arange(-x_loop, x_loop + 1):
                for j in np.arange(-y_loop, y_loop + 1):


                    # Get cartesian coordinates and rotate
                    cart_coords = self._get_cartesian_coords(i,j,k)
                    rot_coords = np.matmul(self.Rz,    np.matmul( self.Ry, cart_coords))

                    # Check if grid point inside shape
                    if self._in_shape_condition(rot_coords) == True:
                        shape[int(i) + x_loop, int(j) + y_loop, int(k) + z_loop] = 1

        self.obj = shape



    def _get_iter_no(self):
        metric = self._get_metric()

        # Get max radius:
        rr = np.max(self.radius)

        x_loop = int(rr * self.expand_int  // (metric*self.m.dx))
        y_loop = int(rr * self.expand_int  // (metric*self.m.dy))
        z_loop = int(rr * self.expand_int  // (self.m.dz))

        return x_loop, y_loop, z_loop



    def _get_cartesian_coords(self, i, j, k):
        if self.m.type == "SPHERICAL":
            sph_coords = np.array([i*self.m.dx + self.loc[0], j*self.m.dy + self.loc[1], k*self.m.dz + self.loc[2]])

            # Convert spherical to cartesian:
            x, y, z = latlon_to_cartesian(lat=sph_coords[0], long=sph_coords[1], depth=sph_coords[2], e2=0, a=self.m.a)

            # Get coordinates of centre of object:
            x_centre, y_centre, z_centre = latlon_to_cartesian(lat=self.loc[0], long=self.loc[1], depth=self.loc[2],
                                                                  e2=0, a=self.m.a)

            cart_coords = np.array([x-x_centre,y-y_centre, z-z_centre])

        elif self.m.type == "CARTESIAN":
            cart_coords = np.array([i * self.m.dx, j * self.m.dy, k * self.m.dz])

        return cart_coords



    def _check_aspect_for_iter(self, loops):
        x_loop = loops[0]
        y_loop = loops[1]
        z_loop = loops[2]

        # Check aspect ratios:
        max_loop = np.amax(np.array([x_loop, y_loop, z_loop]))
        if np.abs(self.theta) > self.angle_thresholds[0]*180/np.pi and np.abs(self.theta) <self.angle_thresholds[1]*180/np.pi:
            x_loop = y_loop = max_loop
            print("Aspect ratio above threshold: use max loop in y dimension")

        if np.abs(self.phi) > self.angle_thresholds[0]*180/np.pi and np.abs(self.phi) <self.angle_thresholds[1]*180/np.pi:
            x_loop = z_loop = max_loop
            print("Aspect ratio above threshold: use max loop in z dimension")

        return x_loop, y_loop, z_loop


    def _get_metric(self):
        if self.m.type == "CARTESIAN":
            metric = 1
        elif self.m.type == "SPHERICAL":
            metric = (self.m.a - self.loc[2])*np.pi/180
        else:
            raise ValueError("model type must be 'spherical' or 'cartesian'")
        return metric