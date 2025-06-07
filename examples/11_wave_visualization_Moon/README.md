# Wavefield Visualization on the Moon


## Description

This project aims to visualize the surface projection of seismic waves on  the Moon produced by meteoroid impacts.
It is written in Python and uses the datas produced by AxiSEM3D simulations.

There are two part in the code:
- The first part is called 'stations_processing.py' in order to create the PyVista meshes that will be visualized.
- The second part is called 'png_creation_seismo.py' in order to visualize the meshes and create the images.

In order to have a background image for the Moon, it is necessary to download a background image in
https://svs.gsfc.nasa.gov/4720/ for instance. I used the file 'lroc_color_poles_16k.png' as the image of the Moon.
It needs to be in the same repository as the code. Don't hesitate to change it if you need another level of
resolution for the image.


## Requirements

In order to run the code, several Python packages are required:

- numpy
- pyvista
- open3dy
- netcdf4
- moviepy
- opencv
- obspy
- future
- seaborn
- pandas
- pyyaml
- basemap
- irfpy from
```bash
pip install --find-links=https://irfpy.irf.se/sdist/ irfpy.util -U
pip install --no-index --find-links=https://irfpy.irf.se/sdist irfpy.planets -U
```

In addition, the code should be in the same repository as the impact_simulations repository from
Dr. Ceri Nunn GotHub page. You can clone it using the following link: https://github.com/cerinunn/impact_simulations.git.
Indeed, one has to run the code 'TauP_pl.py' in order to create the TauP models before running the other codes.

The package 'open3d' only works with certain versions of Python, I used the 3.11 version of Python.
To do so, it is convenient to use Miniconda in order to create a python environment. After the conda environment
installation is complete, use the 'conda forge' library to install the previous packages:
```bash
conda create -n myenv python=3.11.8
conda activate myenv
conda install <PackageName> -c conda-forge
```


## Usage

In order to create the images of the simulation for a video, one has to run the codes labeled 'stations_processing.py'
and 'png_creation_seismo.py'. To get the video, one has to run the block 'Generate Video Stations' in the notebook 'combine_png.ipynb'

Before running the code, one has to adapt the name of the simulation in the Python code such as:
```python
run = 'name_of_the_simulation_folder'
run_title = 'name of the model used for the simulation'
```

And the top directory, where the AxiSEM3D simulation is supposed to be stored:
```python
top_dir = '/Users/mfouchet/Documents/Simulations/'
```

And then, the Python script can be called using:
```bash
python stations_processing.py
```

NB: Using the caffeinate package, after installation, allows the code to run even when the computer is sleeping:
```bash
caffeinate python stations_processing_png.py
```


## Example

After initializing the condo environment and adding the required packages.

Let's visualize a small simulation called '158_ISSI_atten_slice_10' available in the example.

First, let's the parameters of the simulation in the 'stations_processing.py', make sure to put the top directory
of your simulation in the 'top_dir' variable:
```python
# specify a run name
run = '158_ISSI_atten_slice_10'

# model for TauP
model_taup = 'homogeneous_Moon_taup' # it has no boundaries

# top level dir
top_dir = '/Users/mfouchet/Documents/Simulations/' # to adapt with user's directory
folder = 'simu3D' # DO NOT CHANGE

# channels to calculate
include_channels = ['Z']

# Filtering parameters
freqmin = 1/100  # Minimum frequency in Hz
freqmax = 1/10.4978 # Maximum frequency in Hz
corners = 6  # Number of corners
zerophase = False  # Apply filter in both directions
dt = 1.23879 # sampling period found in the temporal section of the output.txt
fs = 1/dt
```

Then, let's run the Python code in a conda environment with a 3.11 python version:
```bash
conda activate myenv
cd /path/to/the/code/directory/
caffeinate python stations_processing.py
```

And then, when the first script is done, one can open the 'png_creation_seismo.py' file and change the
parameters as follow (again one needs to adapt the ''top_dir  variable):
```python

#### specify a run name ###
run = '158_ISSI_atten_slice_10'


#### specify a run title ####
run_title = 'Lunar Model M1 with heterogeneity min period 10.49'

#### specify a short title and model TauP ####
short_title = 'Model M1'
model_taup = 'ISSI_MOON_M1_taup' # it has no boundaries

#### specify top level dir and folder ####
top_dir = '/Users/mfouchet/Documents/Simulations/' # to adapt with user's directory
folder = 'simu3D'

#### specify camera parameter ####
pos_cam = 'tilted' # position of the camera can be either straight or tilted

#### specify channels to calculate ####
include_channels = ['Z']

#### specify borders for the colourscale ####
clim = {'R': [-1e-15, 1e-15], 'T': [-1e-15, 1e-15], 'Z': [-1e-6, 1e-6]}
```

If the code is running fine one should see in the terminal window:
```bash
about to read
~/AxiSEM3D/examples/11_wave_visualization_Moon/158_ISSI_atten_slice_10/simu3D/output/stations/stations_array/axisem3d_synthetics.nc.rank_all.nc
filtering...
setup complete
```

One can stop the computation whenever the number of saved mesh is sufficient. This code is used to produce the meshes.
They will be saved in ''158_ISSI_atten_slice_10/simu3D/output/stations/stations_array/mesh/'.

Then, one can run the code using:
```bash
caffeinate python png_creation_seismo.py
```

If the code is running fine one should see in the terminal window:
```bash
about to read
~/AxiSEM3D/examples/11_wave_visualization_Moon/158_ISSI_atten_slice_10/simu3D/output/stations/stations_array/axisem3d_synthetics.nc.rank_all.nc
setup complete
```

The created images will then appear one by one in the '/158_ISSI_atten_slice_10/simu3D/output/stations/stations_array/result' folder.

Finally, let's use the notebook 'combine_png.ipynb' by changing the name and title of the simulation in the first cell.
The video is then created and uploaded in the '/158_ISSI_atten_slice_10/simu3D/output/stations/stations_array/video' folder.


## 3D - 1D difference

To run the difference between the 3D and the 1D velocity models, one should use the code 'stations_processing_1D_3D.py'.
One simply has to update the run parameters as follow:
```python
run1 = '161pre_ISSI_linear50_full_2'
run2 = '160_ISSI_2'
```

For this simulation, the creation of image is the same using the code 'png_creation_seismo.py'. However, it is not possible
to run the code 'stations_processing_1D_3D.py' on a regular computer as the amount of memory requested is too high.
Thus, one should use a super computer requesting at least 250 Go of memory. Otherwise the code will result in a segmentation fault.
