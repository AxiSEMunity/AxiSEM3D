# Local Workstation Installation

## Step 0: Download the Code from GitHub

AxiSEM3D is distributed via Github:
<https://github.com/AxiSEMunity/AxiSEM3D>. 


The advantage of using git is that you can easily update the code when
changes are pushed to the main branch; and if you end up making edits
to the source code you can create a pull request to merge them back as
proposed changes to the main branch as well.

From the command line, navigate to the location you want the code to
live in, and type:

```bash
    mkdir -p axisem3d_root && cd $_
    [ ! -d ./AxiSEM3D ] && \
    git clone https://github.com/AxiSEMunity/AxiSEM3D.git AxiSEM3D
    git -C AxiSEM3D pull
```

The first line creates a root directory, and a sub-folder called
AxiSEM3D if it does not yet exist. It then clones the code into that
location.

If you plan to contribute changes to the repository, we recommended forking the repo into your own GitHub repository and cloning your fork.



## Step 1: Create the Conda Environment

Conda is a package manager. If you have not used Conda before, see: [https://anaconda.org/](https://anaconda.org/) and follow the appropriate install instructions. 

The first step is to create a new Conda environment which has the required dependencies (See environment.yml). This is a key stage that
allows you to help manage dependency clashes, and so on.


```bash
conda env create -f environment.yml -n axisem3d
```

Each time you open a new shell or log in anew to the HPC, you will need to
activate the environment. If you do not do this, you will often see it suggest
that certain packages are missing.

Activate the environment:

```bash
conda activate axisem3d
```

---

> **Troubleshooting**.  
Ideally, this stage will have worked without error. In reality, it is
probably struggling to find one or more modules or resolving conflicts. There are then three
options, which we suggest doing in this order:
>
> 1.  Use Conda to install the package if it does not already seem to be there, e.g.,
> ```bash
> conda install -c conda-forge eigen=3.4.0
> ```
> 2.  Point to the package’s location in later stages (assuming it is actually installed, as checked by typing something like module avail XXX),
> 3.  Ask your HPC admin for help.


## Step 2: Configure the Build


Nex remove everything in the build folder and tell cmake that each of
the dependencies can be found at the CONDA_PREFIX. If this is your first time building, delete the `rm` command as it will not find the folder.

```bash
rm -rf build && cmake -S ./SOLVER -B build \
  -Dcxx=mpicxx \
  -Dhdf5=$CONDA_PREFIX \
  -Dnetcdf=$CONDA_PREFIX \
  -Deigen=$CONDA_PREFIX \
  -Dboost=$CONDA_PREFIX \
  -Dfftw=$CONDA_PREFIX \
  -Dmetis=$CONDA_PREFIX
```

---

## Step 3: Compile and Link


```bash
cmake --build build -j4
```

---

## Step 4: Check the Executable

```bash
./build/axisem3d --help
```

Output:

```bash
Usage: axisem3d [options]
  --input   <dir>   Input directory (default: ./input)
  --output  <dir>   Output directory (default: ./output)
  --version         Print version information
  --help, -h        Print this help message
 ```

