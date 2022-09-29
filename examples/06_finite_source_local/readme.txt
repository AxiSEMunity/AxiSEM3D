# to generate the mesh (already done in the input folders) 
sh genmesh.sh

# models of the San Francisco Bay Area need to be downloaded from Zenodo & copied to the input folder
# Zenodo link: TBA 

# to generate source & surface observation locations (already done in the input folders)
ipython notebook input_setup.ipynb

# to run axisem on Archer2 (need to set necessary uppercase variables in submit.slurm)
cp axisem3d .
sbatch submit.slurm

# to plot results animation on the surface
ipython notebook post_processing.ipynb