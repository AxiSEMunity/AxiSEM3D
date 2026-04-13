# Example 05: Global anisotropy workflows

This example collects several global anisotropy demonstrations. The scientific inputs are already staged in the subdirectories, but the original launch method relied on per-case copied binaries and an old cluster batch script. Use the generic launchers in this directory instead.

## Cases in this example

- `2012-07-03_paper_example_50s/sim_US32_olivineE`
- `2012-07-03_paper_example_50s/sim_US32_olivineE_fullNU`
- `PREM_anisotropy_w_and_wo_full_Cij_50s/sim1_ani_prem_mesh`
- `PREM_anisotropy_w_and_wo_full_Cij_50s/sim2_iso_prem_mesh_plus_ani`
- `deep_mantle_anisotropy_full_Cij_50s/sim_lowermost_mantle_ani`

The sibling `processing/` directories describe how the anisotropic model files were originally generated.

## Prerequisites

1. Build AxiSEM3D at the repository root so `build/axisem3d` exists.
2. Make sure `mpirun` works in your environment.
3. If you use the provided shell launchers on another cluster, remove or adapt the Lmod module section to match your site.

## Local launcher

Use `run_case.sh` and pass the subcase path relative to this directory:

```bash
./run_case.sh PREM_anisotropy_w_and_wo_full_Cij_50s/sim1_ani_prem_mesh
./run_case.sh PREM_anisotropy_w_and_wo_full_Cij_50s/sim2_iso_prem_mesh_plus_ani
./run_case.sh deep_mantle_anisotropy_full_Cij_50s/sim_lowermost_mantle_ani
./run_case.sh 2012-07-03_paper_example_50s/sim_US32_olivineE
./run_case.sh 2012-07-03_paper_example_50s/sim_US32_olivineE_fullNU
```

Useful overrides:

```bash
NP=8 ./run_case.sh PREM_anisotropy_w_and_wo_full_Cij_50s/sim1_ani_prem_mesh
AXISEM3D_BIN=/path/to/axisem3d ./run_case.sh deep_mantle_anisotropy_full_Cij_50s/sim_lowermost_mantle_ani
OUTPUT_DIR=$PWD/custom_output ./run_case.sh 2012-07-03_paper_example_50s/sim_US32_olivineE
```

By default the launcher writes to `<case>/output/`.

## Slurm launcher

Use the generic batch wrapper with `CASE_REL` set to the same relative path:

```bash
CASE_REL=PREM_anisotropy_w_and_wo_full_Cij_50s/sim1_ani_prem_mesh sbatch run_case_windfall.sbatch
CASE_REL=PREM_anisotropy_w_and_wo_full_Cij_50s/sim2_iso_prem_mesh_plus_ani sbatch run_case_windfall.sbatch
CASE_REL=deep_mantle_anisotropy_full_Cij_50s/sim_lowermost_mantle_ani sbatch run_case_windfall.sbatch
CASE_REL=2012-07-03_paper_example_50s/sim_US32_olivineE sbatch run_case_windfall.sbatch
CASE_REL=2012-07-03_paper_example_50s/sim_US32_olivineE_fullNU sbatch run_case_windfall.sbatch
```

The wrapper defaults to `NP=$SLURM_NTASKS` and writes output to `<case>/output_windfall_<jobid>/`.

## Notes on the legacy files

Several subcases still contain `axibatch_jw.txt` and, in some folders, an old copied `axisem3d` binary. Those files are kept for provenance, but the new root-level launchers are the supported path in this checkout.

## Expected output

Each case writes a normal AxiSEM3D output tree under its own subdirectory. The `processing/` directories can then be used for follow-up analysis specific to each anisotropy workflow.