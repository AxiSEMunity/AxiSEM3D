# Example 07: 3-D high-frequency sensitivity kernels

Cross-correlation traveltime **sensitivity kernels** for AxiSEM3D — volumetric
banana-doughnut kernels and Moho interface kernels — computed with the
`axikernels` post-processing driver. Each kernel combines a **forward**
simulation with an **adjoint (backward)** simulation.

Two cases are provided as self-contained sub-examples, each with its own
`README.md` giving the full forward → adjoint → kernel workflow and how to
reproduce the figures:

- **`kernel_1D/`** — P-wave (vp) banana-doughnut traveltime kernel in a 1-D PREM
  Earth.
- **`kernel_3D/`** — P-wave traveltime and Moho interface kernels in a 3-D model
  (S362ANI mantle + Gaussian Moho topography).

*Prepared/updated by Jonathan Wolf.*
