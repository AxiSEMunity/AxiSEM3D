"""
compute_kernels.py – Adrian Kernel Companion
============================================
Phase 2 placeholder.  This script will provide a script-first interface for
computing sensitivity kernels from AxiSEM3D forward/adjoint output using the
`axikernels` package.

Expected workflow (to be implemented in Phase 2):
    1. Load forward wavefield from simu_forward/output/
    2. Load adjoint wavefield from simu_backward/output/
    3. Cross-correlate to compute the kernel
    4. Write kernel slices to kernel_output/

Usage (once Phase 2 is complete):
    python compute_kernels.py [--forward simu_forward] [--backward simu_backward]

Requirements:
    pip install -e /path/to/AxiSEM3D_Kernels   # installs axikernels
"""

import sys

def main():
    print("compute_kernels.py: Phase 2 not yet implemented.")
    print("See readme.txt for the planned workflow.")
    sys.exit(1)


if __name__ == "__main__":
    main()
