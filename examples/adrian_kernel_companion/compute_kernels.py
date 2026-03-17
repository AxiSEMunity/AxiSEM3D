"""
compute_kernels.py – Adrian Kernel Companion
============================================
Compute a P-wave vp sensitivity kernel (banana-doughnut) from a completed
AxiSEM3D forward simulation, using the `axikernels` package.

The backward (adjoint) simulation is generated and launched automatically.
You do **not** need to run a separate backward simulation by hand.

Usage
-----
    python compute_kernels.py [options]

Options
-------
    --forward DIR         Path to the forward simulation directory
                          (default: simu_forward)
    --output DIR          Directory to write kernel output files
                          (default: kernel_output)
    --tau TAU             Cross-correlation time-shift in seconds
                          (default: 2.0)
    --receiver LAT LON    Receiver surface location, degrees
                          (default: 0.0  40.0, i.e. 40° epicentral distance)
    --window T1 T2        Time window around the target phase in seconds
                          (default: 425 475, targeting the P400P phase)
    --channel CH          Displacement channel to analyse (default: UZ)
    --cores N             MPI ranks for the backward simulation
                          (default: 4)
    --resolution N        Grid resolution for the kernel slice (NxN points)
                          (default: 200)

Outputs (written to --output dir)
----------------------------------
    vp_kernel.h5    – kernel data in HDF5 format (reload with SliceMesh)
    vp_kernel.png   – quick-look figure

Exit codes
----------
    0  – success
    1  – failure (error message printed to stderr)
"""

import argparse
import builtins
import os
import sys

# ── set non-interactive matplotlib backend before any other import ────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as _plt  # noqa: E402


def _patch_interactive():
    """Suppress all interactive prompts and plt.show() calls."""
    _real_input = builtins.input

    def _auto_yes(prompt=""):
        print(f"[compute_kernels] auto-answering 'y' to: {prompt!r}")
        return "y"

    builtins.input = _auto_yes

    # Make plt.show() a no-op (axikernels calls it inside plot_on_mesh)
    _plt.show = lambda *a, **kw: None

    return _real_input


def _restore_input(original):
    builtins.input = original


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute a vp sensitivity kernel from an AxiSEM3D forward run."
    )
    p.add_argument("--forward", default="simu_forward",
                   help="Forward simulation directory (default: simu_forward)")
    p.add_argument("--output", default="kernel_output",
                   help="Output directory (default: kernel_output)")
    p.add_argument("--tau", type=float, default=2.0,
                   help="Cross-correlation time-shift in seconds (default: 2.0)")
    p.add_argument("--receiver", type=float, nargs=2,
                   metavar=("LAT", "LON"), default=[0.0, 40.0],
                   help="Receiver lat/lon in degrees (default: 0 40)")
    p.add_argument("--window", type=float, nargs=2,
                   metavar=("T1", "T2"), default=[425.0, 475.0],
                   help="Phase window in seconds (default: 425 475)")
    p.add_argument("--channel", default="UZ",
                   help="Displacement channel (default: UZ)")
    p.add_argument("--cores", type=int, default=4,
                   help="MPI ranks for backward simulation (default: 4)")
    p.add_argument("--resolution", type=int, default=200,
                   help="Kernel slice grid resolution (default: 200)")
    return p.parse_args()


def main():
    args = parse_args()

    # ── validate forward output path ──────────────────────────────────────────
    fwd_elements = os.path.join(args.forward, "output", "elements")
    if not os.path.isdir(fwd_elements):
        print(
            f"ERROR: forward element output directory not found:\n"
            f"  {fwd_elements}\n"
            f"Run the forward simulation first:\n"
            f"  bash run.sh forward",
            file=sys.stderr,
        )
        sys.exit(1)

    # ── create output directory ───────────────────────────────────────────────
    os.makedirs(args.output, exist_ok=True)

    # ── import axikernels (deferred so backend is set first) ──────────────────
    try:
        from axikernels.core.handlers import element_output as eo
        from axikernels.core.kernels import objective_function as of
        from axikernels.aux import mesher
        import numpy as np
    except ImportError as exc:
        print(
            f"ERROR: axikernels is not installed or cannot be imported:\n"
            f"  {exc}\n"
            f"Install it with:\n"
            f"  pip install -e /path/to/AxiSEM3D_Kernels",
            file=sys.stderr,
        )
        sys.exit(1)

    # ── load forward simulation ───────────────────────────────────────────────
    print(f"Loading forward element output from: {fwd_elements}")
    forward_sim = eo.ElementOutput(path_to_element_output=fwd_elements)

    # ── build receiver point [radius_m, lat_deg, lon_deg] ────────────────────
    receiver_lat, receiver_lon = args.receiver
    receiver_point = [6_371_000.0, receiver_lat, receiver_lon]

    # ── run backward simulation (auto-answers all prompts) ────────────────────
    print("Computing adjoint source and launching backward simulation...")
    print(f"  tau         = {args.tau} s")
    print(f"  receiver    = lat={receiver_lat}, lon={receiver_lon}")
    print(f"  window      = {args.window[0]}–{args.window[1]} s")
    print(f"  channel     = {args.channel}")
    print(f"  MPI ranks   = {args.cores}")

    original_input = _patch_interactive()
    try:
        time_shift_obj = of.XObjectiveFunction(forward_sim)
        time_shift_obj.compute_backward_field(
            tau=args.tau,
            receiver_point=receiver_point,
            window=args.window,
            channel=args.channel,
            cores=args.cores,
        )
    finally:
        _restore_input(original_input)

    if time_shift_obj.backward_simulation is None:
        print("ERROR: backward simulation was not produced.", file=sys.stderr)
        sys.exit(1)

    print("Backward simulation complete.")

    # ── build kernel ──────────────────────────────────────────────────────────
    from axikernels.core.kernels import kernel as kernel_mod
    ker = kernel_mod.Kernel(forward_sim, time_shift_obj.backward_simulation)

    # ── build slice mesh ──────────────────────────────────────────────────────
    source_location = [
        forward_sim.Domain_Radius - forward_sim.source_depth,
        forward_sim.source_lat,
        forward_sim.source_lon,
    ]
    domains = [[4_000_000.0, 6_371_000.0, -np.pi, np.pi]]
    print(f"Building slice mesh (resolution={args.resolution})...")
    slc = mesher.SliceMesh(
        point1=source_location,
        point2=receiver_point,
        domains=domains,
        resolution=args.resolution,
    )

    # ── evaluate vp kernel ────────────────────────────────────────────────────
    print("Evaluating vp sensitivity kernel (this may take a few minutes)...")
    kernel_values = ker.evaluate_vp(slc.points)

    # ── save outputs ──────────────────────────────────────────────────────────
    h5_path = os.path.join(args.output, "vp_kernel")
    png_path = os.path.join(args.output, "vp_kernel")

    print(f"Saving HDF5 to {h5_path}.h5 ...")
    slc.save_data(data=kernel_values, filename=h5_path)

    print(f"Saving quick-look figure to {png_path}.png ...")
    slc.plot_on_mesh(data=kernel_values, high_range=0.1, filename=png_path)

    print("")
    print("Done.  Outputs:")
    print(f"  {h5_path}.h5")
    print(f"  {png_path}.png")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        import traceback
        traceback.print_exc()
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)

