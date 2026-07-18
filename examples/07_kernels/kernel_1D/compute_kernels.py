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
import os
import subprocess
import sys

# ── set non-interactive matplotlib backend before any other import ────────────
import matplotlib
matplotlib.use("Agg")


def _setup_mpi_env(binary_path):
    """Auto-detect the MPI runtime from the binary's linked libmpi and
    prepend the matching mpirun directory to PATH and LD_LIBRARY_PATH.

    Also adds the binary's directory to PATH so that ``mpirun -np N axisem3d``
    can locate the executable (mpirun resolves via PATH, not cwd).

    This mirrors the auto-detection logic in run.sh so that axikernels'
    internal ``mpirun -np N axisem3d`` call uses the correct MPI runtime.
    """
    # Make the binary findable by mpirun via PATH
    binary_dir = os.path.dirname(os.path.abspath(binary_path))
    os.environ["PATH"] = (
        binary_dir + os.pathsep + os.environ.get("PATH", "")
    )

    try:
        result = subprocess.run(
            ["ldd", binary_path],
            capture_output=True, text=True, timeout=10,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return

    for line in result.stdout.splitlines():
        if "libmpi" in line and "=>" in line:
            parts = line.split("=>")
            if len(parts) < 2:
                continue
            lib_path = parts[1].strip().split()[0]
            if not lib_path or lib_path == "not":
                continue
            conda_lib = os.path.dirname(lib_path)
            candidate = os.path.join(os.path.dirname(conda_lib), "bin", "mpirun")
            if os.path.isfile(candidate):
                os.environ["PATH"] = (
                    os.path.dirname(candidate)
                    + os.pathsep
                    + os.environ.get("PATH", "")
                )
                os.environ["LD_LIBRARY_PATH"] = (
                    conda_lib
                    + os.pathsep
                    + os.environ.get("LD_LIBRARY_PATH", "")
                )
                print(f"  MPI auto-detected from binary:")
                print(f"    mpirun: {candidate}")
                print(f"    LD_LIBRARY_PATH prepend: {conda_lib}")
            break


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
    p.add_argument("--cores", type=int, default=8,
                   help="MPI ranks for backward simulation (default: 8)")
    p.add_argument("--resolution", type=int, default=200,
                   help="Kernel slice grid resolution (default: 200)")
    return p.parse_args()


def main():
    args = parse_args()
    forward_dir = os.path.normpath(args.forward)

    # ── validate forward simulation path ──────────────────────────────────────
    fwd_elements = os.path.join(forward_dir, "output", "elements")
    fwd_input = os.path.join(forward_dir, "input")
    fwd_binary = os.path.join(forward_dir, "axisem3d")

    if not os.path.isdir(fwd_elements):
        print(
            f"ERROR: forward element output directory not found:\n"
            f"  {fwd_elements}\n"
            f"Run the forward simulation first:\n"
            f"  bash run.sh forward",
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.isdir(fwd_input):
        print(
            f"ERROR: forward simulation input directory not found:\n"
            f"  {fwd_input}\n"
            f"The --forward argument must point to the full simulation directory,\n"
            f"not only to the element output directory.",
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.isfile(fwd_binary):
        print(
            f"ERROR: axisem3d binary not found in the forward simulation directory:\n"
            f"  {fwd_binary}\n"
            f"Run the forward simulation with bash run.sh forward so axikernels can\n"
            f"reuse the copied solver binary for the backward run.",
            file=sys.stderr,
        )
        sys.exit(1)

    # ── auto-detect MPI runtime from the binary's linked libmpi ─────────────
    # The axikernels backward simulation launcher uses bare "mpirun" from PATH.
    # We must ensure PATH and LD_LIBRARY_PATH point to the same MPI runtime
    # the binary was compiled against (mirrors the logic in run.sh).
    _setup_mpi_env(fwd_binary)

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
    receiver_point = [forward_sim.Domain_Radius, receiver_lat, receiver_lon]

    backward_dir = os.path.join(
        os.path.dirname(forward_dir),
        f"backward_{os.path.basename(forward_dir)}",
    )
    bwd_elements = os.path.join(backward_dir, "output", "elements")

    time_shift_obj = of.XObjectiveFunction(forward_sim, interactive=False)

    # ── reuse or re-run backward simulation ──────────────────────────────────
    # Skip the (expensive) simulation if the backward element output already
    # exists. This lets you re-run kernel computation without re-running
    # AxiSEM3D. Delete the backward directory to force a fresh simulation.
    _bwd_has_output = os.path.isdir(bwd_elements) and any(
        f.endswith(".nc") or f.endswith(".nc.rank0")
        for root, _, files in os.walk(bwd_elements)
        for f in files
    )

    if _bwd_has_output:
        print(f"Backward simulation output found at: {bwd_elements}")
        print("  Skipping simulation — loading existing output.")
        time_shift_obj.backward_simulation = eo.ElementOutput(
            path_to_element_output=bwd_elements
        )
    else:
        print("Computing adjoint source and launching backward simulation...")
        print(f"  tau         = {args.tau} s")
        print(f"  receiver    = lat={receiver_lat}, lon={receiver_lon}")
        print(f"  window      = {args.window[0]}–{args.window[1]} s")
        print(f"  channel     = {args.channel}")
        print(f"  MPI ranks   = {args.cores}")

        time_shift_obj.compute_backward_field(
            tau=args.tau,
            receiver_point=receiver_point,
            window=args.window,
            channel=args.channel,
            cores=args.cores,
        )

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
    domains = [[4_000_000.0, forward_sim.Domain_Radius, -np.pi, np.pi]]
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

    # ─── Moho interface kernel ────────────────────────────────────────────────
    # PREM Moho radius verified from simu_forward/input/prem_iso_elastic.bm
    # (Discontinuity 2, depth 24.40 km = 6371000 - 24400 m)
    MOHO_RADIUS = 6_346_600.0

    import matplotlib.pyplot as plt
    import h5py

    # Build a regular lat/lon grid covering the source–receiver great-circle
    # path. Default: source at (0°, 0°) → receiver at (0°, 40°) along equator.
    moho_n_lat = max(args.resolution // 4, 25)
    moho_n_lon = max(args.resolution // 2, 50)
    moho_lat_deg = np.linspace(-15.0, 15.0, moho_n_lat)
    moho_lon_deg = np.linspace(-5.0, 45.0, moho_n_lon)
    moho_LON_deg, moho_LAT_deg = np.meshgrid(moho_lon_deg, moho_lat_deg)
    # evaluate_K_dn / evaluate_K_dv expect (N, 2) [lat, lon] in radians
    moho_points = np.column_stack(
        [np.deg2rad(moho_LAT_deg.ravel()),
         np.deg2rad(moho_LON_deg.ravel())]
    )

    print(f"\nBuilding Moho interface kernel on {len(moho_points)} shell points ...")
    print(f"  radius     = {MOHO_RADIUS:.0f} m  "
          f"({(6_371_000 - MOHO_RADIUS) / 1000.0:.1f} km depth)")
    print(f"  lat range  = {moho_lat_deg[0]:.1f}°  –  {moho_lat_deg[-1]:.1f}°")
    print(f"  lon range  = {moho_lon_deg[0]:.1f}°  –  {moho_lon_deg[-1]:.1f}°")

    print("Evaluating K_dn component ...")
    kd_dn = ker.evaluate_K_dn(moho_points, MOHO_RADIUS)
    print("Evaluating K_dv component ...")
    kd_dv = ker.evaluate_K_dv(moho_points, MOHO_RADIUS)
    kd_total = kd_dn + kd_dv

    # ── diagnostics ───────────────────────────────────────────────────────────
    for name, arr in [("K_dn", kd_dn), ("K_dv", kd_dv), ("Kd", kd_total)]:
        n_nan = int(np.sum(~np.isfinite(arr)))
        n_zero = int(np.sum(arr == 0.0))
        print(f"  {name:5s}: max={np.nanmax(np.abs(arr)):.4e}, "
              f"zeros={n_zero}/{len(arr)}, non-finite={n_nan}/{len(arr)}")

    # ── save HDF5 (Kd, K_dn, K_dv + lat/lon grids) ───────────────────────────
    moho_h5_path = os.path.join(args.output, "moho_kd.h5")
    print(f"Saving HDF5 to {moho_h5_path} ...")
    with h5py.File(moho_h5_path, "w") as hf:
        hf.create_dataset("lat_deg",  data=moho_LAT_deg)
        hf.create_dataset("lon_deg",  data=moho_LON_deg)
        hf.create_dataset("Kd",       data=kd_total.reshape(moho_LAT_deg.shape))
        hf.create_dataset("K_dn",     data=kd_dn.reshape(moho_LAT_deg.shape))
        hf.create_dataset("K_dv",     data=kd_dv.reshape(moho_LAT_deg.shape))
        hf.attrs["moho_radius_m"] = MOHO_RADIUS
        hf.attrs["description"] = (
            "Moho topography sensitivity kernel. "
            "Kd = K_dn + K_dv (total interface kernel); "
            "K_dn = normal-displacement (traction-jump) part; "
            "K_dv = velocity-contrast (material-jump) part."
        )

    # ── save PNG ──────────────────────────────────────────────────────────────
    moho_png_path = os.path.join(args.output, "moho_kd.png")
    print(f"Saving figure to {moho_png_path} ...")
    kd_grid = kd_total.reshape(moho_LAT_deg.shape)
    # Use 99th percentile of finite, non-zero absolute values for color scale
    finite_mask = np.isfinite(kd_grid) & (kd_grid != 0.0)
    abs_vals = np.abs(kd_grid[finite_mask]) if np.any(finite_mask) \
        else np.array([1.0])
    cbar_max = float(np.percentile(abs_vals, 99))

    fig, ax = plt.subplots(figsize=(8, 4))
    pcm = ax.pcolormesh(
        moho_LON_deg, moho_LAT_deg, kd_grid,
        cmap="RdBu_r", vmin=-cbar_max, vmax=cbar_max,
        shading="auto",
    )
    # Mark source and receiver
    src_lat_deg = float(np.rad2deg(source_location[1]))
    src_lon_deg = float(np.rad2deg(source_location[2]))
    ax.scatter([src_lon_deg], [src_lat_deg],
               c="white", s=80, marker="*", zorder=5, label="Source")
    ax.scatter([receiver_lon], [receiver_lat],
               c="cyan",  s=80, marker="^", zorder=5, label="Receiver")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar_ticks = np.linspace(-cbar_max, cbar_max, 5)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{t:.2e}" for t in cbar_ticks])
    cbar.set_label("Kd  (interface sensitivity)")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(
        f"Moho topography kernel (Kd),  r = {MOHO_RADIUS / 1000.0:.1f} km"
        f"  ({(6_371_000 - MOHO_RADIUS) / 1000.0:.1f} km depth)"
    )
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(moho_png_path, dpi=150)
    plt.close(fig)

    print("")
    print("Done.  Outputs:")
    print(f"  {moho_h5_path}")
    print(f"  {moho_png_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        import traceback
        traceback.print_exc()
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)

