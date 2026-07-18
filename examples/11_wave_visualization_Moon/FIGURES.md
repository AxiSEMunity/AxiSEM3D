# Reproducing the Moon-impact figures

Self-contained Python scripts that turn the impact-wavefield data into the
publication figures and movies. Each reads only the local `data/`, `assets/`, and
`pub_style.py` in this folder — nothing is hard-wired to a particular machine.

## Requirements

```
pip install numpy matplotlib scipy cartopy pillow pyvista imageio-ffmpeg
```

- **numpy, matplotlib, scipy** — core.
- **cartopy** — map projections for the 2-D globe figures.
- **Pillow** — reads the lunar-texture JPEG.
- **pyvista + vtk** — the true 3-D textured globe; off-screen/headless rendering works
  out of the box (`pv.OFF_SCREEN = True`).
- **imageio-ffmpeg** — writes MP4s; it bundles its own ffmpeg binary, so **no system
  ffmpeg is needed**.

## Inputs (already in this folder)

- `data/moon_surface.npz` — the impact surface field `U(Δ, t)` sampled along one
  meridian (`dist`, `time`, `uz`) plus one 3-component seismogram. This is the compact
  extract the plots consume; regenerate it from a full AxiSEM3D run if you change the
  simulation. *(The high-res `moon_8k.jpg` texture ships only in the example_gallery
  copy; the slim copy has `moon_2k.jpg`.)*
- `assets/moon_2k.jpg` (+ `moon_8k.jpg` in the gallery) — real lunar imagery.
- `pub_style.py` — shared plotting style (Helvetica, colour-blind-safe palette).

## Scripts

| Script | Output | What it shows |
|---|---|---|
| `make_figures.py` | record section, seismogram, overview, plain snapshots | the base publication set |
| `make_surface_snapshots.py` | `moon_surface_snapshots.pdf` | 8 snapshots of the wavefield on the **real Moon surface** (2-D orthographic) |
| `make_topo_overlay.py [--movie]` | `moon_wavefield_on_surface.png` / `.mp4` | wavefront over the real surface — still, or an MP4 |
| `make_movie.py` | `movies/moon_wavefield.mp4` | diverging-colormap globe animation |
| `make_3d_globe.py [--spin|--time]` | `moon_3d_globe.png` / `movies/moon_3d_*.mp4` | wavefield draped on a **true 3-D textured globe** (still, spinning, or time-evolving) |
| `make_3d_snapshots.py` | `moon_3d_snapshots.pdf` | 8 snapshots on the 3-D globe (smooth/axisymmetric impact run) |
| `make_3d_snapshots_scatter.py` | `moon_3d_snapshots_scatter.pdf` | 8 snapshots on the 3-D globe, for the **megaregolith-scattering** run — genuinely 2-D `U(lat, lon, t)` on a receiver grid, not a single meridian |

Run any of them with `python make_<name>.py`. Movies are written to a local `movies/`
folder. The 3-D scripts prefer `assets/moon_8k.jpg` and fall back to `moon_2k.jpg`.
All figures embed **editable Helvetica** text (Illustrator-friendly).

## Scattering run

`make_3d_snapshots_scatter.py` needs the megaregolith-scattering simulation's own
input/output, which is too large to ship in the slim `examples/` copy (the mesh
alone is ~39 MB): the **example_gallery** copy of this example ships
`sim/input_scatter/` (config, receiver grid, and the stochastic scattering model
`moon_scatter.nc`, ~57 MB total) so you can run the solver yourself —
```
axisem3d --input sim/input_scatter --output sim/output_scatter
```
— then `python make_3d_snapshots_scatter.py`. `moon_scatter.nc` holds fractional
vp/vs/rho perturbations for a stochastic megaregolith layer (top few tens of km,
~8% RMS, Gaussian-correlated) injected as a `StructuredGridV3D` 3-D model on top of
the 1-D `ISSI_MOON_M1` reference model; the azimuthal scattering this produces is
why the simulation needs `type_Nr: CONSTANT` at a high value (256) rather than the
usual `POINTWISE`/depth-adaptive Nr used by the smooth run.
