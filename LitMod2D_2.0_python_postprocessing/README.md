# LitMod2D_2.0 Python post-processing

Three scripts for plotting completed LitMod2D profile models.

| Script | Figure |
|---|---|
| `001_plot_model_fit_and_section.py` | Observed/calculated elevation, surface heat flow, geoid, Bouguer anomaly, and temperature or density section |
| `002_plot_density_velocity_temperature.py` | Density section, temperature section, and density/Vp/Vs/temperature columns |
| `003_plot_vp_vs_and_anomalies.py` | Topography, Vp, Vp anomaly, Vs, and Vs anomaly |

Script 003 uses `support/ak135.txt` for its ak135 reference. Keep the `support`
directory with the scripts.

## Installation

Create a conda environment with Python, GMT, and the plotting libraries:

```bash
conda create -n litmod-post -c conda-forge python=3.12 numpy pandas scipy matplotlib pygmt
conda activate litmod-post
```

Check the GMT/PyGMT installation:

```bash
python -c "import pygmt; pygmt.show_versions()"
```

`requirements.txt` lists the Python packages for pip installations. GMT is a
separate library; scripts 001 and 003 call it through PyGMT. Script 002 uses
Matplotlib and SciPy.

## Running the scripts

Interactive mode asks for the model directory and the figure-specific choice:

```bash
python 001_plot_model_fit_and_section.py
python 002_plot_density_velocity_temperature.py
python 003_plot_vp_vs_and_anomalies.py
```

The same runs can be written as commands:

```bash
python 001_plot_model_fit_and_section.py MODEL_DIR --style new --plot-mode both
python 002_plot_density_velocity_temperature.py MODEL_DIR --references 200,500,910
python 003_plot_vp_vs_and_anomalies.py MODEL_DIR --reference-x 910
```

`MODEL_DIR` is the completed LitMod2D model directory. Each script also
provides `--help`.

## Input files

### 001: model fit and sections

Standard filenames:

- observations: `0Topo.dat`, `0HF.dat`, `0Geoid.dat`, `0Med_CBA.dat`;
- calculated profiles: `topo_out.dat`, `SHF_out.dat`, `geoid_out.dat`,
  `bouguer_out.dat`;
- sections and boundaries: `tempout.dat`, `dens.dat`, `bodies.out`.

Observation tables contain `distance_km value [uncertainty]`. Calculated
profiles contain `distance_km value`. Section files contain
`distance_km vertical_coordinate_km value`. Body boundaries use GMT segments
beginning with `>`.

`FILE_CANDIDATES` in script 001 lists the accepted filename aliases. Matching
is case-insensitive. `--style new` draws uncertainty bands; `--style old`
draws points and error bars. `--plot-mode both` writes separate temperature
and density figures.

### 002: density, temperature, Vp, and Vs

Required files: `dens.dat`, `velocities.dat`, and `tempout.dat`.
`bodies.out` supplies the boundary lines. `--references` takes comma-separated
profile distances in kilometres. Interactive mode proposes the three interior
quarter points.

`002_density_overview` uses density for the 2-D section.
`002_temperature_velocity_overview` uses temperature for the 2-D section and
plots density, Vp, Vs, and temperature at the selected reference columns.

### 003: Vp, Vs, and anomalies

`velocities.dat` contains `X Z Vp Vs`. A six-column
`post_processing_*.dat` file is also accepted, with Vp and Vs in columns 5 and
6. Boundary lines are read from `bodies_GUI.out` or `bodies.out`.

The anomaly calculation is:

```text
(model value - reference value) / reference value * 100%
```

Reference choices:

- interactive option 1: `support/ak135.txt`;
- `--reference-model FILE`: a table containing depth, Vp, and Vs;
- `--reference-x DISTANCE_KM`: the nearest model column.

Record the selected reference with each anomaly figure.
The plotting boundary is extended to the next 10 km interval so the GMT region
matches the 10 km by 5 km anomaly grid.

## Outputs and units

Outputs are written to `MODEL_DIR/python_postprocessing` unless `--output` or
`--output-dir` is supplied. Repeating a run replaces figures with the same
name.

Distances and depths are in kilometres. Density is in kg/m3, temperature in
degrees Celsius, and Vp/Vs in km/s. The section grids are interpolated for
plotting.

## Systems checked

- Windows: all three scripts completed with a LitMod2D model.
- Ubuntu under WSL: all three scripts compiled with Python 3.12.
- Windows, macOS, and Linux packages for GMT/PyGMT are available through
  conda-forge.

