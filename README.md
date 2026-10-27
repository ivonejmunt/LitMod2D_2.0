# LitMod2D_2.0 — Python 3 / PyQt5 Update

This repository provides an updated version of **LitMod2D_2.0**, a geophysical–petrological tool for integrated modeling of the thermal, compositional, density, and seismic structure of the lithosphere and sublithospheric upper mantle along two-dimensional profiles.

This update modernizes the Python interface for **Python 3 and PyQt5** and adapts the related GUI dependencies and workflows. The original **Fortran numerical solver is retained without changes**.

## Included components

- Python 3 / PyQt5 interface for building, loading, editing, running, and visualizing 2D lithospheric models.
- The original LitMod2D Fortran numerical solver, using platform-specific executables for Windows, Linux, and macOS.
- Runtime libraries and interface assets required by the corresponding platform.
- Python dependency list (`requirements.txt`).

**Note:** Thermodynamic lookup tables, study-specific input files, and other data required for numerical runs should be provided according to the original LitMod2D workflow.

## Installation

1. Install a supported Python 3 environment, then install dependencies:

   ```bash
   python -m pip install -r requirements.txt
   ```

2. Start the application from the `GUI` directory so the existing relative paths for images, inputs, and the solver work as expected:

   ```bash
   cd GUI
   python main.py
   ```

3. Use **Build Model** or **Load Model** in the interface. For numerical runs, ensure that the required model inputs and thermodynamic data tables are available in the working directory and retain their original LitMod2D formats.

Place the solver executable for your operating system in the `GUI/` directory, using the filenames expected by the existing application: `LITMOD_V4_WINDOWS.exe` (Windows), `LITMOD_V4_LINUX` (Linux), or `LITMOD_V4_MAC` (macOS). On Linux and macOS, grant execution permission if needed (`chmod +x`); on Windows, keep any required runtime DLLs alongside the executable.

## Repository structure

```text
LitMod2D_2.0-Python3-PyQt5/
├── GUI/
│   ├── main.py
│   ├── my_lib.py
│   ├── plot_lib.py
│   ├── ...                         # Other GUI and model editing modules
│   ├── Images/
│   ├── LITMOD_V4_*                 # Platform-specific Fortran solvers
│   └── *.dll                       # Platform-specific runtime libraries, if required
├── requirements.txt
├── README.md
├── LICENSE
└── .gitignore
```

## Scientific basis and citation

This version is based on **LitMod2D_2.0** developed by the original authors. The scientific methodology and numerical solver remain attributable to the original research team. Please cite the relevant original publications when using the software:

- Afonso, J. C., Fernández, M., Ranalli, G., Griffin, W. L., & Connolly, J. A. D. (2008). Integrated geophysical-petrological modeling of the lithosphere and sublithospheric upper mantle: Methodology and applications. *Geochemistry, Geophysics, Geosystems*. https://doi.org/10.1029/2007GC001834
- Kumar, A. et al. (2020). LitMod2D_2.0. *Geochemistry, Geophysics, Geosystems*. https://doi.org/10.1029/2019GC008777

Upstream distribution: https://github.com/ajay6763/LitMod2D_2.0_package_dist_users
