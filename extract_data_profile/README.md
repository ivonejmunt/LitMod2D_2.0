# Spatial profile extractor

`extract_profile.py` extracts values along a profile defined by two longitude-latitude endpoints.

## Inputs

- Geographic grids: `.nc`, `.nc4`, `.cdf`, `.grd`
- Point tables: `.xyz`, `.dat`, `.txt`, `.csv`, `.tsv`, `.xls`, `.xlsx`

Grid coordinates must be longitude and latitude in degrees. Headerless tables use the first three columns as longitude, latitude and value. Named columns can be selected with `--lon-col`, `--lat-col` and `--value-col`.

Install the Python packages with:

```powershell
pip install -r requirements.txt
```

GMT is used to convert legacy binary `.grd` files. NetCDF grids and point tables do not use GMT.

## Methods

`corridor` is the default method. It projects source points onto the profile, keeps points within the cross-track corridor and calculates a distance-weighted moving mean. The defaults are:

- 5 km profile sampling
- 20 km cross-track limit
- ±20 km along-profile moving window

`value_std` is the population standard deviation of the values inside each moving window.

`centerline` directly interpolates a regular grid along the profile.

## Run

Interactive mode:

```powershell
python extract_profile.py
```

The corridor plot prompt is:

```text
Spatial variability display [1=shaded band, 2=short SD bars] [1]:
```

Run with endpoints on the command line:

```powershell
python extract_profile.py input.nc 81.5 40.5 90 48.7
```

Draw the spatial SD with short bars:

```powershell
python extract_profile.py input.nc 81.5 40.5 90 48.7 --spread-style bars
```

Read a named value column from a table:

```powershell
python extract_profile.py observations.csv 81.5 40.5 90 48.7 --lon-col Longitude --lat-col Latitude --value-col HeatFlow
```

Select a variable and depth from a multidimensional grid:

```powershell
python extract_profile.py model.nc 81.5 40.5 90 48.7 --variable dVp --select depth=100
```

Outputs are written to `profile_output` unless `--output-dir` is supplied:

- profile CSV with all calculated columns;
- headerless, space-delimited LitMod `.dat` file;
- PNG and PDF plots.

The corridor-mode `.dat` file follows the observed-profile format used by LitMod:

```text
distance_km  value_mean  value_std
```

Centre-line mode writes `distance_km` and the interpolated value because it has no moving-window spatial SD.

Use `--save-details` to also save the projected source points and a JSON file containing the input hash and calculation settings:

```powershell
python extract_profile.py input.nc 81.5 40.5 90 48.7 --save-details
```

Run `python extract_profile.py --help` for all options.
