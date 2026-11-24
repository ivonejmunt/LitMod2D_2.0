"""Extract a profile from a geographic grid or longitude-latitude-value table.

The default corridor mode applies the distance-weighted moving window used in
the original PROJBANDA workflow. Centerline mode directly interpolates a grid.
Run ``python extract_profile.py --help`` for examples and options.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from math import atan2, cos, radians, sin
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


EARTH_RADIUS_KM = 6371.0088
GRID_SUFFIXES = {".nc", ".nc4", ".cdf", ".grd"}
TABLE_SUFFIXES = {".csv", ".tsv", ".txt", ".xyz", ".dat", ".xls", ".xlsx"}
LON_ALIASES = ("longitude", "lon", "long", "x")
LAT_ALIASES = ("latitude", "lat", "y")
VALUE_ALIASES = ("value", "z", "data", "elevation", "depth", "temperature", "density")


@dataclass(frozen=True)
class Profile:
    start_lon: float
    start_lat: float
    end_lon: float
    end_lat: float


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def haversine_distance_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    lon1r, lat1r, lon2r, lat2r = np.radians([lon1, lat1, lon2, lat2])
    dlon = lon2r - lon1r
    dlat = lat2r - lat1r
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2.0) ** 2
    return float(EARTH_RADIUS_KM * 2.0 * np.arctan2(np.sqrt(a), np.sqrt(max(0.0, 1.0 - a))))


def initial_bearing_rad(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    lon1r, lat1r, lon2r, lat2r = map(radians, (lon1, lat1, lon2, lat2))
    dlon = lon2r - lon1r
    y = sin(dlon) * cos(lat2r)
    x = cos(lat1r) * sin(lat2r) - sin(lat1r) * cos(lat2r) * cos(dlon)
    return atan2(y, x)


def sample_distances(profile: Profile, step_km: float) -> np.ndarray:
    length = haversine_distance_km(
        profile.start_lon, profile.start_lat, profile.end_lon, profile.end_lat
    )
    if length <= 0.0:
        raise ValueError("The profile endpoints are identical")
    distances = np.arange(0.0, length, step_km, dtype=float)
    if distances.size == 0 or not np.isclose(distances[-1], length):
        distances = np.append(distances, length)
    return distances


def great_circle_track(profile: Profile, distances_km: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon1, lat1, lon2, lat2 = np.radians(
        [profile.start_lon, profile.start_lat, profile.end_lon, profile.end_lat]
    )
    start = np.array([np.cos(lat1) * np.cos(lon1), np.cos(lat1) * np.sin(lon1), np.sin(lat1)])
    end = np.array([np.cos(lat2) * np.cos(lon2), np.cos(lat2) * np.sin(lon2), np.sin(lat2)])
    omega = float(np.arccos(np.clip(np.dot(start, end), -1.0, 1.0)))
    if omega == 0.0:
        raise ValueError("The profile endpoints are identical")
    fractions = np.clip(distances_km / (omega * EARTH_RADIUS_KM), 0.0, 1.0)
    vectors = (
        np.sin((1.0 - fractions) * omega)[:, None] * start[None, :]
        + np.sin(fractions * omega)[:, None] * end[None, :]
    ) / np.sin(omega)
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]
    lon = np.degrees(np.arctan2(vectors[:, 1], vectors[:, 0]))
    lat = np.degrees(np.arctan2(vectors[:, 2], np.hypot(vectors[:, 0], vectors[:, 1])))
    return lon, lat


def project_arrays(profile: Profile, lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon1 = radians(profile.start_lon)
    lat1 = radians(profile.start_lat)
    lonr = np.radians(lon.astype(float))
    latr = np.radians(lat.astype(float))
    dlon = lonr - lon1
    dlat = latr - lat1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(latr) * np.sin(dlon / 2.0) ** 2
    angular = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(np.maximum(0.0, 1.0 - a)))
    bearing = np.arctan2(
        np.sin(dlon) * np.cos(latr),
        np.cos(lat1) * np.sin(latr) - np.sin(lat1) * np.cos(latr) * np.cos(dlon),
    )
    profile_bearing = initial_bearing_rad(
        profile.start_lon, profile.start_lat, profile.end_lon, profile.end_lat
    )
    delta = bearing - profile_bearing
    cross = np.arcsin(np.clip(np.sin(angular) * np.sin(delta), -1.0, 1.0)) * EARTH_RADIUS_KM
    along = np.arctan2(np.sin(angular) * np.cos(delta), np.cos(angular)) * EARTH_RADIUS_KM
    return along, cross


def parse_selectors(items: Sequence[str] | None) -> dict[str, str | float]:
    selectors: dict[str, str | float] = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"Invalid selector {item!r}; use DIM=VALUE")
        name, raw = (part.strip() for part in item.split("=", 1))
        if not name or not raw:
            raise ValueError(f"Invalid selector {item!r}; use DIM=VALUE")
        try:
            value: str | float = float(raw)
        except ValueError:
            value = raw
        selectors[name] = value
    return selectors


def choose_name(names: Sequence[str], explicit: str | None, aliases: Sequence[str], kind: str) -> str:
    lookup = {str(name).lower(): str(name) for name in names}
    if explicit:
        if explicit in names:
            return explicit
        if explicit.lower() in lookup:
            return lookup[explicit.lower()]
        raise ValueError(f"{kind} coordinate {explicit!r} not found; available: {list(names)}")
    for alias in aliases:
        if alias in lookup:
            return lookup[alias]
    raise ValueError(f"Could not identify the {kind} coordinate; available: {list(names)}")


def open_grid(
    path: Path,
    variable: str | None,
    lon_name: str | None,
    lat_name: str | None,
    selectors: dict[str, str | float],
) -> tuple[xr.DataArray, str, str, str]:
    dataset: xr.Dataset | None = None
    temporary_grid: tempfile.TemporaryDirectory[str] | None = None
    try:
        try:
            dataset = xr.open_dataset(path)
        except (OSError, ValueError) as xarray_error:
            if path.suffix.lower() != ".grd":
                raise
            # Legacy GMT binary grids are converted to temporary NetCDF.
            gmt = shutil.which("gmt")
            if gmt is None:
                raise ValueError(
                    "This .grd file is not NetCDF. Run the script in a PyGMT environment "
                    "that provides the GMT executable, or convert the grid to NetCDF."
                ) from xarray_error
            temporary_grid = tempfile.TemporaryDirectory(prefix="profile_grid_")
            converted = Path(temporary_grid.name) / "converted.nc"
            process = subprocess.run(
                [gmt, "grdconvert", str(path), str(converted)],
                text=True,
                capture_output=True,
                check=False,
            )
            if process.returncode != 0:
                raise ValueError(
                    "GMT could not convert the legacy .grd file: "
                    + (process.stderr.strip() or process.stdout.strip())
                )
            dataset = xr.open_dataset(converted)
        if variable:
            if variable not in dataset.data_vars:
                raise ValueError(f"Variable {variable!r} not found; available: {list(dataset.data_vars)}")
            data = dataset[variable]
            variable_name = variable
        else:
            candidates = [name for name, value in dataset.data_vars.items() if value.ndim >= 2]
            if len(candidates) != 1:
                raise ValueError(
                    "Select a grid variable with --variable; available data variables: "
                    + ", ".join(dataset.data_vars)
                )
            variable_name = candidates[0]
            data = dataset[variable_name]
        for dimension, value in selectors.items():
            if dimension not in data.dims and dimension not in data.coords:
                raise ValueError(f"Selector {dimension!r} not found in {variable_name}; dimensions: {data.dims}")
            if isinstance(value, float):
                data = data.sel({dimension: value}, method="nearest")
            else:
                data = data.sel({dimension: value})
        data = data.squeeze(drop=True)
        names = list(dict.fromkeys([*data.dims, *data.coords]))
        longitude = choose_name(names, lon_name, LON_ALIASES, "longitude")
        latitude = choose_name(names, lat_name, LAT_ALIASES, "latitude")
        if data[longitude].ndim != 1 or data[latitude].ndim != 1:
            raise ValueError("The current grid reader requires one-dimensional longitude and latitude coordinates")
        extra = [dim for dim in data.dims if dim not in {longitude, latitude}]
        if extra:
            raise ValueError(
                f"Grid still has unsliced dimensions {extra}; select them with repeated --select DIM=VALUE"
            )
        data = data.transpose(latitude, longitude).sortby(latitude).sortby(longitude).load()
        lon_values = data[longitude].values.astype(float)
        lat_values = data[latitude].values.astype(float)
        if np.nanmin(lon_values) < -360 or np.nanmax(lon_values) > 360:
            raise ValueError("Longitude coordinates are outside a geographic degree range")
        if np.nanmin(lat_values) < -90 or np.nanmax(lat_values) > 90:
            raise ValueError("Latitude coordinates are outside -90 to 90 degrees")
        return data, longitude, latitude, variable_name
    finally:
        if dataset is not None:
            dataset.close()
        if temporary_grid is not None:
            temporary_grid.cleanup()


def subset_grid(
    data: xr.DataArray,
    lon_name: str,
    lat_name: str,
    profile: Profile,
    buffer_km: float,
) -> xr.DataArray:
    latitude_limit = max(abs(profile.start_lat), abs(profile.end_lat))
    lat_buffer = buffer_km / 110.574
    lon_scale = max(1e-6, 111.320 * cos(radians(min(89.0, latitude_limit + lat_buffer))))
    lon_buffer = buffer_km / lon_scale
    lon_min = min(profile.start_lon, profile.end_lon) - lon_buffer
    lon_max = max(profile.start_lon, profile.end_lon) + lon_buffer
    lat_min = min(profile.start_lat, profile.end_lat) - lat_buffer
    lat_max = max(profile.start_lat, profile.end_lat) + lat_buffer
    subset = data.sel({lon_name: slice(lon_min, lon_max), lat_name: slice(lat_min, lat_max)})
    if subset.size == 0:
        raise ValueError("The profile and corridor do not overlap the grid")
    return subset


def grid_to_points(
    data: xr.DataArray,
    lon_name: str,
    lat_name: str,
    profile: Profile | None = None,
    max_perp_km: float | None = None,
    along_buffer_km: float = 0.0,
) -> pd.DataFrame:
    lon = data[lon_name].values.astype(float)
    lat = data[lat_name].values.astype(float)
    values = data.values.astype(float)
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    valid = np.isfinite(values)
    if profile is not None and max_perp_km is not None:
        along, cross = project_arrays(profile, lon_grid[valid], lat_grid[valid])
        length = haversine_distance_km(
            profile.start_lon, profile.start_lat, profile.end_lon, profile.end_lat
        )
        corridor = (
            (along >= -along_buffer_km)
            & (along <= length + along_buffer_km)
            & (np.abs(cross) <= max_perp_km)
        )
        valid_indices = np.flatnonzero(valid)
        filtered = np.zeros(valid.size, dtype=bool)
        filtered[valid_indices[corridor]] = True
        valid = filtered.reshape(valid.shape)
    return pd.DataFrame(
        {
            "longitude": lon_grid[valid],
            "latitude": lat_grid[valid],
            "value": values[valid],
            "point_id": np.char.add("grid_", np.arange(valid.sum()).astype(str)),
        }
    )


def first_data_line(path: Path) -> str:
    with path.open("r", encoding="utf-8-sig", errors="replace") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                return stripped
    raise ValueError(f"No data rows found in {path}")


def line_looks_numeric(line: str, delimiter: str | None) -> bool:
    fields = line.split(delimiter) if delimiter else line.replace(",", " ").split()
    if len(fields) < 3:
        return False
    try:
        [float(value) for value in fields[:3]]
    except ValueError:
        return False
    return True


def resolve_column(
    table: pd.DataFrame,
    specification: str | None,
    aliases: Sequence[str],
    fallback_index: int,
    kind: str,
) -> object:
    columns = list(table.columns)
    if specification:
        if specification in columns:
            return specification
        if specification.isdigit():
            index = int(specification) - 1
            if 0 <= index < len(columns):
                return columns[index]
        lookup = {str(column).lower(): column for column in columns}
        if specification.lower() in lookup:
            return lookup[specification.lower()]
        raise ValueError(f"{kind} column {specification!r} not found; available: {columns}")
    lookup = {str(column).strip().lower(): column for column in columns}
    for alias in aliases:
        if alias in lookup:
            return lookup[alias]
    if fallback_index < len(columns):
        return columns[fallback_index]
    raise ValueError(f"Could not identify the {kind} column; available: {columns}")


def read_table_points(
    path: Path,
    header: str,
    delimiter: str | None,
    sheet_name: str,
    lon_col: str | None,
    lat_col: str | None,
    value_col: str | None,
) -> tuple[pd.DataFrame, dict[str, str]]:
    actual_delimiter = delimiter.encode().decode("unicode_escape") if delimiter else None
    if path.suffix.lower() in {".xls", ".xlsx"}:
        header_row = None if header == "no" else 0
        sheet: str | int = int(sheet_name) if sheet_name.isdigit() else sheet_name
        table = pd.read_excel(path, sheet_name=sheet, header=header_row)
    else:
        if header == "auto":
            header_row = None if line_looks_numeric(first_data_line(path), actual_delimiter) else 0
        else:
            header_row = 0 if header == "yes" else None
        if actual_delimiter:
            table = pd.read_csv(
                path, sep=actual_delimiter, comment="#", header=header_row, encoding="utf-8-sig"
            )
        else:
            table = pd.read_csv(
                path,
                sep=None,
                engine="python",
                comment="#",
                header=header_row,
                encoding="utf-8-sig",
            )
    longitude = resolve_column(table, lon_col, LON_ALIASES, 0, "longitude")
    latitude = resolve_column(table, lat_col, LAT_ALIASES, 1, "latitude")
    value = resolve_column(table, value_col, VALUE_ALIASES, 2, "value")
    points = pd.DataFrame(
        {
            "longitude": pd.to_numeric(table[longitude], errors="coerce"),
            "latitude": pd.to_numeric(table[latitude], errors="coerce"),
            "value": pd.to_numeric(table[value], errors="coerce"),
            "point_id": [f"row_{index + 1}" for index in range(len(table))],
        }
    ).dropna(subset=["longitude", "latitude", "value"])
    if points.empty:
        raise ValueError("The selected table columns contain no valid numeric points")
    return points.reset_index(drop=True), {
        "longitude_column": str(longitude),
        "latitude_column": str(latitude),
        "value_column": str(value),
    }


def corridor_profile(
    points: pd.DataFrame,
    profile: Profile,
    sample_step_km: float,
    max_perp_km: float,
    half_width_km: float,
    coefficient: float,
    minimum_weight: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    along, cross = project_arrays(
        profile, points["longitude"].to_numpy(), points["latitude"].to_numpy()
    )
    length = haversine_distance_km(
        profile.start_lon, profile.start_lat, profile.end_lon, profile.end_lat
    )
    keep = (
        (along >= -half_width_km)
        & (along <= length + half_width_km)
        & (np.abs(cross) <= max_perp_km)
    )
    projected = points.loc[keep].copy().reset_index(drop=True)
    projected["distance_km"] = along[keep]
    projected["perp_km"] = cross[keep]
    projected["abs_perp_km"] = np.abs(cross[keep])
    distances = sample_distances(profile, sample_step_km)
    track_lon, track_lat = great_circle_track(profile, distances)
    rows: list[dict[str, float | int]] = []
    along_values = projected["distance_km"].to_numpy(dtype=float)
    cross_values = projected["perp_km"].to_numpy(dtype=float)
    source_values = projected["value"].to_numpy(dtype=float)
    for index, sample in enumerate(distances):
        selected = np.abs(along_values - sample) <= half_width_km
        count = int(selected.sum())
        if count == 0:
            rows.append(
                {
                    "distance_km": sample,
                    "longitude": track_lon[index],
                    "latitude": track_lat[index],
                    "value_mean": np.nan,
                    "value_std": np.nan,
                    "minimum": np.nan,
                    "maximum": np.nan,
                    "n_points": 0,
                }
            )
            continue
        selected_along = along_values[selected]
        selected_cross = cross_values[selected]
        values = source_values[selected]
        radial_squared = (selected_along - sample) ** 2 + selected_cross**2
        # Original PROJBANDA radial weight.
        weights = np.maximum(1.0 - coefficient * radial_squared, minimum_weight)
        mean = float(np.average(values, weights=weights))
        rows.append(
            {
                "distance_km": sample,
                "longitude": track_lon[index],
                "latitude": track_lat[index],
                "value_mean": mean,
                "value_std": float(np.std(values, ddof=0)),
                "minimum": float(np.min(values)),
                "maximum": float(np.max(values)),
                "n_points": count,
            }
        )
    return pd.DataFrame(rows), projected


def centerline_profile(
    data: xr.DataArray,
    lon_name: str,
    lat_name: str,
    profile: Profile,
    sample_step_km: float,
    method: str,
) -> pd.DataFrame:
    distances = sample_distances(profile, sample_step_km)
    track_lon, track_lat = great_circle_track(profile, distances)
    query_lon = xr.DataArray(track_lon, dims="sample")
    query_lat = xr.DataArray(track_lat, dims="sample")
    sampled = data.interp({lon_name: query_lon, lat_name: query_lat}, method=method).values.astype(float)
    return pd.DataFrame(
        {
            "distance_km": distances,
            "longitude": track_lon,
            "latitude": track_lat,
            "value_mean": sampled,
            "value_std": np.nan,
            "minimum": sampled,
            "maximum": sampled,
            "n_points": np.where(np.isfinite(sampled), 1, 0),
        }
    )


def write_litmod_dat(result: pd.DataFrame, mode: str, path: Path) -> list[str]:
    columns = ["distance_km", "value_mean"]
    if mode == "corridor":
        columns.append("value_std")
    rows = result.loc[result[columns].notna().all(axis=1), columns]
    np.savetxt(path, rows.to_numpy(dtype=float), fmt="%.4f")
    return columns


def plot_result(
    result: pd.DataFrame,
    label: str,
    mode: str,
    spread_style: str,
    output_stem: Path,
) -> list[Path]:
    valid = result["value_mean"].notna()
    fig, ax = plt.subplots(figsize=(11.0, 4.5), constrained_layout=True)
    if valid.any():
        x = result.loc[valid, "distance_km"].to_numpy()
        mean = result.loc[valid, "value_mean"].to_numpy()
        spread = result.loc[valid, "value_std"].to_numpy()
        if mode == "corridor" and np.isfinite(spread).any():
            if spread_style == "band":
                ax.fill_between(
                    x,
                    mean - spread,
                    mean + spread,
                    color="#4c78a8",
                    alpha=0.20,
                    label="±1 spatial SD",
                )
            else:
                ax.errorbar(
                    x,
                    mean,
                    yerr=spread,
                    fmt="none",
                    ecolor="#4c78a8",
                    elinewidth=0.7,
                    capsize=2.5,
                    capthick=0.7,
                    alpha=0.75,
                    label="±1 spatial SD",
                )
        ax.plot(x, mean, color="#1f4e79", linewidth=1.5, marker="o", markersize=2.5, label="Profile value")
    ax.set(xlabel="Distance along profile (km)", ylabel=label, title=f"{label} profile ({mode})")
    ax.grid(color="0.88", linewidth=0.5)
    ax.legend(frameon=False)
    png = output_stem.with_suffix(".png")
    pdf = output_stem.with_suffix(".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return [png, pdf]


def detect_input_type(path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    suffix = path.suffix.lower()
    if suffix in GRID_SUFFIXES:
        return "grid"
    if suffix in TABLE_SUFFIXES:
        return "table"
    raise ValueError(f"Cannot infer the input type from {suffix!r}; use --input-type")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_path", nargs="?", type=Path)
    parser.add_argument("start_lon", nargs="?", type=float)
    parser.add_argument("start_lat", nargs="?", type=float)
    parser.add_argument("end_lon", nargs="?", type=float)
    parser.add_argument("end_lat", nargs="?", type=float)
    parser.add_argument("--input-type", choices=("auto", "grid", "table"), default="auto")
    parser.add_argument("--mode", choices=("corridor", "centerline"), default="corridor")
    parser.add_argument("--sample-step-km", type=float, default=5.0)
    parser.add_argument("--max-perp-km", type=float, default=20.0)
    parser.add_argument("--window-half-width-km", type=float, default=20.0)
    parser.add_argument("--weight-coefficient", type=float, default=0.00036)
    parser.add_argument("--minimum-weight", type=float, default=0.1)
    parser.add_argument("--variable", help="Grid variable; required when a file has several 2-D variables.")
    parser.add_argument("--lon-name", help="Grid longitude coordinate name.")
    parser.add_argument("--lat-name", help="Grid latitude coordinate name.")
    parser.add_argument("--select", action="append", help="Select an extra grid dimension, for example depth=100.")
    parser.add_argument("--interpolation", choices=("linear", "nearest"), default="linear")
    parser.add_argument("--header", choices=("auto", "yes", "no"), default="auto")
    parser.add_argument("--delimiter", help=r"Table delimiter; use \t for a tab.")
    parser.add_argument("--sheet-name", default="0")
    parser.add_argument("--lon-col", help="Table longitude column name or 1-based number.")
    parser.add_argument("--lat-col", help="Table latitude column name or 1-based number.")
    parser.add_argument("--value-col", help="Table value column name or 1-based number.")
    parser.add_argument("--value-label", help="Y-axis label; default is the variable/value column name.")
    parser.add_argument(
        "--spread-style",
        choices=("band", "bars"),
        help="Show spatial SD as a shaded band or short capped bars; default is band.",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--output-prefix")
    parser.add_argument(
        "--save-details",
        action="store_true",
        help="Also save projected source points and calculation provenance.",
    )
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


def prompt_missing(args: argparse.Namespace) -> None:
    interactive = args.input_path is None or any(
        value is None for value in (args.start_lon, args.start_lat, args.end_lon, args.end_lat)
    )
    if args.input_path is None:
        args.input_path = Path(input("Grid, NetCDF, XYZ, CSV or Excel file: ").strip().strip('"'))
    coordinates = [args.start_lon, args.start_lat, args.end_lon, args.end_lat]
    if any(value is None for value in coordinates):
        start = input("Profile start longitude latitude: ").replace(",", " ").split()
        end = input("Profile end longitude latitude: ").replace(",", " ").split()
        if len(start) != 2 or len(end) != 2:
            raise ValueError("Enter two numbers for each endpoint")
        args.start_lon, args.start_lat = map(float, start)
        args.end_lon, args.end_lat = map(float, end)
    if args.spread_style is None:
        args.spread_style = "band"
        if interactive and args.mode == "corridor" and not args.no_plot:
            choice = input(
                "Spatial variability display [1=shaded band, 2=short SD bars] [1]: "
            ).strip()
            if choice in ("", "1"):
                args.spread_style = "band"
            elif choice == "2":
                args.spread_style = "bars"
            else:
                raise ValueError("Enter 1 for a shaded band or 2 for short SD bars")


def validate_args(args: argparse.Namespace) -> None:
    if not args.input_path.exists():
        raise FileNotFoundError(args.input_path)
    if args.sample_step_km <= 0 or args.max_perp_km <= 0 or args.window_half_width_km <= 0:
        raise ValueError("Sampling and distance parameters must be positive")
    if args.weight_coefficient < 0:
        raise ValueError("--weight-coefficient cannot be negative")
    if not 0 <= args.minimum_weight <= 1:
        raise ValueError("--minimum-weight must be between 0 and 1")


def main() -> None:
    args = parse_args()
    prompt_missing(args)
    args.input_path = args.input_path.expanduser().resolve()
    validate_args(args)
    profile = Profile(args.start_lon, args.start_lat, args.end_lon, args.end_lat)
    source_type = detect_input_type(args.input_path, args.input_type)
    selectors = parse_selectors(args.select)
    output_dir = (args.output_dir or Path.cwd() / "profile_output").expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix or f"{args.input_path.stem}_profile"
    output_stem = output_dir / prefix
    metadata: dict[str, object] = {
        "input_path": str(args.input_path),
        "input_type": source_type,
        "mode": args.mode,
        "profile": asdict(profile),
        "sample_step_km": args.sample_step_km,
    }
    if args.save_details:
        metadata["input_sha256"] = sha256_file(args.input_path)
    points: pd.DataFrame | None = None
    projected = pd.DataFrame()
    if source_type == "grid":
        grid, lon_name, lat_name, variable = open_grid(
            args.input_path, args.variable, args.lon_name, args.lat_name, selectors
        )
        metadata.update(
            {
                "grid_variable": variable,
                "longitude_coordinate": lon_name,
                "latitude_coordinate": lat_name,
                "selectors": selectors,
            }
        )
        if args.mode == "centerline":
            result = centerline_profile(
                grid, lon_name, lat_name, profile, args.sample_step_km, args.interpolation
            )
            metadata["interpolation"] = args.interpolation
        else:
            points = grid_to_points(
                subset_grid(
                    grid,
                    lon_name,
                    lat_name,
                    profile,
                    args.max_perp_km + args.window_half_width_km,
                ),
                lon_name,
                lat_name,
                profile,
                args.max_perp_km,
                args.window_half_width_km,
            )
            metadata["valid_grid_points_in_bounding_box"] = len(points)
    else:
        if args.mode == "centerline":
            raise ValueError("--mode centerline requires a regular grid; use corridor for XYZ/point tables")
        points, columns = read_table_points(
            args.input_path,
            args.header,
            args.delimiter,
            args.sheet_name,
            args.lon_col,
            args.lat_col,
            args.value_col,
        )
        metadata.update(columns)
        metadata["valid_source_points"] = len(points)

    if args.mode == "corridor":
        assert points is not None
        result, projected = corridor_profile(
            points,
            profile,
            args.sample_step_km,
            args.max_perp_km,
            args.window_half_width_km,
            args.weight_coefficient,
            args.minimum_weight,
        )
        metadata.update(
            {
                "maximum_cross_track_km": args.max_perp_km,
                "moving_window_half_width_km": args.window_half_width_km,
                "weight_formula": "max(1 - coefficient * radial_distance_km^2, minimum_weight)",
                "weight_coefficient": args.weight_coefficient,
                "minimum_weight": args.minimum_weight,
                "value_std_definition": "unweighted population standard deviation of selected spatial values",
                "spread_display": args.spread_style,
                "projected_points": len(projected),
            }
        )
    csv_path = output_stem.with_suffix(".csv")
    result.to_csv(csv_path, index=False, float_format="%.6f", encoding="utf-8-sig")
    dat_path = output_stem.with_suffix(".dat")
    metadata["litmod_dat_columns"] = write_litmod_dat(result, args.mode, dat_path)
    outputs = [csv_path, dat_path]
    # Optional files for point-level checks and reproducibility.
    if args.save_details and not projected.empty:
        projected_path = output_stem.with_name(output_stem.name + "_projected_points").with_suffix(".csv")
        projected.to_csv(projected_path, index=False, float_format="%.6f", encoding="utf-8-sig")
        outputs.append(projected_path)
    label = args.value_label or str(metadata.get("grid_variable") or metadata.get("value_column") or "Value")
    if not args.no_plot:
        outputs.extend(plot_result(result, label, args.mode, args.spread_style, output_stem))
    if args.save_details:
        metadata_path = output_stem.with_name(output_stem.name + "_provenance").with_suffix(".json")
        metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
        outputs.append(metadata_path)
    for path in outputs:
        print(f"Saved {path}")


if __name__ == "__main__":
    main()
