#!/usr/bin/env python
"""Plot surface-data fits and LitMod temperature or density sections.

Run interactively with ``python 001_plot_model_fit_and_section.py`` or use::

    python 001_plot_model_fit_and_section.py MODEL_DIR --style new --plot-mode both

Observation tables contain ``distance_km value [uncertainty]``. Calculated
profiles contain ``distance_km value``. Section files contain
``distance_km vertical_coordinate_km value``. Body boundaries use GMT segments
beginning with ``>``.

The standard filenames are ``0Topo.dat``, ``0HF.dat``, ``0Geoid.dat``,
``0Med_CBA.dat``, ``topo_out.dat``, ``SHF_out.dat``, ``geoid_out.dat``,
``bouguer_out.dat``, ``tempout.dat``, ``dens.dat``, and ``bodies.out``.
Filename aliases are listed in ``FILE_CANDIDATES``.
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Iterable, NamedTuple

import numpy as np


class PanelSpec(NamedTuple):
    region_y: tuple[float, float]
    frame: list[str]


class PlotConfig(NamedTuple):
    style: str
    line_pen: str
    uncoupled_pen: str
    observed_fill: str
    body_pen: str
    top_panel_yshift: str
    heatflow: PanelSpec
    geoid: PanelSpec
    bouguer: PanelSpec
    topography: PanelSpec


TEMPERATURE_CPT = """\
0 26/102/240 100 26/102/240
100 25/175/255 200 25/175/255
200 50/190/255 300 50/190/255
300 97/225/240 400 97/225/240
400 106/235/225 500 106/235/225
500 138/236/174 600 138/236/174
600 172/245/168 700 172/245/168
700 223/245/141 800 223/245/141
800 247/215/104 900 247/215/104
900 255/189/87 1000 255/189/87
1000 244/117/75 1100 244/117/75
1100 238/80/78 1200 238/80/78
1200 250/80/80 1300 250/110/110
1300 255/158/158 1400 255/180/180
1400 255/196/196 1500 255/196/196
1500 255/215/215 1600 255/215/215
B 26/102/240
F 255/215/215
N 128/128/128
"""

FILE_CANDIDATES = {
    "topography_observed": (
        "0Topo.dat",
        "ETOPO_profile_with_error_Total.dat",
        "ETOPO_profile_with_error_SN.dat",
    ),
    "bouguer_observed": (
        "0Med_CBA.dat",
        "Bouger_profile_with_error_Total.dat",
        "Bouguer_profile_with_error_Total.dat",
        "Bouger_profile_with_error_SN.dat",
        "Bouguer_Alpine_Himalaya_with_error_SN.dat",
        "csBouguer_GECO_ICGEM_with_error_SN.dat",
    ),
    "geoid_observed": (
        "0Geoid.dat",
        "GEOID_profile_with_error_Total.dat",
        "Geoid_profile_with_error.dat",
        "Geoid_profile_with_error_SN.dat",
    ),
    "heatflow_observed": (
        "0HF.dat",
        "South_ProHF.dat",
    ),
    "topography_calculated": ("topo_out.dat",),
    "bouguer_calculated": ("bouguer_out.dat",),
    "geoid_calculated": ("geoid_out.dat",),
    "heatflow_calculated": ("SHF_out.dat",),
    "temperature": ("tempout.dat",),
    "density": ("dens.dat",),
    "bodies": ("bodies.out", "bodies_GUI.out"),
}

PROFILE_PANEL_GAP_CM = 0.8


def build_config(style: str) -> PlotConfig:
    if style == "old":
        return PlotConfig(
            style="old",
            line_pen="3p,0/100/240",
            uncoupled_pen="2p,0/0/130",
            observed_fill="255/0/0",
            body_pen="0.5p,black",
            top_panel_yshift="25.6c",
            heatflow=PanelSpec(
                (40.0, 180.0),
                [
                    "x100f50+lDistance (km)",
                    "ya10f5g70+lsurface heat flow (mW/m@+2@+)",
                    "WNEs",
                ],
            ),
            geoid=PanelSpec((-8.0, 17.0), ["x100f50g100", "y5f1g10+lGeoid (m)", "WnEs"]),
            bouguer=PanelSpec(
                (-110.0, 300.0),
                ["x100f50g100+lDistance (km)", "y50f10g100+lBouguer (mGal)", "WnEs"],
            ),
            topography=PanelSpec(
                (-4000.0, 2000.0),
                ["x100f50g100", "y1000f500g2000+lElevation (m)", "WnEs"],
            ),
        )

    if style != "new":
        raise ValueError(f"Unknown style: {style}")

    return PlotConfig(
        style="new",
        line_pen="2.0p,0/82/204@35",
        uncoupled_pen="1.6p,0/82/204@55",
        observed_fill="255/158/190@35",
        body_pen="0.5p,120/120/120",
        top_panel_yshift="25c",
        heatflow=PanelSpec(
            (40.0, 150.0),
            [
                "x100f50+lDistance (km)",
                "ya30f100+lsurface heat flow (mW/m@+2@+)",
                "WNEs",
            ],
        ),
        geoid=PanelSpec((-11.0, 7.0), ["x100f50", "y10f5+lGeoid (m)", "WnEs"]),
        bouguer=PanelSpec(
            (-10.0, 290.0),
            ["x100f50+lDistance (km)", "y60f20g500+lBouguer (mGal)", "WnEs"],
        ),
        topography=PanelSpec(
            (-4200.0, 1700.0),
            ["x100f50", "y1000f500g5000+lElevation (m)", "WnEs"],
        ),
    )


def read_numeric_table(path: Path, min_cols: int = 2) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith(">"):
                continue
            values: list[float] = []
            for token in stripped.replace(",", " ").split():
                try:
                    values.append(float(token))
                except ValueError:
                    values = []
                    break
            if len(values) >= min_cols:
                rows.append(values)
    if not rows:
        raise ValueError(f"No numeric rows with at least {min_cols} columns in {path}")
    width = max(len(row) for row in rows)
    return np.array([row + [math.nan] * (width - len(row)) for row in rows], dtype=float)


def optional_numeric_table(path: Path, min_cols: int = 2) -> np.ndarray | None:
    if not path.exists():
        print(f"Warning: missing {path.name}; skipping it.")
        return None
    return read_numeric_table(path, min_cols=min_cols)


def profile_x_range(observed_topography: Path) -> tuple[float, float]:
    table = read_numeric_table(observed_topography, min_cols=1)
    return float(table[0, 0]), float(table[-1, 0])


def error_band_polygon(table: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if table.shape[1] < 3:
        raise ValueError("Error-band table needs distance, value, and error columns.")
    valid = table[:, :3]
    valid = valid[~np.isnan(valid).any(axis=1)]
    x = valid[:, 0]
    center = valid[:, 1]
    err = valid[:, 2]
    poly_x = np.concatenate([x, x[::-1], x[:1]])
    poly_y = np.concatenate([center + err, (center - err)[::-1], (center + err)[:1]])
    return poly_x, poly_y


def read_bodies_segments(path: Path) -> list[np.ndarray]:
    if not path.exists():
        return []
    segments: list[np.ndarray] = []
    current: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if current:
                    segments.append(np.array(current, dtype=float))
                    current = []
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            try:
                current.append((float(parts[0]), -float(parts[1])))
            except ValueError:
                continue
    if current:
        segments.append(np.array(current, dtype=float))
    return segments


def resolve_file(data_dir: Path, role: str, required: bool = True, quiet: bool = False) -> Path:
    """Resolve a model file from FILE_CANDIDATES."""
    candidates = FILE_CANDIDATES[role]
    existing_by_casefold = {
        path.name.casefold(): path
        for path in data_dir.iterdir()
        if path.is_file()
    }
    for index, name in enumerate(candidates):
        path = existing_by_casefold.get(name.casefold())
        if path is not None:
            if index and not quiet:
                print(f"Using {path.name} for {role.replace('_', ' ')}.")
            return path
    expected = data_dir / candidates[0]
    if required:
        tried = ", ".join(candidates)
        raise FileNotFoundError(f"Missing {role.replace('_', ' ')} file. Tried: {tried}")
    return expected


def require_files(data_dir: Path, roles: Iterable[str]) -> None:
    missing = [role for role in roles if not resolve_file(data_dir, role, required=False, quiet=True).exists()]
    if missing:
        details = [f"{role}: {', '.join(FILE_CANDIDATES[role])}" for role in missing]
        raise FileNotFoundError("Missing required inputs:\n  - " + "\n  - ".join(details))


def plot_observed_markers(fig, table: np.ndarray | None, style: str, fill: str) -> None:
    if table is None or table.shape[1] < 2:
        return
    valid = table[:, :3] if table.shape[1] >= 3 else table[:, :2]
    valid = valid[~np.isnan(valid[:, :2]).any(axis=1)]
    if valid.shape[1] >= 3:
        for x, y, err in valid[:, :3]:
            if not math.isnan(err):
                fig.plot(x=[x, x], y=[y - err, y + err], pen="1p,130/130/130")
    fig.plot(x=valid[:, 0], y=valid[:, 1], style=style, fill=fill, pen="0.25p,120/120/120")


def plot_error_band(fig, table: np.ndarray | None, fill: str) -> None:
    if table is None:
        return
    x, y = error_band_polygon(table)
    fig.plot(x=x, y=y, fill=fill, pen="0.25p,255/158/190@45")


def plot_calculated_line(fig, table: np.ndarray | None, pen: str, y_col: int = 1) -> None:
    if table is None or table.shape[1] <= y_col:
        return
    valid = table[:, [0, y_col]]
    valid = valid[~np.isnan(valid).any(axis=1)]
    if len(valid):
        fig.plot(x=valid[:, 0], y=valid[:, 1], pen=pen)


def observed_y_values(table: np.ndarray | None) -> np.ndarray:
    """Collect observations and uncertainty limits for axis scaling."""
    if table is None or table.shape[1] < 2:
        return np.array([], dtype=float)
    valid = table[:, :3] if table.shape[1] >= 3 else table[:, :2]
    valid = valid[~np.isnan(valid[:, :2]).any(axis=1)]
    if not len(valid):
        return np.array([], dtype=float)
    if valid.shape[1] >= 3:
        with_err = valid[~np.isnan(valid[:, 2])]
        if len(with_err):
            return np.concatenate([with_err[:, 1] - with_err[:, 2], with_err[:, 1] + with_err[:, 2]])
    return valid[:, 1]


def calculated_y_values(table: np.ndarray | None, y_cols: Iterable[int] = (1,)) -> np.ndarray:
    values: list[np.ndarray] = []
    if table is None:
        return np.array([], dtype=float)
    for y_col in y_cols:
        if table.shape[1] <= y_col:
            continue
        col = table[:, y_col]
        col = col[~np.isnan(col)]
        if len(col):
            values.append(col)
    if not values:
        return np.array([], dtype=float)
    return np.concatenate(values)


def nice_axis_bounds(values: Iterable[np.ndarray], step: float, default: tuple[float, float]) -> tuple[float, float]:
    merged = [np.asarray(value, dtype=float) for value in values if len(value)]
    if not merged:
        return default
    all_values = np.concatenate(merged)
    all_values = all_values[np.isfinite(all_values)]
    if not len(all_values):
        return default
    y_min = float(np.nanmin(all_values))
    y_max = float(np.nanmax(all_values))
    if math.isclose(y_min, y_max):
        y_min -= step
        y_max += step
    padding = max((y_max - y_min) * 0.08, step * 0.5)
    lower = math.floor((y_min - padding) / step) * step
    upper = math.ceil((y_max + padding) / step) * step
    return lower, upper


def make_nearneighbor_grid(pygmt, table: np.ndarray, region: list[float], spacing: float):
    valid = table[:, :3]
    valid = valid[~np.isnan(valid).any(axis=1)]
    return pygmt.nearneighbor(
        x=valid[:, 0],
        y=-valid[:, 1],
        z=valid[:, 2],
        region=region,
        spacing=spacing,
        search_radius=10,
        sectors=4,
    )


def profile_panel_yshift(profile_height_cm: float) -> str:
    """Return the vertical shift below the surface-data panels."""
    return f"{-(profile_height_cm + PROFILE_PANEL_GAP_CM):.3f}c"


def write_temperature_cpt(path: Path) -> None:
    path.write_text(TEMPERATURE_CPT, encoding="utf-8")


def colorbar_frame(plot_mode: str) -> list[str]:
    if plot_mode == "density":
        return ["x50f10+lMantle Density (kg/m3)"]
    return ["x200f100+lTemperature (deg C)"]


def move_finished_file(temp_output: Path, output: Path) -> Path:
    """Move a PyGMT result from temporary storage to the output path."""
    if output.exists():
        output.unlink()
    shutil.move(str(temp_output), str(output))
    return output


def plot_profile_panels(
    fig,
    pygmt,
    config: PlotConfig,
    data_dir: Path,
    region_x: tuple[float, float],
    width_cm: float,
    plot_mode: str,
    temp_cpt: Path,
    dens_cpt: Path,
) -> None:
    temp = read_numeric_table(resolve_file(data_dir, "temperature"), min_cols=3)
    density_path = resolve_file(data_dir, "density", required=False)
    density = optional_numeric_table(density_path, min_cols=3)
    if plot_mode == "density" and density is None:
        raise FileNotFoundError(f"Density section requested but {density_path.name} was not found.")
    region = [region_x[0], region_x[1], -5, 400]
    profile_height_cm = 405.0 * width_cm / max(region_x[1] - region_x[0], 1.0)
    projection = f"X{width_cm}c/-{profile_height_cm:.3f}c"

    temp_grid = make_nearneighbor_grid(pygmt, temp, region, spacing=5)
    density_grid = make_nearneighbor_grid(pygmt, density, region, spacing=5) if density is not None else None
    if plot_mode == "density" and density_grid is not None:
        grid = density_grid
        cmap = str(dens_cpt)
        cbar_frame = colorbar_frame("density")
    else:
        grid = temp_grid
        cmap = str(temp_cpt)
        cbar_frame = colorbar_frame("temperature")

    fig.shift_origin(yshift=profile_panel_yshift(profile_height_cm))
    fig.basemap(region=region, projection=projection, frame=["x100f50+lDistance (km)", "ya50f10+lDepth (km)", "WneS"])
    fig.grdimage(grid=grid, cmap=cmap)
    fig.colorbar(position="JMR+o0.5c/0c+w7c/0.2c", frame=cbar_frame)

    for segment in read_bodies_segments(resolve_file(data_dir, "bodies")):
        fig.plot(x=segment[:, 0], y=segment[:, 1], pen=config.body_pen)

    fig.grdcontour(grid=temp_grid, levels=[1320], pen="4p,10/10/10")

    if plot_mode == "temperature" and density_grid is not None:
        fig.grdcontour(
            grid=density_grid,
            levels=40,
            annotation="160+f6p",
            label_placement="d12c",
            limit=[3270, 3620],
            pen="1p,white",
        )
    else:
        fig.grdcontour(grid=temp_grid, levels=200, annotation="400+f6p", limit=[0, 1300], pen="1p,white")


def plot_litmod_output(
    data_dir: Path,
    output: Path,
    style: str = "new",
    plot_mode: str = "temperature",
    width_cm: float = 16.6,
) -> Path:
    import pygmt

    config = build_config(style)
    require_files(
        data_dir,
        [
            "topography_observed",
            "topography_calculated",
            "bouguer_calculated",
            "geoid_calculated",
            "heatflow_calculated",
            "temperature",
            "bodies",
        ],
    )

    topography_observed_path = resolve_file(data_dir, "topography_observed")
    region_x = profile_x_range(topography_observed_path)
    projection_top = f"X{width_cm}c/3c"
    marker_style = "c0.07c"

    pygmt.config(
        FONT_LABEL="11p",
        FONT_ANNOT_PRIMARY="11p",
        MAP_FRAME_PEN="0.5p",
        PS_MEDIA="A4",
    )

    with tempfile.TemporaryDirectory(prefix="litmod_pygmt_") as tmp:
        tmpdir = Path(tmp)
        temp_cpt = tmpdir / "temperature.cpt"
        dens_cpt = tmpdir / "density.cpt"
        write_temperature_cpt(temp_cpt)
        pygmt.makecpt(cmap="haxby", reverse=True, background=True, series=[3250, 3650, 40], output=str(dens_cpt))

        fig = pygmt.Figure()
        fig.shift_origin(xshift="2c", yshift=config.top_panel_yshift)

        surface_panels = (
            ("heatflow", config.heatflow, 10.0, False, (1,), (config.line_pen,)),
            ("geoid", config.geoid, 3.0, True, (1,), (config.line_pen,)),
            ("bouguer", config.bouguer, 50.0, True, (1,), (config.line_pen,)),
            ("topography", config.topography, 1000.0, True, (1, 2), (config.line_pen, config.uncoupled_pen)),
        )
        for index, (role, panel, step, has_error_band, y_columns, pens) in enumerate(surface_panels):
            if index:
                fig.shift_origin(yshift="-3.5c")
            observed_path = topography_observed_path if role == "topography" else resolve_file(
                data_dir, f"{role}_observed", required=False
            )
            observed = (
                read_numeric_table(observed_path, min_cols=2)
                if role == "topography"
                else optional_numeric_table(observed_path, min_cols=2)
            )
            calculated = read_numeric_table(resolve_file(data_dir, f"{role}_calculated"), min_cols=2)
            y_range = nice_axis_bounds(
                [observed_y_values(observed), calculated_y_values(calculated, y_cols=y_columns)],
                step=step,
                default=panel.region_y,
            )
            fig.basemap(
                region=[region_x[0], region_x[1], *y_range],
                projection=projection_top,
                frame=panel.frame,
            )
            if config.style == "new" and has_error_band:
                plot_error_band(fig, observed, config.observed_fill)
            else:
                plot_observed_markers(fig, observed, marker_style, config.observed_fill)
            for y_column, pen in zip(y_columns, pens):
                plot_calculated_line(fig, calculated, pen, y_col=y_column)

        plot_profile_panels(
            fig=fig,
            pygmt=pygmt,
            config=config,
            data_dir=data_dir,
            region_x=region_x,
            width_cm=width_cm,
            plot_mode=plot_mode,
            temp_cpt=temp_cpt,
            dens_cpt=dens_cpt,
        )

        # Use an ASCII temporary path for Ghostscript on Windows.
        temp_output = tmpdir / output.name
        fig.savefig(str(temp_output))
        output = move_finished_file(temp_output, output)
    return output


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot LitMod output with PyGMT.")
    parser.add_argument("model_dir", nargs="?", type=Path, help="Completed LitMod model directory.")
    parser.add_argument(
        "--style",
        choices=["new", "old"],
        help="Output style. If omitted during interactive use: 1=new, 2=legacy scatter.",
    )
    parser.add_argument(
        "--plot-mode",
        choices=["temperature", "density", "both"],
        default="both",
        help="Color-filled profile panel: temperature, density, or both as separate figures.",
    )
    parser.add_argument("--output", type=Path, default=None, help="Output figure path, e.g. Data_fit_New.pdf or .png.")
    parser.add_argument("--width", type=float, default=16.6, help="Panel width in centimeters.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    interactive = args.model_dir is None
    if args.model_dir is None:
        value = input("Completed LitMod model directory: ").strip().strip('"')
        if not value:
            print("ERROR: a model directory is required.", file=sys.stderr)
            return 2
        args.model_dir = Path(value)
    if args.style is None:
        if interactive:
            choice = input("Output style [1=new error bands, 2=legacy scatter] [1]: ").strip().lower()
            if choice in {"", "1", "new"}:
                args.style = "new"
            elif choice in {"2", "old"}:
                args.style = "old"
            else:
                print("ERROR: style must be 1 or 2.", file=sys.stderr)
                return 2
        else:
            args.style = "new"
    data_dir = args.model_dir.expanduser().resolve()
    if not data_dir.is_dir():
        print(f"ERROR: model directory not found: {data_dir}", file=sys.stderr)
        return 2
    modes = ["temperature", "density"] if args.plot_mode == "both" else [args.plot_mode]
    for mode in modes:
        if args.output is None:
            output = data_dir / "python_postprocessing" / f"001_model_fit_{mode}.pdf"
        elif len(modes) == 1:
            output = args.output
        else:
            output = args.output.with_name(f"{args.output.stem}_{mode}{args.output.suffix}")
        output = output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        try:
            result = plot_litmod_output(
                data_dir=data_dir,
                output=output,
                style=args.style,
                plot_mode=mode,
                width_cm=args.width,
            )
        except Exception as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1
        print(f"Saved {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


