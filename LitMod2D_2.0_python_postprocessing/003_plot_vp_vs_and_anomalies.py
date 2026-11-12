"""Plot LitMod Vp, Vs, and relative anomalies in five PyGMT panels.

Run ``python 003_plot_vp_vs_and_anomalies.py`` and choose either the bundled
ak135 reference or one model column for both Vp and Vs. Command-line use::

    python 003_plot_vp_vs_and_anomalies.py MODEL_DIR --reference-x 910

Anomalies are ``(model - reference) / reference * 100``. The velocity input is
``velocities.dat`` or ``post_processing_*.dat``. The output directory defaults
to ``MODEL_DIR/python_postprocessing``.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd

try:
    import pygmt
except ModuleNotFoundError:
    pygmt = None


def read_velocity_file(path):
    """Read a four-column velocities file or six-column post-processing file."""
    raw = pd.read_csv(path, sep=r"\s+", header=None, comment="#")
    if raw.shape[1] >= 6:
        # Six-column layout: X, Z, column 3, column 4, Vp, Vs.
        out = pd.DataFrame(
            {
                "X": raw.iloc[:, 0].astype(float),
                "Depth": np.abs(raw.iloc[:, 1].astype(float)),
                "Vp": raw.iloc[:, 4].astype(float),
                "Vs": raw.iloc[:, 5].astype(float),
            }
        )
    elif raw.shape[1] >= 4:
        # Four-column layout: X, Z, Vp, Vs.
        z = raw.iloc[:, 1].astype(float)
        out = pd.DataFrame(
            {
                "X": raw.iloc[:, 0].astype(float),
                "Depth": np.abs(z),
                "Vp": raw.iloc[:, 2].astype(float),
                "Vs": raw.iloc[:, 3].astype(float),
            }
        )
    else:
        raise ValueError(f"Velocity file needs at least 4 columns: {path}")

    out["Vp"] = out["Vp"].replace(0, np.nan)
    out["Vs"] = out["Vs"].replace(0, np.nan)
    return out.replace([np.inf, -np.inf], np.nan).dropna(subset=["X", "Depth"])


def read_reference_model(ref_file):
    """Read a reference velocity model as depth, Vp, Vs."""
    header = None
    numeric_rows = []
    with Path(ref_file).open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            try:
                values = [float(part) for part in parts]
            except ValueError:
                if header is None:
                    header = [part.lower() for part in parts]
                continue
            if len(values) >= 3:
                numeric_rows.append(values)

    if not numeric_rows:
        raise ValueError(f"Reference model has no numeric rows: {ref_file}")

    if header and "vp" in header and "vs" in header:
        z_idx = header.index("depth") if "depth" in header else 0
        vp_idx = header.index("vp")
        vs_idx = header.index("vs")
    else:
        max_cols = max(len(row) for row in numeric_rows)
        if max_cols >= 4:
            z_idx, vp_idx, vs_idx = 0, 2, 3
        else:
            z_idx, vp_idx, vs_idx = 0, 1, 2

    records = []
    for row in numeric_rows:
        if len(row) > max(z_idx, vp_idx, vs_idx):
            records.append((abs(row[z_idx]), row[vp_idx], row[vs_idx]))
    return pd.DataFrame(records, columns=["Depth", "Ref_Vp", "Ref_Vs"]).sort_values("Depth")


def read_topography(data_dir):
    for name in ("0Topo.dat", "topo.inp", "PTOT_topo.dat"):
        path = Path(data_dir) / name
        if path.exists():
            return pd.read_csv(
                path,
                sep=r"\s+",
                header=None,
                usecols=[0, 1],
                names=["X", "Elev"],
            ).dropna()
    raise FileNotFoundError(f"No topography file found in {data_dir}")


def read_bodies_segments(data_dir):
    for name in ("bodies_GUI.out", "bodies.out"):
        path = Path(data_dir) / name
        if path.exists():
            break
    else:
        return []

    segments = []
    current = []
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
                current.append((float(parts[0]), abs(float(parts[1]))))
            except ValueError:
                continue
    if current:
        segments.append(np.array(current, dtype=float))
    return segments


def move_finished_pdf(temp_pdf, output_file):
    output_file = Path(output_file)
    if output_file.exists():
        output_file.unlink()
    shutil.move(str(temp_pdf), str(output_file))
    return output_file


def make_velocity_grid(df, column, region):
    clean = df.dropna(subset=[column])
    return pygmt.nearneighbor(
        x=clean["X"],
        y=clean["Depth"],
        z=clean[column],
        region=region,
        spacing="1/0.5",
        search_radius="10",
        sectors="4",
    )


def make_reference_grid(df_ref, column, left, right, region):
    x_vals = np.arange(left, right + 10, 10)
    frames = []
    for xpos in x_vals:
        tmp = df_ref[["Depth", column]].copy()
        tmp["X"] = xpos
        frames.append(tmp[["X", "Depth", column]])
    ref_xyz = pd.concat(frames, ignore_index=True).dropna()
    ref_grid = pygmt.nearneighbor(
        x=ref_xyz["X"],
        y=ref_xyz["Depth"],
        z=ref_xyz[column],
        region=region,
        spacing="10/5",
        search_radius="10",
        sectors="4",
    )
    return pygmt.grdsample(grid=ref_grid, spacing="10/5", region=region, verbose="e")


def aligned_upper_bound(start, stop, spacing):
    """Extend a grid boundary to an integer number of intervals."""
    intervals = math.ceil((stop - start) / spacing - 1e-10)
    return start + intervals * spacing


def extract_reference_column(df, xref):
    x_unique = np.sort(df["X"].dropna().unique())
    x_near = float(x_unique[np.argmin(np.abs(x_unique - xref))])
    if abs(x_near - xref) > 1e-6:
        print(f"Warning: X={xref:g} km not found exactly; using nearest column X={x_near:g} km.")
    ref = df[np.isclose(df["X"], x_near)][["Depth", "Vp", "Vs"]].copy()
    ref.rename(columns={"Vp": "Ref_Vp", "Vs": "Ref_Vs"}, inplace=True)
    return ref.dropna().sort_values("Depth"), x_near


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model_dir", nargs="?", type=Path, help="completed LitMod model directory")
    parser.add_argument("--output-dir", type=Path, help="default: MODEL_DIR/python_postprocessing")
    parser.add_argument("--output", type=Path, help="output PDF path")
    parser.add_argument("--velocity-file", type=Path, help="default: velocities.dat or first post_processing_*.dat")
    reference = parser.add_mutually_exclusive_group()
    reference.add_argument("--reference-x", type=float, help="use the nearest model column at this distance (km)")
    reference.add_argument("--reference-model", type=Path, help="depth/Vp/Vs reference table; default: bundled ak135")
    parser.add_argument("--profile-length", type=float, help="plot length in km; default: data maximum")
    parser.add_argument("--max-depth", type=float, default=400.0, help="maximum plot depth in km")
    parser.add_argument("--anomaly-limit", type=float, default=1.0, help="symmetric anomaly color limit in percent")
    args = parser.parse_args()

    if pygmt is None:
        raise SystemExit(
            "PyGMT is not installed in this Python environment. Install GMT and PyGMT, "
            "then run this script again."
        )
    interactive = args.model_dir is None
    if args.model_dir is None:
        value = input("Completed LitMod model directory: ").strip().strip('"')
        if not value:
            parser.error("a model directory is required")
        args.model_dir = Path(value)
    if interactive and args.reference_x is None and args.reference_model is None:
        choice = input("Vp/Vs anomaly reference [1=ak135, 2=model column] [1]: ").strip().lower()
        if choice in {"2", "column", "model"}:
            value = input("Reference model-column distance (km): ").strip()
            try:
                args.reference_x = float(value)
            except ValueError:
                parser.error("reference model-column distance must be numeric")
        elif choice not in {"", "1", "ak135"}:
            parser.error("reference choice must be 1 or 2")
    data_dir = args.model_dir.expanduser().resolve()
    if not data_dir.is_dir():
        parser.error(f"model directory not found: {data_dir}")
    script_dir = Path(__file__).resolve().parent
    output_dir = (args.output_dir or data_dir / "python_postprocessing").expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    velocity_file = args.velocity_file
    if velocity_file is not None:
        velocity_file = velocity_file.expanduser().resolve()
    else:
        velocity_file = data_dir / "velocities.dat"
    if not velocity_file.exists():
        candidates = sorted(data_dir.glob("post_processing_*.dat"))
        if not candidates:
            parser.error(f"no velocities.dat or post_processing_*.dat found in {data_dir}")
        velocity_file = candidates[0]
        print(f"Using {velocity_file.name} because velocities.dat was not found.")

    df = read_velocity_file(velocity_file)
    df_topo = read_topography(data_dir)
    bodies_segments = read_bodies_segments(data_dir)

    left = 0.0
    requested_right = args.profile_length or float(df["X"].max())
    right = aligned_upper_bound(left, requested_right, 10.0)
    max_depth = aligned_upper_bound(0.0, args.max_depth, 5.0)
    region = [left, right, 0, max_depth]
    anomaly_limit = args.anomaly_limit

    if args.reference_x is not None:
        df_ref, x_near = extract_reference_column(df, args.reference_x)
        ref_label = f"Column{int(round(x_near))}"
    else:
        ref_file = (args.reference_model or script_dir / "support" / "ak135.txt").expanduser().resolve()
        if not ref_file.exists():
            parser.error(f"reference model not found: {ref_file}")
        df_ref = read_reference_model(ref_file)
        ref_label = ref_file.stem

    print(f"Creating velocity grids from {velocity_file.name}")
    grid_vp_fine = make_velocity_grid(df, "Vp", region)
    grid_vs_fine = make_velocity_grid(df, "Vs", region)

    grid_vp_sample = pygmt.grdsample(
        grid=grid_vp_fine, spacing="10/5", region=region, verbose="e"
    )
    grid_vs_sample = pygmt.grdsample(
        grid=grid_vs_fine, spacing="10/5", region=region, verbose="e"
    )
    grid_ref_vp = make_reference_grid(df_ref, "Ref_Vp", left, right, region)
    grid_ref_vs = make_reference_grid(df_ref, "Ref_Vs", left, right, region)

    grid_vp_anom = (grid_vp_sample - grid_ref_vp) / grid_ref_vp * 100
    grid_vs_anom = (grid_vs_sample - grid_ref_vs) / grid_ref_vs * 100

    output_file = args.output or output_dir / f"003_vp_vs_anomalies_{ref_label}.pdf"
    output_file = output_file.expanduser().resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    print("Creating figure")
    fig = pygmt.Figure()
    anomaly_cpt_series = f"-{anomaly_limit}/{anomaly_limit}/0.05"

    with fig.subplot(
        nrows=5,
        ncols=1,
        figsize=("18c", "48c"),
        margins=["0.2c", "1.7c"],
        title="LitMod 2D Profile",
        frame=["WSne"],
        sharex=True,
    ):
        with fig.set_panel(panel=0):
            elev_min, elev_max = df_topo["Elev"].min(), df_topo["Elev"].max()
            y_min = min(-1000, np.floor(elev_min / 500) * 500)
            y_max = max(1500, np.ceil(elev_max / 500) * 500)
            fig.plot(
                x=df_topo["X"],
                y=df_topo["Elev"],
                region=[left, right, y_min, y_max],
                projection="X?",
                pen="1p,red",
            )
            fig.plot(x=[left, right], y=[0, 0], pen="0.5p,black")
            fig.basemap(frame=["WsNE", "xaf+l ", "yaf+lElevation (m)"])

        with fig.set_panel(panel=1):
            pygmt.makecpt(cmap="haxby", series="7.6/9.0/0.01", continuous=True, reverse=True, background=True)
            fig.grdimage(grid=grid_vp_fine, region=region, projection="X?/-?", cmap=True)
            fig.colorbar(position="JMR+o0.5c/0c+w6c/0.5c", frame=["x0.2f0.1+lVelocity Vp (km/s)"])
            for segment in bodies_segments:
                fig.plot(x=segment[:, 0], y=segment[:, 1], pen="0.5p,black")
            fig.grdcontour(grid=grid_vp_fine, levels=0.1, annotation="0.2+f6p,black+gwhite", limit="7.6/9.0", pen="0.25p,white", verbose="e")
            fig.basemap(frame=["Wesn", "xaf+l ", "ya50f10+lDepth (km)"])

        with fig.set_panel(panel=2):
            pygmt.makecpt(cmap="polar", series=anomaly_cpt_series, reverse=True, continuous=True, background=True)
            fig.grdimage(grid=grid_vp_anom, region=region, projection="X?/-?", cmap=True)
            fig.colorbar(position="JMR+o0.5c/0c+w6c/0.5c", frame=["x0.5f0.1+lVp anomaly %"])
            for segment in bodies_segments:
                fig.plot(x=segment[:, 0], y=segment[:, 1], pen="0.5p,black")
            fig.grdcontour(grid=grid_vp_anom, levels=0.5, annotation="0.5+f6p,black+gwhite", limit=f"-{anomaly_limit}/-0.2", pen="0.35p,grey30", verbose="e")
            fig.grdcontour(grid=grid_vp_anom, levels=0.5, annotation="0.5+f6p,black+gwhite", limit=f"0.2/{anomaly_limit}", pen="0.35p,grey30", verbose="e")
            fig.basemap(frame=["Wesn", "xaf+l ", "ya50f10+lDepth (km)"])

        with fig.set_panel(panel=3):
            pygmt.makecpt(cmap="haxby", series="4.25/4.75/0.01", continuous=True, reverse=True, background=True)
            fig.grdimage(grid=grid_vs_fine, region=region, projection="X?/-?", cmap=True)
            fig.colorbar(position="JMR+o0.5c/0c+w6c/0.5c", frame=["x0.05f0.025+lVelocity Vs (km/s)"])
            for segment in bodies_segments:
                fig.plot(x=segment[:, 0], y=segment[:, 1], pen="0.5p,black")
            fig.grdcontour(grid=grid_vs_fine, levels=0.05, annotation="0.1+f6p,black+gwhite", limit="4.25/4.75", pen="0.25p,white", verbose="e")
            fig.basemap(frame=["Wesn", "xaf+l ", "ya50f10+lDepth (km)"])

        with fig.set_panel(panel=4):
            pygmt.makecpt(cmap="polar", series=anomaly_cpt_series, reverse=True, continuous=True, background=True)
            fig.grdimage(grid=grid_vs_anom, region=region, projection="X?/-?", cmap=True)
            fig.colorbar(position="JMR+o0.5c/0c+w6c/0.5c", frame=["x0.5f0.1+lVs anomaly %"])
            for segment in bodies_segments:
                fig.plot(x=segment[:, 0], y=segment[:, 1], pen="0.5p,black")
            fig.grdcontour(grid=grid_vs_anom, levels=0.5, annotation="0.5+f6p,black+gwhite", limit=f"-{anomaly_limit}/-0.2", pen="0.35p,grey30", verbose="e")
            fig.grdcontour(grid=grid_vs_anom, levels=0.5, annotation="0.5+f6p,black+gwhite", limit=f"0.2/{anomaly_limit}", pen="0.35p,grey30", verbose="e")
            fig.basemap(frame=["WeSn", "xaf+lDistance (km)", "ya50f10+lDepth (km)"])

    temp_dir = Path(tempfile.mkdtemp(prefix="litmod_ajay_pygmt_"))
    temp_pdf = temp_dir / output_file.name
    fig.savefig(str(temp_pdf))
    output_file = move_finished_pdf(temp_pdf, output_file)
    shutil.rmtree(temp_dir, ignore_errors=True)
    print(f"\nPlot saved to {output_file}")


if __name__ == "__main__":
    main()
