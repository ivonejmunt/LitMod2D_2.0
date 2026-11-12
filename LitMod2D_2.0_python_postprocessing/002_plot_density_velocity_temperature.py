"""Plot LitMod density and temperature sections with property columns.

Run ``python 002_plot_density_velocity_temperature.py`` and enter the model
directory plus the reference-column distances. Press Enter to use the three
interior quarter points. Command-line use::

    python 002_plot_density_velocity_temperature.py MODEL_DIR --references 200,500,910

Required files: ``dens.dat``, ``velocities.dat``, and ``tempout.dat``.
``bodies.out`` supplies the boundary lines. Figures are written to
``MODEL_DIR/python_postprocessing`` unless ``--output-dir`` is set.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from scipy.interpolate import griddata


plt.rcParams.update(
    {
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


@dataclass
class Settings:
    model_dir: Path
    output_dir: Path
    mode: str
    references: list[float]
    max_depth: float
    horizontal_step: float
    vertical_step: float
    crust_density_limit: float
    lab_temperature: float


def read_numeric_table(path: Path, minimum_columns: int) -> np.ndarray:
    """Read numeric rows, skipping headers, comments, and count lines."""
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", ">")):
                continue
            fields = stripped.replace(",", " ").split()
            if len(fields) < minimum_columns:
                continue
            try:
                rows.append([float(value) for value in fields[:minimum_columns]])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No rows with {minimum_columns} numeric columns were found in {path}")
    return np.asarray(rows, dtype=float)


def litmod_depth(vertical_coordinate: np.ndarray) -> np.ndarray:
    """Convert the LitMod vertical coordinate to positive depth."""
    return -np.asarray(vertical_coordinate, dtype=float)


def read_body_segments(path: Path) -> list[np.ndarray]:
    """Read LitMod body boundaries as distance-depth segments."""
    if not path.is_file():
        print(f"Body boundaries not plotted; file not found: {path}")
        return []
    segments: list[np.ndarray] = []
    current: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if current:
                    segments.append(np.asarray(current, dtype=float))
                    current = []
                continue
            fields = stripped.split()
            if len(fields) < 2:
                continue
            try:
                current.append([float(fields[0]), -float(fields[1])])
            except ValueError:
                continue
    if current:
        segments.append(np.asarray(current, dtype=float))
    return segments


def grid_section(
    data: np.ndarray,
    value_column: int,
    x_min: float,
    x_max: float,
    depth_min: float,
    depth_max: float,
    horizontal_step: float,
    vertical_step: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate scattered LitMod values to a regular distance-depth grid."""
    x = data[:, 0]
    depth = litmod_depth(data[:, 1])
    values = data[:, value_column]
    valid = (
        np.isfinite(x)
        & np.isfinite(depth)
        & np.isfinite(values)
        & (x >= x_min)
        & (x <= x_max)
        & (depth >= depth_min)
        & (depth <= depth_max)
    )
    if np.count_nonzero(valid) < 3:
        raise ValueError("Too few valid LitMod samples remain inside the requested plot region.")
    x_axis = np.arange(x_min, x_max + horizontal_step * 0.5, horizontal_step)
    depth_axis = np.arange(depth_min, depth_max + vertical_step * 0.5, vertical_step)
    grid_x, grid_depth = np.meshgrid(x_axis, depth_axis)
    points = np.column_stack((x[valid], depth[valid]))
    try:
        grid = griddata(points, values[valid], (grid_x, grid_depth), method="linear")
    except Exception as error:
        print(f"Linear interpolation failed ({error}); nearest-neighbour interpolation will be used.")
        grid = griddata(points, values[valid], (grid_x, grid_depth), method="nearest")
    if np.any(~np.isfinite(grid)):
        nearest = griddata(points, values[valid], (grid_x, grid_depth), method="nearest")
        grid = np.where(np.isfinite(grid), grid, nearest)
    return x_axis, depth_axis, grid


def plot_boundaries(ax: plt.Axes, segments: list[np.ndarray], max_depth: float) -> None:
    for segment in segments:
        visible = (segment[:, 1] >= -5.0) & (segment[:, 1] <= max_depth)
        if np.count_nonzero(visible) >= 2:
            ax.plot(segment[visible, 0], segment[visible, 1], color="black", linewidth=0.65)


def nearest_column(data: np.ndarray, reference_x: float) -> tuple[float, np.ndarray]:
    """Extract the model column nearest to a profile distance."""
    unique_x = np.unique(data[:, 0][np.isfinite(data[:, 0])])
    if unique_x.size == 0:
        raise ValueError("The model table contains no finite distance values.")
    actual_x = float(unique_x[np.argmin(np.abs(unique_x - reference_x))])
    column = data[np.isclose(data[:, 0], actual_x, rtol=0.0, atol=1e-7)]
    order = np.argsort(litmod_depth(column[:, 1]))
    return actual_x, column[order]


def one_column_crustal_average(column: np.ndarray, density_limit: float) -> float:
    depth = litmod_depth(column[:, 1])
    density = column[:, 2]
    valid = np.isfinite(depth) & np.isfinite(density) & (density < density_limit)
    if np.count_nonzero(valid) < 2:
        return np.nan
    depth = depth[valid]
    density = density[valid]
    order = np.argsort(depth)
    depth = depth[order]
    density = density[order]
    thickness = np.diff(depth)
    total = depth[-1] - depth[0]
    if total <= 0:
        return np.nan
    layer_density = 0.5 * (density[:-1] + density[1:])
    return float(np.sum(layer_density * thickness) / total)


def crustal_average_curve(data: np.ndarray, density_limit: float) -> tuple[np.ndarray, np.ndarray]:
    """Calculate 5 km crustal means and a five-column moving average."""
    distances = np.unique(data[:, 0][np.isfinite(data[:, 0])])
    sample_step = 5.0
    sample_x = np.arange(float(distances.min()), float(distances.max()) + sample_step * 0.5, sample_step)
    averages = np.full(sample_x.shape, np.nan, dtype=float)
    for index, target_x in enumerate(sample_x):
        distance = float(distances[np.argmin(np.abs(distances - target_x))])
        column = data[np.isclose(data[:, 0], distance, rtol=0.0, atol=1e-7)]
        averages[index] = one_column_crustal_average(column, density_limit)
    finite = np.isfinite(averages).astype(float)
    values = np.where(np.isfinite(averages), averages, 0.0)
    kernel = np.ones(5, dtype=float)
    numerator = np.convolve(values, kernel, mode="same")
    denominator = np.convolve(finite, kernel, mode="same")
    smoothed = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=denominator > 0,
    )
    return sample_x, smoothed


def reference_colors(number: int) -> list[tuple[float, float, float, float]]:
    return [plt.get_cmap("Dark2")(value) for value in np.linspace(0.0, 0.85, number)]


def add_reference_lines(ax: plt.Axes, references: list[float], colors: list, top: float = -5.0) -> None:
    for reference, color in zip(references, colors):
        ax.axvline(reference, color=color, linewidth=1.1, linestyle="--")
        ax.text(reference, top + 2.0, f"{reference:g} km", color=color, ha="center", va="bottom", fontsize=7)


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> list[Path]:
    outputs = [output_dir / f"{stem}.png", output_dir / f"{stem}.pdf"]
    fig.savefig(outputs[0], dpi=300, bbox_inches="tight")
    fig.savefig(outputs[1], dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outputs


def plot_density_overview(
    density: np.ndarray,
    boundaries: list[np.ndarray],
    settings: Settings,
) -> list[Path]:
    x_min = float(np.nanmin(density[:, 0]))
    x_max = float(np.nanmax(density[:, 0]))
    crust_max_depth = min(80.0, settings.max_depth)
    colors = reference_colors(len(settings.references))

    crust_x, crust_depth, crust_grid = grid_section(
        density, 2, x_min, x_max, -5.0, crust_max_depth,
        settings.horizontal_step, settings.vertical_step,
    )
    full_x, full_depth, full_grid = grid_section(
        density, 2, x_min, x_max, -5.0, settings.max_depth,
        settings.horizontal_step, settings.vertical_step,
    )
    crust_grid = np.ma.masked_where(crust_grid >= settings.crust_density_limit, crust_grid)
    curve_x, curve_density = crustal_average_curve(density, settings.crust_density_limit)

    fig, axes = plt.subplots(3, 1, figsize=(11.0, 10.0), constrained_layout=True)
    image = axes[0].pcolormesh(
        crust_x, crust_depth, crust_grid, shading="auto", cmap="turbo",
        norm=Normalize(vmin=2300.0, vmax=settings.crust_density_limit), rasterized=True,
    )
    fig.colorbar(image, ax=axes[0], label="Density (kg/m3)", pad=0.015)
    plot_boundaries(axes[0], boundaries, crust_max_depth)
    add_reference_lines(axes[0], settings.references, colors)
    axes[0].set(title=f"Crustal density below {settings.crust_density_limit:g} kg/m3", ylabel="Depth (km)")
    axes[0].set_ylim(crust_max_depth, -5.0)

    axes[1].plot(curve_x, curve_density, color="0.25", linewidth=1.4)
    if np.any(np.isfinite(curve_density)):
        mean_density = float(np.nanmean(curve_density))
        axes[1].axhline(mean_density, color="0.45", linestyle="--", linewidth=1.0, label=f"Mean = {mean_density:.0f} kg/m3")
        axes[1].legend(frameon=False)
    axes[1].set(title="Thickness-weighted average crustal density", ylabel="Density (kg/m3)")

    image = axes[2].pcolormesh(full_x, full_depth, full_grid, shading="auto", cmap="turbo", rasterized=True)
    fig.colorbar(image, ax=axes[2], label="Density (kg/m3)", pad=0.015)
    plot_boundaries(axes[2], boundaries, settings.max_depth)
    add_reference_lines(axes[2], settings.references, colors)
    axes[2].set(title="Full density section", xlabel="Distance (km)", ylabel="Depth (km)")
    axes[2].set_ylim(settings.max_depth, -5.0)
    for ax in axes:
        ax.set_xlim(x_min, x_max)
    return save_figure(fig, settings.output_dir, "002_density_overview")


def plot_temperature_velocity_overview(
    density: np.ndarray,
    velocity: np.ndarray,
    temperature: np.ndarray,
    boundaries: list[np.ndarray],
    settings: Settings,
) -> list[Path]:
    x_min = float(np.nanmin(temperature[:, 0]))
    x_max = float(np.nanmax(temperature[:, 0]))
    colors = reference_colors(len(settings.references))
    x_axis, depth_axis, temperature_grid = grid_section(
        temperature, 2, x_min, x_max, -5.0, settings.max_depth,
        settings.horizontal_step, settings.vertical_step,
    )

    fig = plt.figure(figsize=(13.0, 7.2), constrained_layout=True)
    layout = fig.add_gridspec(2, 4, height_ratios=(1.55, 1.0))
    section_ax = fig.add_subplot(layout[0, :])
    image = section_ax.pcolormesh(
        x_axis, depth_axis, temperature_grid, shading="auto", cmap="turbo",
        norm=Normalize(vmin=0.0, vmax=max(1600.0, float(np.nanpercentile(temperature[:, 2], 99)))),
        rasterized=True,
    )
    fig.colorbar(image, ax=section_ax, label="Temperature (degrees C)", pad=0.012)
    finite_grid = temperature_grid[np.isfinite(temperature_grid)]
    if finite_grid.size:
        contour_max = min(1600.0, float(np.nanmax(finite_grid)))
        contour_levels = np.arange(200.0, contour_max + 1.0, 200.0)
        if contour_levels.size:
            contours = section_ax.contour(
                x_axis, depth_axis, temperature_grid, levels=contour_levels,
                colors="white", linewidths=0.55,
            )
            section_ax.clabel(contours, fontsize=6, inline=True, fmt="%.0f")
        if float(np.nanmin(finite_grid)) <= settings.lab_temperature <= float(np.nanmax(finite_grid)):
            section_ax.contour(
                x_axis, depth_axis, temperature_grid, levels=[settings.lab_temperature],
                colors="black", linewidths=1.2, linestyles="--",
            )
    plot_boundaries(section_ax, boundaries, settings.max_depth)
    add_reference_lines(section_ax, settings.references, colors)
    section_ax.set(
        title=f"Temperature section; dashed black line = {settings.lab_temperature:g} degrees C",
        xlabel="Distance (km)", ylabel="Depth (km)",
        xlim=(x_min, x_max), ylim=(settings.max_depth, -5.0),
    )

    curve_axes = [fig.add_subplot(layout[1, index]) for index in range(4)]
    specifications = [
        (curve_axes[0], density, 2, "Density (kg/m3)"),
        (curve_axes[1], velocity, 2, "Vp (km/s)"),
        (curve_axes[2], velocity, 3, "Vs (km/s)"),
        (curve_axes[3], temperature, 2, "Temperature (degrees C)"),
    ]
    for ax, table, value_column, label in specifications:
        for requested_x, color in zip(settings.references, colors):
            actual_x, column = nearest_column(table, requested_x)
            depth = litmod_depth(column[:, 1])
            values = column[:, value_column]
            valid = np.isfinite(depth) & np.isfinite(values) & (depth >= -5.0) & (depth <= settings.max_depth)
            if table is velocity:
                valid &= values > 0.0
            ax.plot(values[valid], depth[valid], color=color, linewidth=1.25, label=f"{actual_x:g} km")
            if not np.isclose(requested_x, actual_x):
                print(f"{label}: requested {requested_x:g} km; used nearest column {actual_x:g} km.")
        if label.startswith("Temperature"):
            ax.axvline(settings.lab_temperature, color="0.45", linestyle="--", linewidth=0.9)
        ax.set(title=label, xlabel=label, ylim=(settings.max_depth, -5.0))
        ax.grid(color="0.88", linewidth=0.5)
    curve_axes[0].set_ylabel("Depth (km)")
    curve_axes[-1].legend(frameon=False, loc="best")
    return save_figure(fig, settings.output_dir, "002_temperature_velocity_overview")


def parse_references(value: str) -> list[float]:
    fields = value.replace(";", ",").replace(" ", ",").split(",")
    references = [float(field) for field in fields if field.strip()]
    if not references:
        raise ValueError("At least one reference distance is required.")
    return references


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot density, temperature, Vp, and Vs from LitMod-2D output.")
    parser.add_argument("model_dir", nargs="?", type=Path, help="Completed LitMod model directory.")
    parser.add_argument("--output-dir", type=Path, help="Directory for generated PNG and PDF files.")
    parser.add_argument("--mode", choices=("density", "thermal", "both"), default="both")
    parser.add_argument("--references", help="Comma-separated distances in km; default: profile quarter points.")
    parser.add_argument("--max-depth", type=float, default=400.0)
    parser.add_argument("--horizontal-step", type=float, default=2.0)
    parser.add_argument("--vertical-step", type=float, default=2.0)
    parser.add_argument("--crust-density-limit", type=float, default=3050.0)
    parser.add_argument("--lab-temperature", type=float, default=1320.0)
    return parser.parse_args()


def configure() -> Settings:
    args = parse_args()
    interactive = args.model_dir is None
    if args.model_dir is None:
        value = input("Completed LitMod model directory: ").strip().strip('"')
        if not value:
            raise ValueError("A model directory is required.")
        args.model_dir = Path(value)
    model_dir = args.model_dir.expanduser().resolve()
    if not model_dir.is_dir():
        raise FileNotFoundError(f"Model directory does not exist: {model_dir}")
    for filename in ("dens.dat", "velocities.dat", "tempout.dat"):
        path = model_dir / filename
        if not path.is_file():
            raise FileNotFoundError(f"Required LitMod file not found: {path}")

    density = read_numeric_table(model_dir / "dens.dat", 3)
    x_min = float(np.nanmin(density[:, 0]))
    x_max = float(np.nanmax(density[:, 0]))
    default_references = ",".join(f"{value:g}" for value in np.linspace(x_min, x_max, 5)[1:-1])
    if interactive and args.references is None:
        value = input(
            f"Reference columns for density/Vp/Vs/temperature (km) [{default_references}]: "
        ).strip()
        args.references = value or default_references
    references = parse_references(args.references or default_references)
    output_dir = (args.output_dir or (model_dir / "python_postprocessing")).expanduser().resolve()

    if args.max_depth <= 0.0 or args.horizontal_step <= 0.0 or args.vertical_step <= 0.0:
        raise ValueError("Maximum depth and grid steps must be positive.")
    outside = [value for value in references if value < x_min or value > x_max]
    if outside:
        raise ValueError(f"Reference distances outside the model range {x_min:g}-{x_max:g} km: {outside}")
    output_dir.mkdir(parents=True, exist_ok=True)
    return Settings(
        model_dir=model_dir, output_dir=output_dir, mode=args.mode, references=references,
        max_depth=args.max_depth, horizontal_step=args.horizontal_step, vertical_step=args.vertical_step,
        crust_density_limit=args.crust_density_limit, lab_temperature=args.lab_temperature,
    )


def main() -> None:
    settings = configure()
    density = read_numeric_table(settings.model_dir / "dens.dat", 3)
    velocity = read_numeric_table(settings.model_dir / "velocities.dat", 4)
    temperature = read_numeric_table(settings.model_dir / "tempout.dat", 3)
    boundaries = read_body_segments(settings.model_dir / "bodies.out")
    outputs: list[Path] = []
    if settings.mode in {"density", "both"}:
        print("Creating density overview")
        outputs.extend(plot_density_overview(density, boundaries, settings))
    if settings.mode in {"thermal", "both"}:
        print("Creating temperature and velocity overview")
        outputs.extend(plot_temperature_velocity_overview(density, velocity, temperature, boundaries, settings))
    for path in outputs:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
