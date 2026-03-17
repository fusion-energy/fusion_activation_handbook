#!/usr/bin/env python3
"""Generate periodic table heatmaps of activation metrics.

Reads element_data.json from generate_data.py and produces color-coded
periodic tables showing specific activity (Bq/kg) and photon contact dose
rate (Sv/hr) at each cooling time.

Usage:
    python plot_periodic_table.py
    python plot_periodic_table.py --input results/element_data.json
    python plot_periodic_table.py --metric dose
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm, BoundaryNorm

# ============================================================================
# Periodic table layout  (row, col)  — 0-indexed, row 0 = top
# ============================================================================

PERIODIC_TABLE = {
    "H": (0, 0), "He": (0, 17),
    "Li": (1, 0), "Be": (1, 1),
    "B": (1, 12), "C": (1, 13), "N": (1, 14), "O": (1, 15), "F": (1, 16), "Ne": (1, 17),
    "Na": (2, 0), "Mg": (2, 1),
    "Al": (2, 12), "Si": (2, 13), "P": (2, 14), "S": (2, 15), "Cl": (2, 16), "Ar": (2, 17),
    "K": (3, 0), "Ca": (3, 1), "Sc": (3, 2), "Ti": (3, 3), "V": (3, 4), "Cr": (3, 5),
    "Mn": (3, 6), "Fe": (3, 7), "Co": (3, 8), "Ni": (3, 9), "Cu": (3, 10), "Zn": (3, 11),
    "Ga": (3, 12), "Ge": (3, 13), "As": (3, 14), "Se": (3, 15), "Br": (3, 16), "Kr": (3, 17),
    "Rb": (4, 0), "Sr": (4, 1), "Y": (4, 2), "Zr": (4, 3), "Nb": (4, 4), "Mo": (4, 5),
    "Tc": (4, 6), "Ru": (4, 7), "Rh": (4, 8), "Pd": (4, 9), "Ag": (4, 10), "Cd": (4, 11),
    "In": (4, 12), "Sn": (4, 13), "Sb": (4, 14), "Te": (4, 15), "I": (4, 16), "Xe": (4, 17),
    "Cs": (5, 0), "Ba": (5, 1),
    "Hf": (5, 3), "Ta": (5, 4), "W": (5, 5), "Re": (5, 6), "Os": (5, 7), "Ir": (5, 8),
    "Pt": (5, 9), "Au": (5, 10), "Hg": (5, 11), "Tl": (5, 12), "Pb": (5, 13), "Bi": (5, 14),
    # Lanthanides — separate row with visual gap
    "La": (7, 3), "Ce": (7, 4), "Pr": (7, 5), "Nd": (7, 6), "Pm": (7, 7),
    "Sm": (7, 8), "Eu": (7, 9), "Gd": (7, 10), "Tb": (7, 11), "Dy": (7, 12),
    "Ho": (7, 13), "Er": (7, 14), "Tm": (7, 15), "Yb": (7, 16), "Lu": (7, 17),
}

ELEMENT_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
    "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
    "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
}


# ============================================================================
# Drawing helpers
# ============================================================================

def _brightness(rgba):
    return 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]


def draw_periodic_table(ax, values, title, units,
                        cmap_name="YlOrRd", vmin=None, vmax=None,
                        show_colorbar=True):
    """Draw one periodic-table heatmap on *ax*.

    Parameters
    ----------
    values : dict[str, float]
        Element symbol -> metric value.
    show_colorbar : bool
        Whether to draw a colorbar on this axes.
    """
    n_levels = 10
    pos = {k: v for k, v in values.items() if isinstance(v, (int, float)) and v > 0}

    # Determine range
    if vmin is not None and vmax is not None:
        lo, hi = max(vmin, 1e-30), vmax
    elif pos:
        lo, hi = max(min(pos.values()), 1e-30), max(pos.values())
    else:
        lo, hi = 1.0, 10.0

    # Discrete colormap with n_levels bins (log-spaced boundaries)
    bounds = np.logspace(np.log10(lo), np.log10(hi), n_levels + 1)
    base_cmap = plt.get_cmap(cmap_name, n_levels)
    norm = BoundaryNorm(bounds, ncolors=base_cmap.N)
    cmap = base_cmap

    pad = 0.06
    for elem, (row, col) in PERIODIC_TABLE.items():
        x, y = col, -row
        val = values.get(elem)
        has_data = val is not None

        if has_data and val > 0:
            color = cmap(norm(val))
            tc = "white" if _brightness(color) < 0.5 else "black"
        elif has_data:
            color = cmap(0.0)  # bottom of colormap for zero activation
            tc = "black"
        else:
            color = "#f5f5f5"  # light gray — no data
            tc = "#aaaaaa"

        rect = mpatches.FancyBboxPatch(
            (x + pad, y - 1 + pad), 1 - 2 * pad, 1 - 2 * pad,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="#666666", linewidth=0.5,
        )
        ax.add_patch(rect)

        ax.text(x + 0.5, y - 0.4, elem,
                ha="center", va="center", fontsize=8,
                fontweight="bold", color=tc)
        ax.text(x + 0.5, y - 0.15, str(ELEMENT_Z.get(elem, "")),
                ha="center", va="center", fontsize=5, color=tc, alpha=0.7)

    # Lanthanide pointer
    ax.text(2.5, -5.5, "*", ha="center", va="center", fontsize=14, color="#666")
    ax.text(2.5, -7.5, "*", ha="center", va="center", fontsize=14, color="#666")

    ax.set_xlim(-0.3, 18.3)
    ax.set_ylim(-8.3, 0.8)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)

    if show_colorbar and (vmin is not None or pos):
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.55, aspect=20, pad=0.02)
        cbar.set_label(units, fontsize=9)
        cbar.ax.tick_params(labelsize=7)


# ============================================================================
# Multi-panel and single-panel figures
# ============================================================================

def _extract_values(data, metric_key, cooling_idx):
    """Pull one cooling-time slice out of the results dict."""
    vals = {}
    for elem, info in data["elements"].items():
        series = info[metric_key]
        if cooling_idx < len(series):
            v = series[cooling_idx]
            vals[elem] = float(v) if v is not None else 0.0
    return vals


def _global_range(data, metric_key):
    """Return (vmin, vmax) across all cooling times for consistent colour scale."""
    pos = []
    for info in data["elements"].values():
        for v in info[metric_key]:
            if isinstance(v, (int, float)) and v > 0:
                pos.append(v)
    if pos:
        return min(pos), max(pos)
    return 1.0, 10.0


def plot_single_table(data, metric, cooling_idx, output_path):
    """Single periodic-table figure for one metric at one cooling time."""
    cooling = data["metadata"]["cooling_times"]
    key = "activity_Bq_per_kg" if metric == "activity" else "contact_dose_Sv_per_hr"
    units = "Bq/kg" if metric == "activity" else "Sv/hr"
    label = cooling[cooling_idx]

    # Use global range across all cooling times so colour scale is consistent
    vmin, vmax = _global_range(data, key)
    vals = _extract_values(data, key, cooling_idx)
    metric_label = "Specific Activity" if metric == "activity" else "Contact Dose Rate"
    irr = data["metadata"]["irradiation_fpy"]
    wl = data["metadata"]["wall_loading_MW_m2"]
    title = (
        f"{metric_label} at {label}\n"
        f"{irr} FPY DT First-Wall Irradiation ({wl} MW/m²)"
    )

    cmap_name = "viridis" if metric == "dose" else "YlOrRd"
    fig, ax = plt.subplots(figsize=(18, 9))
    draw_periodic_table(ax, vals, title, units, cmap_name=cmap_name,
                        vmin=vmin, vmax=vmax)
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ============================================================================
# "Ever in top N" summary plots
# ============================================================================

def _find_top_n_elements(data, metric_key, n=10):
    """Return the set of elements that appear in the top *n* at any cooling time."""
    cooling_times = data["metadata"]["cooling_times"]
    top_set = set()
    for idx in range(len(cooling_times)):
        vals = _extract_values(data, metric_key, idx)
        ranked = sorted(vals, key=lambda e: vals[e], reverse=True)
        top_set.update(ranked[:n])
    return top_set


def plot_top_n_table(ax, data, metric_key, n, title):
    """Draw a red/green periodic table: red = ever in top *n*, green = never."""
    top = _find_top_n_elements(data, metric_key, n)
    all_elems = set(data["elements"].keys())

    pad = 0.06
    for elem, (row, col) in PERIODIC_TABLE.items():
        x, y = col, -row

        if elem in top:
            color = "#e53935"   # red
            tc = "white"
        elif elem in all_elems:
            color = "#43a047"   # green
            tc = "white"
        else:
            color = "#f5f5f5"   # light gray — no data
            tc = "#aaaaaa"

        rect = mpatches.FancyBboxPatch(
            (x + pad, y - 1 + pad), 1 - 2 * pad, 1 - 2 * pad,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="#666666", linewidth=0.5,
        )
        ax.add_patch(rect)

        ax.text(x + 0.5, y - 0.4, elem,
                ha="center", va="center", fontsize=8,
                fontweight="bold", color=tc)
        ax.text(x + 0.5, y - 0.15, str(ELEMENT_Z.get(elem, "")),
                ha="center", va="center", fontsize=5, color=tc, alpha=0.7)

    ax.text(2.5, -5.5, "*", ha="center", va="center", fontsize=14, color="#666")
    ax.text(2.5, -7.5, "*", ha="center", va="center", fontsize=14, color="#666")

    ax.set_xlim(-0.3, 18.3)
    ax.set_ylim(-8.3, 0.8)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)

    # Legend
    legend_patches = [
        mpatches.Patch(facecolor="#e53935", edgecolor="#666", label=f"Top {n} at any cooling time"),
        mpatches.Patch(facecolor="#43a047", edgecolor="#666", label="Never in top group"),
        mpatches.Patch(facecolor="#f5f5f5", edgecolor="#666", label="No data"),
    ]
    ax.legend(handles=legend_patches, loc="lower left", fontsize=9,
              framealpha=0.9)


def plot_top_n_summaries(data, output_dir, n=10):
    """Generate three summary periodic tables: activity, dose, combined."""
    irr = data["metadata"]["irradiation_fpy"]
    wl = data["metadata"]["wall_loading_MW_m2"]
    subtitle = f"{irr} FPY DT First-Wall ({wl} MW/m²)"

    act_key = "activity_Bq_per_kg"
    dose_key = "contact_dose_Sv_per_hr"

    # Activity
    fig, ax = plt.subplots(figsize=(18, 9))
    plot_top_n_table(ax, data, act_key, n,
                     f"Elements Ever in Top {n} by Specific Activity\n{subtitle}")
    path = os.path.join(output_dir, f"top{n}_activity.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")

    # Dose
    fig, ax = plt.subplots(figsize=(18, 9))
    plot_top_n_table(ax, data, dose_key, n,
                     f"Elements Ever in Top {n} by Contact Dose Rate\n{subtitle}")
    path = os.path.join(output_dir, f"top{n}_dose.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")

    # Combined: union of top-N activity and top-N dose
    top_act = _find_top_n_elements(data, act_key, n)
    top_dose = _find_top_n_elements(data, dose_key, n)
    combined = top_act | top_dose

    fig, ax = plt.subplots(figsize=(18, 9))
    # Draw manually since combined comes from two sources
    all_elems = set(data["elements"].keys())
    pad = 0.06
    for elem, (row, col) in PERIODIC_TABLE.items():
        x, y = col, -row
        in_act = elem in top_act
        in_dose = elem in top_dose

        if in_act and in_dose:
            color = "#b71c1c"   # dark red — both
            tc = "white"
            marker = "A+D"
        elif in_act:
            color = "#e53935"   # red — activity only
            tc = "white"
            marker = "A"
        elif in_dose:
            color = "#f57c00"   # orange — dose only
            tc = "white"
            marker = "D"
        elif elem in all_elems:
            color = "#43a047"   # green
            tc = "white"
            marker = ""
        else:
            color = "#f5f5f5"
            tc = "#aaaaaa"
            marker = ""

        rect = mpatches.FancyBboxPatch(
            (x + pad, y - 1 + pad), 1 - 2 * pad, 1 - 2 * pad,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="#666666", linewidth=0.5,
        )
        ax.add_patch(rect)
        ax.text(x + 0.5, y - 0.4, elem,
                ha="center", va="center", fontsize=8,
                fontweight="bold", color=tc)
        ax.text(x + 0.5, y - 0.15, str(ELEMENT_Z.get(elem, "")),
                ha="center", va="center", fontsize=5, color=tc, alpha=0.7)
        if marker:
            ax.text(x + 0.5, y - 0.7, marker,
                    ha="center", va="center", fontsize=5, color=tc, alpha=0.8)

    ax.text(2.5, -5.5, "*", ha="center", va="center", fontsize=14, color="#666")
    ax.text(2.5, -7.5, "*", ha="center", va="center", fontsize=14, color="#666")
    ax.set_xlim(-0.3, 18.3)
    ax.set_ylim(-8.3, 0.8)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        f"Elements Ever in Top {n} by Activity and/or Dose\n{subtitle}",
        fontsize=11, fontweight="bold", pad=10,
    )
    legend_patches = [
        mpatches.Patch(facecolor="#b71c1c", edgecolor="#666",
                       label=f"Top {n} activity AND dose"),
        mpatches.Patch(facecolor="#e53935", edgecolor="#666",
                       label=f"Top {n} activity only"),
        mpatches.Patch(facecolor="#f57c00", edgecolor="#666",
                       label=f"Top {n} dose only"),
        mpatches.Patch(facecolor="#43a047", edgecolor="#666",
                       label="Never in top group"),
        mpatches.Patch(facecolor="#f5f5f5", edgecolor="#666",
                       label="No data"),
    ]
    ax.legend(handles=legend_patches, loc="lower left", fontsize=9,
              framealpha=0.9)

    path = os.path.join(output_dir, f"top{n}_combined.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Periodic table heatmaps of fusion activation data"
    )
    parser.add_argument("--input", default="results/element_data.json")
    parser.add_argument("--output-dir", default="results")
    parser.add_argument("--metric", choices=["activity", "dose", "both"],
                        default="both")
    parser.add_argument("--top-n", type=int, default=10,
                        help="Number of worst elements to flag in summary plots")
    args = parser.parse_args()

    with open(args.input) as f:
        data = json.load(f)

    os.makedirs(args.output_dir, exist_ok=True)
    metrics = ["activity", "dose"] if args.metric == "both" else [args.metric]

    cooling_tags = [
        "shutdown", "1_hour", "1_day", "1_week",
        "1_month", "1_year", "10_years", "100_years",
    ]
    for m in metrics:
        for idx, tag in enumerate(cooling_tags):
            path = os.path.join(args.output_dir, f"periodic_table_{m}_{tag}.png")
            plot_single_table(data, m, idx, path)

    # Summary "ever in top N" plots
    plot_top_n_summaries(data, args.output_dir, n=args.top_n)


if __name__ == "__main__":
    main()
