#!/usr/bin/env python3
"""Generate periodic table heatmaps of activation metrics.

Reads element_data.json from generate_data.py and produces color-coded
periodic tables showing waste disposal rating (Fetter limits) and photon
contact dose rate (Sv/hr) at each cooling time.

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
import matplotlib.ticker as ticker
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

    # Discrete colormap with n_levels bins, log-scaled
    cmap = plt.get_cmap(cmap_name, n_levels)
    norm = LogNorm(vmin=lo, vmax=hi)

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
        cbar.ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
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
    key = "waste_disposal_rating" if metric == "wdr" else "contact_dose_Sv_per_hr"
    units = "WDR (Fetter)" if metric == "wdr" else "Sv/hr"
    label = cooling[cooling_idx]

    vals = _extract_values(data, key, cooling_idx)
    metric_label = "Waste Disposal Rating (Fetter)" if metric == "wdr" else "Contact Dose Rate"
    irr = data["metadata"]["irradiation_fpy"]
    wl = data["metadata"]["wall_loading_MW_m2"]
    title = (
        f"{metric_label} at {label}\n"
        f"{irr} FPY DT First-Wall Irradiation ({wl} MW/m²)"
    )

    cmap_name = "viridis" if metric == "dose" else "YlOrRd"
    fig, ax = plt.subplots(figsize=(18, 9))
    draw_periodic_table(ax, vals, title, units, cmap_name=cmap_name)
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ============================================================================
# Summary plots — problematic elements by category
# ============================================================================

# Cooling time categories mapped to data array indices
# (indices match COOLING_LABELS in generate_data.py)
COOLING_INDEX = {
    "shutdown": 0, "1 hour": 1, "1 day": 2, "1 week": 3,
    "1 month": 4, "1 year": 5, "5 years": 6, "10 years": 7,
    "50 years": 8, "100 years": 9,
}

SUMMARY_CATEGORIES = [
    # (label, metric_key, cooling_time, filename_tag)
    ("Problematic for Short-Term Maintenance (1 day)",
     "contact_dose_Sv_per_hr", "1 day", "dose_short"),
    ("Problematic for Medium-Term Maintenance (1 week)",
     "contact_dose_Sv_per_hr", "1 week", "dose_medium"),
    ("Problematic for Long-Term Maintenance (1 month)",
     "contact_dose_Sv_per_hr", "1 month", "dose_long"),
    ("Problematic for Short-Term Waste (1 year)",
     "waste_disposal_rating", "1 year", "waste_short"),
    ("Problematic for Medium-Term Waste (10 years)",
     "waste_disposal_rating", "10 years", "waste_medium"),
    ("Problematic for Long-Term Waste (100 years)",
     "waste_disposal_rating", "100 years", "waste_long"),
]


def _find_top_n_at_time(data, metric_key, cooling_idx, n=10):
    """Return the set of top *n* elements at a specific cooling time."""
    vals = _extract_values(data, metric_key, cooling_idx)
    ranked = sorted(vals, key=lambda e: vals[e], reverse=True)
    return set(ranked[:n])


def _draw_red_green_table(ax, data, top_set, title, n):
    """Draw a red/green periodic table: red = in top set, green = not."""
    all_elems = set(data["elements"].keys())

    pad = 0.06
    for elem, (row, col) in PERIODIC_TABLE.items():
        x, y = col, -row

        if elem in top_set:
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

    legend_patches = [
        mpatches.Patch(facecolor="#e53935", edgecolor="#666",
                       label=f"Problematic (worst {n})"),
        mpatches.Patch(facecolor="#43a047", edgecolor="#666",
                       label="Not in problematic group"),
        mpatches.Patch(facecolor="#f5f5f5", edgecolor="#666",
                       label="No data"),
    ]
    ax.legend(handles=legend_patches, loc="lower left", fontsize=9,
              framealpha=0.9)


def plot_top_n_summaries(data, output_dir, n=10):
    """Generate 6 category summary tables plus a combined waste/dose table."""
    irr = data["metadata"]["irradiation_fpy"]
    wl = data["metadata"]["wall_loading_MW_m2"]
    subtitle = f"{irr} FPY DT First-Wall ({wl} MW/m²)"

    # Collect top sets for the combined plot
    all_top_dose = set()
    all_top_waste = set()

    for label, metric_key, cooling_time, tag in SUMMARY_CATEGORIES:
        cooling_idx = COOLING_INDEX[cooling_time]
        top_set = _find_top_n_at_time(data, metric_key, cooling_idx, n)

        if "dose" in metric_key:
            all_top_dose |= top_set
        else:
            all_top_waste |= top_set

        fig, ax = plt.subplots(figsize=(18, 9))
        _draw_red_green_table(ax, data, top_set,
                              f"{label}\n{subtitle}", n)
        path = os.path.join(output_dir, f"top{n}_{tag}.png")
        plt.savefig(path, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"Saved: {path}")

    # --- Combined: problematic waste / dose / both ---
    fig, ax = plt.subplots(figsize=(18, 9))
    all_elems = set(data["elements"].keys())
    pad = 0.06
    for elem, (row, col) in PERIODIC_TABLE.items():
        x, y = col, -row
        in_waste = elem in all_top_waste
        in_dose = elem in all_top_dose

        if in_waste and in_dose:
            color = "#b71c1c"   # dark red — both
            tc = "white"
        elif in_waste:
            color = "#e53935"   # red — waste only
            tc = "white"
        elif in_dose:
            color = "#f57c00"   # orange — dose only
            tc = "white"
        elif elem in all_elems:
            color = "#43a047"   # green
            tc = "white"
        else:
            color = "#f5f5f5"
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
    ax.set_title(
        f"Problematic Elements for Waste and/or Maintenance\n{subtitle}",
        fontsize=11, fontweight="bold", pad=10,
    )
    legend_patches = [
        mpatches.Patch(facecolor="#b71c1c", edgecolor="#666",
                       label="Problematic for waste AND dose"),
        mpatches.Patch(facecolor="#e53935", edgecolor="#666",
                       label="Problematic for waste only"),
        mpatches.Patch(facecolor="#f57c00", edgecolor="#666",
                       label="Problematic for dose only"),
        mpatches.Patch(facecolor="#43a047", edgecolor="#666",
                       label="Not in problematic group"),
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
    parser.add_argument("--metric", choices=["wdr", "dose", "both"],
                        default="both")
    parser.add_argument("--top-n", type=int, default=10,
                        help="Number of worst elements to flag as problematic in summary plots")
    args = parser.parse_args()

    with open(args.input) as f:
        data = json.load(f)

    os.makedirs(args.output_dir, exist_ok=True)
    metrics = ["wdr", "dose"] if args.metric == "both" else [args.metric]

    # All cooling tags with their indices into the data arrays
    all_cooling_tags = [
        "shutdown", "1_hour", "1_day", "1_week",
        "1_month", "1_year", "5_years", "10_years", "50_years", "100_years",
    ]

    # WDR: long-term waste-relevant timescales
    # Dose: maintenance-relevant timescales
    metric_cooling_tags = {
        "wdr": ["1_year", "10_years", "100_years"],
        "dose": ["1_day", "1_week", "1_month"],
    }

    for m in metrics:
        for tag in metric_cooling_tags[m]:
            idx = all_cooling_tags.index(tag)
            path = os.path.join(args.output_dir, f"periodic_table_{m}_{tag}.png")
            plot_single_table(data, m, idx, path)

    # Summary problematic element plots
    plot_top_n_summaries(data, args.output_dir, n=args.top_n)


if __name__ == "__main__":
    main()
