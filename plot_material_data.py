"""Plot waste disposal rating vs cooling time from activation screening JSON results.

Reads material_data.json (multi-element materials) and produces one plot
per material showing WDR vs cooling time with the clearance limit.

Usage:
    python plot_material_data.py
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS_FILE = Path("results/material_data.json")

# Cooling time labels → approximate days for plotting
UNIT_TO_DAYS = {
    "hour": 1 / 24,
    "hours": 1 / 24,
    "day": 1,
    "days": 1,
    "week": 7,
    "weeks": 7,
    "month": 30.4,
    "months": 30.4,
    "year": 365.25,
    "years": 365.25,
}


def label_to_days(label):
    """Convert a cooling time label like '5 years' or 'shutdown' to days."""
    if label == "shutdown":
        return 1 / (24 * 60)  # ~1 minute placeholder
    parts = label.split()
    value = float(parts[0])
    unit = parts[1]
    return value * UNIT_TO_DAYS[unit]


with open(RESULTS_FILE) as f:
    data = json.load(f)

cooling_labels = data["metadata"]["cooling_times"]
cooling_years = np.array([label_to_days(label) / 365.25 for label in cooling_labels])
waste_limits = data["metadata"]["waste_limits"]
wall_loading = data["metadata"]["wall_loading_MW_m2"]
availability = data["metadata"].get("availability", None)
irradiation_fpy = data["metadata"]["irradiation_fpy"]

subtitle = f"{wall_loading} MW/m², {irradiation_fpy} FPY"
if availability is not None:
    subtitle += f", {availability*100:.0f}% availability"

def format_mat_title(name):
    """Format material name for plot title, e.g. 'eurofer-Lindau-spec' -> 'Eurofer (Lindau spec)'."""
    name = name.replace("-", " ")
    parts = name.split()
    if parts and parts[0].lower() == "eurofer":
        source = " ".join(parts[1:])
        if source:
            return f"Eurofer ({source})"
        return "Eurofer"
    return name

materials = data["materials"]

for mat_key, mat_data in materials.items():
    mat_name = mat_data["name"]
    wdr = mat_data["waste_disposal_rating"]
    dominant = mat_data["dominant_nuclide_wdr"]

    for xscale in ["log", "linear"]:
        fig, ax = plt.subplots(figsize=(8, 8))

        ax.plot(cooling_years, wdr, "o-", linewidth=2, color="black", label="WDR")

        # Annotate dominant nuclide at each point
        for i, (yr, val, nuc) in enumerate(zip(cooling_years, wdr, dominant)):
            if nuc and val > 0:
                ax.annotate(
                    nuc, (yr, val),
                    textcoords="offset points", xytext=(5, 8),
                    fontsize=7, color="grey",
                )

        ax.axhline(y=1.0, color="red", linestyle="--", linewidth=1.5, label="Clearance limit (WDR = 1)")

        ax.set_xscale(xscale)
        ax.set_yscale("log")
        if xscale == "log":
            ax.set_xlim(left=1e-5)
        ax.set_xlabel("Cooling time [years]")
        ax.set_ylabel(f"Waste disposal rating ({waste_limits.replace('_', ' ')})")
        ax.set_title(f"Waste disposal rating {format_mat_title(mat_name)}\n{subtitle}")
        ax.legend(loc="upper right", fontsize=9)
        ax.yaxis.set_minor_locator(plt.LogLocator(base=10, subs='auto', numticks=100))
        ax.grid(True, which="major", alpha=0.3)
        ax.grid(True, which="minor", alpha=0.15, linestyle=":")

        if xscale == "log":
            ax2 = ax.twiny()
            ax2.set_xscale("log")
            ax2.set_xlim(ax.get_xlim())
            tick_labels = [
                (1 / (365.25 * 24), "1h"),
                (1 / 365.25, "1d"),
                (7 / 365.25, "1wk"),
                (1 / 12, "1mo"),
                (1, "1yr"),
                (10, "10yr"),
                (100, "100yr"),
                (1000, "1000yr"),
                (10000, "10000yr"),
            ]
            t_min, t_max = ax.get_xlim()
            ticks = [(d, l) for d, l in tick_labels if t_min <= d <= t_max * 1.2]
            ax2.set_xticks([d for d, _ in ticks])
            ax2.set_xticklabels([l for _, l in ticks])
            ax2.minorticks_off()

        fig.tight_layout()
        suffix = f"_xlinear" if xscale == "linear" else ""
        output = f"waste_disposal_rating_{mat_key}{suffix}.png"
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Saved {output}")
        plt.close(fig)

# ── Combined comparison plots (log x and linear x) ──────────────────────────

if len(materials) > 1:
  for xscale in ["log", "linear"]:
    fig, ax = plt.subplots(figsize=(10, 8))

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i, (mat_key, mat_data) in enumerate(materials.items()):
        wdr = mat_data["waste_disposal_rating"]
        label = format_mat_title(mat_data["name"])
        ax.plot(
            cooling_years, wdr, "o-", linewidth=2,
            color=color_cycle[i % len(color_cycle)],
            label=label, markersize=4,
        )

    ax.axhline(y=1.0, color="red", linestyle="--", linewidth=1.5, label="Clearance limit (WDR = 1)")

    ax.set_xscale(xscale)
    ax.set_yscale("log")
    if xscale == "log":
        ax.set_xlim(left=1e-5)
    ax.set_xlabel("Cooling time [years]")
    ax.set_ylabel(f"Waste disposal rating ({waste_limits.replace('_', ' ')})")
    ax.set_title(f"Waste disposal rating comparison\n{subtitle}")
    ax.legend(loc="best", fontsize=8)
    ax.yaxis.set_minor_locator(plt.LogLocator(base=10, subs='auto', numticks=100))
    ax.grid(True, which="major", alpha=0.3)
    ax.grid(True, which="minor", alpha=0.15, linestyle=":")

    if xscale == "log":
        ax2 = ax.twiny()
        ax2.set_xscale("log")
        ax2.set_xlim(ax.get_xlim())
        tick_labels = [
            (1 / (365.25 * 24), "1h"),
            (1 / 365.25, "1d"),
            (7 / 365.25, "1wk"),
            (1 / 12, "1mo"),
            (1, "1yr"),
            (10, "10yr"),
            (100, "100yr"),
            (1000, "1000yr"),
            (10000, "10000yr"),
        ]
        t_min, t_max = ax.get_xlim()
        ticks = [(d, l) for d, l in tick_labels if t_min <= d <= t_max * 1.2]
        ax2.set_xticks([d for d, _ in ticks])
        ax2.set_xticklabels([l for _, l in ticks])
        ax2.minorticks_off()

    fig.tight_layout()
    suffix = f"_xlinear" if xscale == "linear" else ""
    output = f"waste_disposal_rating_comparison{suffix}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)

# ── Per-material activity plots ──────────────────────────────────────────────

NAT_U_BQ_PER_KG = 2.51e7  # natural uranium specific activity

for mat_key, mat_data in materials.items():
    if "activity_Bq_per_kg" not in mat_data:
        continue
    mat_name = mat_data["name"]
    activity = mat_data["activity_Bq_per_kg"]
    dominant = mat_data.get("dominant_nuclide_activity", [None] * len(activity))

    for xscale in ["log", "linear"]:
        fig, ax = plt.subplots(figsize=(8, 8))

        ax.plot(cooling_years, activity, "o-", linewidth=2, color="black", label="Total activity")
        ax.axhline(y=NAT_U_BQ_PER_KG, color="red", linestyle="--", linewidth=1.5, label="Natural uranium")

        for i, (yr, val, nuc) in enumerate(zip(cooling_years, activity, dominant)):
            if nuc and val > 0:
                ax.annotate(
                    nuc, (yr, val),
                    textcoords="offset points", xytext=(5, 8),
                    fontsize=7, color="grey",
                )

        ax.set_xscale(xscale)
        ax.set_yscale("log")
        if xscale == "log":
            ax.set_xlim(left=1e-5)
        ax.set_xlabel("Cooling time [years]")
        ax.set_ylabel("Specific activity [Bq/kg]")
        ax.set_title(f"Specific activity {format_mat_title(mat_name)}\n{subtitle}")
        ax.legend(loc="upper right", fontsize=9)
        ax.yaxis.set_minor_locator(plt.LogLocator(base=10, subs='auto', numticks=100))
        ax.grid(True, which="major", alpha=0.3)
        ax.grid(True, which="minor", alpha=0.15, linestyle=":")

        if xscale == "log":
            ax2 = ax.twiny()
            ax2.set_xscale("log")
            ax2.set_xlim(ax.get_xlim())
            tick_labels = [
                (1 / (365.25 * 24), "1h"),
                (1 / 365.25, "1d"),
                (7 / 365.25, "1wk"),
                (1 / 12, "1mo"),
                (1, "1yr"),
                (10, "10yr"),
                (100, "100yr"),
                (1000, "1000yr"),
                (10000, "10000yr"),
            ]
            t_min, t_max = ax.get_xlim()
            ticks = [(d, l) for d, l in tick_labels if t_min <= d <= t_max * 1.2]
            ax2.set_xticks([d for d, _ in ticks])
            ax2.set_xticklabels([l for _, l in ticks])
            ax2.minorticks_off()

        fig.tight_layout()
        suffix = f"_xlinear" if xscale == "linear" else ""
        output = f"activity_{mat_key}{suffix}.png"
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Saved {output}")
        plt.close(fig)

# ── Combined activity comparison plot ────────────────────────────────────────

has_activity = any("activity_Bq_per_kg" in m for m in materials.values())
if len(materials) > 1 and has_activity:
  for xscale in ["log", "linear"]:
    fig, ax = plt.subplots(figsize=(10, 8))

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i, (mat_key, mat_data) in enumerate(materials.items()):
        if "activity_Bq_per_kg" not in mat_data:
            continue
        activity = mat_data["activity_Bq_per_kg"]
        label = format_mat_title(mat_data["name"])
        ax.plot(
            cooling_years, activity, "o-", linewidth=2,
            color=color_cycle[i % len(color_cycle)],
            label=label, markersize=4,
        )

    ax.axhline(y=2.51e7, color="red", linestyle="--", linewidth=1.5, label="Natural uranium")

    ax.set_xscale(xscale)
    ax.set_yscale("log")
    if xscale == "log":
        ax.set_xlim(left=1e-5)
    ax.set_xlabel("Cooling time [years]")
    ax.set_ylabel("Specific activity [Bq/kg]")
    ax.set_title(f"Specific activity comparison\n{subtitle}")
    ax.legend(loc="best", fontsize=8)
    ax.yaxis.set_minor_locator(plt.LogLocator(base=10, subs='auto', numticks=100))
    ax.grid(True, which="major", alpha=0.3)
    ax.grid(True, which="minor", alpha=0.15, linestyle=":")

    if xscale == "log":
        ax2 = ax.twiny()
        ax2.set_xscale("log")
        ax2.set_xlim(ax.get_xlim())
        tick_labels = [
            (1 / (365.25 * 24), "1h"),
            (1 / 365.25, "1d"),
            (7 / 365.25, "1wk"),
            (1 / 12, "1mo"),
            (1, "1yr"),
            (10, "10yr"),
            (100, "100yr"),
            (1000, "1000yr"),
            (10000, "10000yr"),
        ]
        t_min, t_max = ax.get_xlim()
        ticks = [(d, l) for d, l in tick_labels if t_min <= d <= t_max * 1.2]
        ax2.set_xticks([d for d, _ in ticks])
        ax2.set_xticklabels([l for _, l in ticks])
        ax2.minorticks_off()

    fig.tight_layout()
    suffix = f"_xlinear" if xscale == "linear" else ""
    output = f"activity_comparison{suffix}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)
