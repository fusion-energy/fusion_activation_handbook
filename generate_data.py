#!/usr/bin/env python3
"""Systematic activation screening of elements under DT fusion first-wall irradiation.

Depletes each element (Z=1-83) under a representative DT fusion neutron spectrum
for 40 full-power years, then tracks activity and contact dose rate through
multiple cooling times. Results saved as JSON for visualization.

Usage:
    python generate_data.py              # Run all 81 elements
    python generate_data.py --test       # Run 10 test elements
    python generate_data.py --elements Fe Cr W Ni
"""

import argparse
import json
import os
import sys
import time
import traceback

import numpy as np
import openmc
from openmc.deplete import Chain

# ============================================================================
# Element data
# ============================================================================

# Z=1-83, excluding Tc (43) and Pm (61) — no stable isotopes
ALL_ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
    'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Sm',
    'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
]

TEST_ELEMENTS = ['Fe', 'Cr', 'W', 'Ni', 'Co', 'Mo', 'Nb', 'V', 'Ti', 'Ta']

ELEMENT_Z = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
    'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Ru': 44,
    'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51,
    'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
    'Pr': 59, 'Nd': 60, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
    'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
    'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
    'Tl': 81, 'Pb': 82, 'Bi': 83,
}

# ============================================================================
# Irradiation parameters
# ============================================================================

HOUR = 3600
DAY = 86400
YEAR = 365.25 * DAY

WALL_LOADING_MW_M2 = 1.0    # MW/m² neutron wall loading
FLUX = 4.4e14                # n/cm²/s for ~1 MW/m² DT wall loading
IRRADIATION_YEARS = 40       # full-power years

# Timesteps: irradiation followed by cooling intervals
# After these steps we have snapshots at:
#   shutdown, +1h, +1d, +1wk, +1mo, +1yr, +5yr, +10yr, +50yr, +100yr
TIMESTEPS = [
    40 * YEAR,     # irradiation
    1 * HOUR,      # cool to shutdown + 1 hour
    23 * HOUR,     # cool to 1 day
    6 * DAY,       # cool to 1 week
    21 * DAY,      # cool to ~1 month (28 d total)
    337 * DAY,     # cool to 1 year
    4 * YEAR,      # cool to 5 years
    5 * YEAR,      # cool to 10 years
    40 * YEAR,     # cool to 50 years
    50 * YEAR,     # cool to 100 years
]

SOURCE_RATES = [FLUX, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

COOLING_LABELS = [
    "shutdown", "1 hour", "1 day", "1 week",
    "1 month", "1 year", "5 years", "10 years", "50 years", "100 years",
]


# ============================================================================
# Spectrum
# ============================================================================

def load_spectrum_file(path):
    """Load a FISPACT-style spectrum file (e.g. HCPB-FW 616-group).

    File format:
      - 617 energy boundaries in eV (descending, high to low)
      - blank line
      - 616 flux values (same high-to-low order)
      - optional trailing normalization and label lines

    Returns (flux, energy_boundaries) both in ascending energy order
    as OpenMC expects.
    """
    with open(path) as f:
        lines = f.read().strip().split("\n")

    # Parse all numeric values
    values = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            values.append(float(line))
        except ValueError:
            break  # stop at non-numeric lines (e.g. filename label)

    # 617 boundaries + 616 flux = 1233, possibly +1 normalization = 1234
    n_boundaries = (len(values) + 1) // 2
    n_groups = n_boundaries - 1

    energies = np.array(values[:n_boundaries])    # descending
    flux = np.array(values[n_boundaries:n_boundaries + n_groups])

    # Reverse to ascending energy order for OpenMC
    energies = energies[::-1]
    flux = flux[::-1]

    peak_idx = np.argmax(flux)
    peak_mid = 0.5 * (energies[peak_idx] + energies[peak_idx + 1]) / 1e6
    print(f"Loaded spectrum: {path}")
    print(f"  {n_groups} groups, {energies[0]:.2e} – {energies[-1]:.2e} eV")
    print(f"  Peak at {peak_mid:.2f} MeV (group {peak_idx})")

    # The HCPB-FW file uses the LLNL-616 group structure
    return flux, "LLNL-616"



# ============================================================================
# Depletion driver
# ============================================================================

def deplete_element(element, spectrum, energy_groups, chain_file=None):
    """Deplete a pure element and extract activation metrics at each cooling time.

    Returns dict with activity (Bq/kg), contact dose (Sv/hr), and dominant
    contributing nuclides at each cooling-time snapshot.
    """
    mat = openmc.Material()
    mat.add_element(element, 1.0)
    mat.set_density("g/cm3", 1.0)
    mat.volume = 1.0
    mat.depletable = True
    mat.temperature = 293.6

    results_list = mat.deplete(
        multigroup_flux=spectrum,
        energy_group_structure=energy_groups,
        timesteps=TIMESTEPS,
        source_rates=SOURCE_RATES,
        timestep_units="s",
        chain_file=chain_file,
    )

    # results_list[0] = initial (pristine), [1] = end of irradiation, …
    depleted = results_list[1:]

    activities = []
    dose_rates = []
    dominant_act = []
    dominant_dose = []

    for mat_step in depleted:
        # --- specific activity ---
        act = mat_step.get_activity(units="Bq/kg")
        act_nuc = mat_step.get_activity(units="Bq/kg", by_nuclide=True)
        activities.append(act)
        if act_nuc:
            dominant_act.append(max(act_nuc, key=act_nuc.get))
        else:
            dominant_act.append(None)

        # --- contact dose rate ---
        dose = mat_step.get_photon_contact_dose_rate(dose_quantity="effective")
        dose_nuc = mat_step.get_photon_contact_dose_rate(
            dose_quantity="effective", by_nuclide=True
        )
        dose_rates.append(dose)
        if isinstance(dose_nuc, dict) and dose_nuc:
            dominant_dose.append(max(dose_nuc, key=dose_nuc.get))
        else:
            dominant_dose.append(None)

    return {
        "Z": ELEMENT_Z[element],
        "symbol": element,
        "activity_Bq_per_kg": activities,
        "contact_dose_Sv_per_hr": dose_rates,
        "dominant_nuclide_activity": dominant_act,
        "dominant_nuclide_dose": dominant_dose,
    }


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Activation screening of elements under DT fusion irradiation"
    )
    parser.add_argument(
        "--test", action="store_true",
        help="Run only 10 test elements (Fe, Cr, W, Ni, Co, Mo, Nb, V, Ti, Ta)",
    )
    parser.add_argument("--elements", nargs="+", default=None,
                        help="Specific elements to run")
    parser.add_argument("--chain", default=None,
                        help="Path to depletion chain XML file")
    parser.add_argument("--cross_sections", default=None,
                        help="Path to cross_sections.xml file")
    parser.add_argument("--spectrum", default="hcpb_fw_616.txt",
                        help="Path to FISPACT-style spectrum file (default: hcpb_fw_616.txt)")
    parser.add_argument("--output", default="results/element_data.json",
                        help="Output JSON path")
    args = parser.parse_args()

    if args.cross_sections:
        openmc.config["cross_sections"] = args.cross_sections
    if args.chain:
        openmc.config["chain_file"] = args.chain

    if args.elements:
        elements = args.elements
    elif args.test:
        elements = TEST_ELEMENTS
    else:
        elements = ALL_ELEMENTS

    print(f"Activation screening: {len(elements)} elements")
    print(f"Irradiation: {IRRADIATION_YEARS} FPY @ {WALL_LOADING_MW_M2} MW/m²")
    print(f"Flux: {FLUX:.2e} n/cm²/s")
    print(f"Cooling times: {', '.join(COOLING_LABELS)}\n")

    spectrum, energy_groups = load_spectrum_file(args.spectrum)
    print()

    # Load the full chain once, then reduce per-element
    full_chain = Chain.from_xml(args.chain) if args.chain else Chain.from_xml(
        openmc.config["chain_file"]
    )
    print(f"Chain: {len(full_chain)} nuclides (will reduce to level 5 per element)\n")

    results = {}
    failed = []

    for i, element in enumerate(elements):
        t0 = time.time()
        print(
            f"[{i+1}/{len(elements)}] {element} (Z={ELEMENT_Z[element]})...",
            end=" ", flush=True,
        )
        try:
            # Reduce chain to nuclides reachable within 5 transmutation
            # steps from this element's natural isotopes
            mat_tmp = openmc.Material()
            mat_tmp.add_element(element, 1.0)
            initial_nuclides = list(mat_tmp.get_nuclides())
            reduced_chain = full_chain.reduce(initial_nuclides, level=5)

            data = deplete_element(element, spectrum, energy_groups,
                                   chain_file=reduced_chain)
            results[element] = data
            dt = time.time() - t0
            act0 = data["activity_Bq_per_kg"][0]
            dose0 = data["contact_dose_Sv_per_hr"][0]
            print(f"OK ({dt:.1f}s)  shutdown: {act0:.2e} Bq/kg  {dose0:.2e} Sv/hr")
        except Exception as exc:
            dt = time.time() - t0
            print(f"FAILED ({dt:.1f}s): {exc}")
            traceback.print_exc()
            failed.append(element)

    # ---- save ----
    output = {
        "metadata": {
            "description": (
                "Systematic activation screening of elements "
                "under DT fusion first-wall irradiation"
            ),
            "spectrum": "Approximate DT first-wall, VITAMIN-J-42",
            "wall_loading_MW_m2": WALL_LOADING_MW_M2,
            "total_flux_n_cm2_s": FLUX,
            "irradiation_fpy": IRRADIATION_YEARS,
            "cooling_times": COOLING_LABELS,
            "density_g_cm3": 1.0,
            "volume_cm3": 1.0,
            "notes": (
                "Density/volume are arbitrary; Bq/kg and Sv/hr are "
                "density-independent in this transport-free framework."
            ),
        },
        "elements": results,
    }

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\nSaved: {args.output}")
    if failed:
        print(f"Failed elements: {', '.join(failed)}")
    print(f"Successful: {len(results)}/{len(elements)}")


if __name__ == "__main__":
    main()
