#!/usr/bin/env python3
"""Systematic activation screening of structural materials under DT fusion first-wall irradiation.

Depletes each material under a representative DT fusion neutron spectrum
for 40 full-power years, then tracks waste disposal rating (Fetter limits) and
contact dose rate through multiple cooling times. Results saved as JSON for
visualization.

Material compositions sourced from:
  https://github.com/fusion-energy/neutronics_material_maker/blob/048ccef/src/neutronics_material_maker/data/structural_materials.json

Usage:
    python generate_materials_data.py                        # Run all materials
    python generate_materials_data.py --test                 # Run 3 test materials
    python generate_materials_data.py --materials eurofer SS_316L_N_IG tungsten
"""

import argparse
import json
import os
import time
import traceback

import numpy as np
import openmc
from openmc.deplete import Chain

# ============================================================================
# Material definitions
# ============================================================================
# Each entry: elements (name -> fraction), density (g/cm3), percent_type (ao/wo)

MATERIALS = {
    # Original composition from neutronics_material_maker
    "eurofer-neutronics-handbook": {
        "elements": {
            "Fe": 0.88821, "B": 1e-05, "C": 0.00105, "N": 0.0004, "O": 1e-05,
            "Al": 4e-05, "Si": 0.00026, "P": 2e-05, "S": 3e-05, "Ti": 1e-05,
            "V": 0.002, "Cr": 0.09, "Mn": 0.0055, "Co": 5e-05, "Ni": 0.0001,
            "Cu": 3e-05, "Nb": 5e-05, "Mo": 3e-05, "Ta": 0.0012, "W": 0.011,
        },
        "density": 7.78,
        "percent_type": "ao",
    },
    # Lindau et al. 2005, Fusion Eng. Des. 75-79, pp.989-994
    # Specification nominal values (wt%)
    "eurofer-Lindau-spec": {
        "elements": {
            "C": 0.11, "Cr": 9.0, "W": 1.1, "Mn": 0.40, "V": 0.20,
            "Ta": 0.12, "N": 0.030, "P": 0.005, "S": 0.005, "B": 0.001,
            "O": 0.01, "Si": 0.05, "Ti": 0.01,
            "Fe": 88.819,  # balance
        },
        "density": 7.78,
        "percent_type": "wo",
    },
    # Lindau et al. 2005, Heat 9741 actual analysis (mid-range values, wt%)
    "eurofer-Lindau-heat9741": {
        "elements": {
            "C": 0.115, "Cr": 8.89, "W": 1.11, "Mn": 0.44, "V": 0.19,
            "Ta": 0.14, "N": 0.026, "P": 0.005, "S": 0.004, "B": 0.0007,
            "O": 0.0015, "Si": 0.05, "Ti": 0.01,
            "Nb": 0.0005, "Mo": 0.002, "Ni": 0.018, "Cu": 0.012, "Co": 0.005,
            "Al": 0.008,
            "Fe": 88.964,  # balance
        },
        "density": 7.78,
        "percent_type": "wo",
    },
    # PMC 2025, Heat 993394 batch 3 (Saarschmiede / EUROfusion, wt%)
    "eurofer-PMC-heat993394": {
        "elements": {
            "C": 0.105, "Si": 0.024, "Mn": 0.56, "P": 0.0025, "S": 0.001,
            "Ni": 0.013, "Cr": 9.08, "Mo": 0.005, "V": 0.235, "W": 1.07,
            "Ti": 0.001, "Cu": 0.005, "Nb": 0.005, "Al": 0.009, "N": 0.039,
            "Ta": 0.125, "Co": 0.007,
            "Fe": 88.718,  # balance
        },
        "density": 7.78,
        "percent_type": "wo",
    },
    # Eurofer 97 at spec limits — no Nb/Mo/Ni/Co/Cu (pure spec, no impurities)
    "eurofer-spec-clean": {
        "elements": {
            "C": 0.11, "Cr": 9.0, "W": 1.1, "Mn": 0.40, "V": 0.20,
            "Ta": 0.12, "N": 0.030, "B": 0.001,
            "Fe": 89.039,  # balance
        },
        "density": 7.78,
        "percent_type": "wo",
    },
    "SS_316L_N_IG": {
        "elements": {
            "Fe": 62.973, "C": 0.030, "Mn": 2.00, "Si": 0.50, "P": 0.03,
            "S": 0.015, "Cr": 18.00, "Ni": 12.50, "Mo": 2.70, "N": 0.080,
            "B": 0.002, "Cu": 1.0, "Co": 0.05, "Nb": 0.01, "Ti": 0.10,
            "Ta": 0.01,
        },
        "density": 7.93,
        "percent_type": "ao",
    },
    "tungsten": {
        "elements": {"W": 1.0},
        "density": 19.3,
        "percent_type": "ao",
    },
    "tungsten_with_impurities": {
        "elements": {
            "W": 0.999595, "Ag": 1e-05, "Al": 1.5e-05, "As": 5e-06,
            "Ba": 5e-06, "Ca": 5e-06, "Cd": 5e-06, "Co": 1e-05, "Cr": 2e-05,
            "Cu": 1e-05, "Fe": 3e-05, "K": 1e-05, "Mg": 5e-06, "Mn": 5e-06,
            "Na": 1e-05, "Nb": 1e-05, "Ni": 5e-06, "Pb": 5e-06, "Ta": 2e-05,
            "Ti": 5e-06, "Zn": 5e-06, "Zr": 5e-06, "Mo": 1e-04, "C": 3e-05,
            "H": 5e-06, "N": 5e-06, "O": 2e-05, "P": 2e-05, "S": 5e-06,
            "Si": 2e-05,
        },
        "density": 19.3,
        "percent_type": "ao",
    },
    "P91": {
        "elements": {
            "Fe": 89, "Cr": 9.1, "Mo": 1, "Mn": 0.5, "Si": 0.4,
        },
        "density": 7.96,
        "percent_type": "ao",
    },
    "SST91": {
        "elements": {
            "C": 0.10, "Mn": 0.45, "P": 0.02, "S": 0.01, "Si": 0.35,
            "Cr": 8.75, "Mo": 0.95, "V": 0.215, "N": 0.05, "Ni": 0.4,
            "Al": 0.04, "Nb": 0.08, "Fe": 88.585,
        },
        "density": 7.77,
        "percent_type": "wo",
    },
    "a36_steel": {
        "elements": {
            "C": 0.0025, "Si": 0.00283, "Mn": 0.0103, "Fe": 0.989,
            "Cu": 0.002, "P": 0.0004, "S": 0.0005,
        },
        "density": 7.85,
        "percent_type": "wo",
    },
    "CuCrZr": {
        "elements": {"Cu": 0.91864, "Cr": 0.08, "Zr": 0.00136},
        "density": 8.9,
        "percent_type": "wo",
    },
    "CuCrZr_with_impurities": {
        "elements": {
            "Cu": 0.9871, "Cr": 0.0075, "Zr": 0.0011, "Co": 0.0005,
            "Ta": 0.0001, "Nb": 0.001, "B": 1e-05, "O": 0.00032,
            "Mg": 0.0004, "Al": 3e-05, "Si": 0.0004, "P": 0.00014,
            "S": 4e-05, "Mn": 2e-05, "Fe": 0.0002, "Ni": 0.0006,
            "Zn": 0.0001, "As": 0.0001, "Sn": 0.0001, "Sb": 0.00011,
            "Pb": 0.0001, "Bi": 3e-05,
        },
        "density": 8.9,
        "percent_type": "ao",
    },
    "Ti-6Al-4V": {
        "elements": {
            "Ti": 89.255, "Al": 6.125, "V": 4, "Fe": 0.3, "O": 0.2,
            "C": 0.08, "N": 0.05, "H": 0.015, "Y": 0.005,
        },
        "density": 4.429,
        "percent_type": "ao",
    },
    "SiC": {
        "elements": {"Si": 1.0, "C": 1.0},
        "density": 3.21,
        "percent_type": "ao",
    },
    "zircaloy_4": {
        "elements": {
            "Zr": 0.9823, "Sn": 0.0145, "Fe": 0.0021, "Cr": 0.001,
            "O": 0.0001,
        },
        "density": 6.56,
        "percent_type": "ao",
    },
    "V-4Cr-4Ti": {
        "elements": {
            "V": 92.0, "Cr": 4.0, "Ti": 4.0,
        },
        "density": 6.05,
        "percent_type": "wo",
    },
}

ALL_MATERIALS = list(MATERIALS.keys())
TEST_MATERIALS = ["eurofer-spec-clean", "eurofer-PMC-heat993394", "SS_316L_N_IG"]

# ============================================================================
# Irradiation parameters
# ============================================================================

HOUR = 3600
DAY = 86400
YEAR = 365.25 * DAY

WALL_LOADING_MW_M2 = 1.0
AVAILABILITY = 0.8
FLUX = 4.4e14 * AVAILABILITY  # n/cm²/s for ~1 MW/m² DT wall loading    
IRRADIATION_YEARS = 5

TIMESTEPS = [
    5 * YEAR,     # irradiation
    1 * HOUR,      # cool to shutdown + 1 hour
    23 * HOUR,     # cool to 1 day
    6 * DAY,       # cool to 1 week
    21 * DAY,      # cool to ~1 month (28 d total)
    337 * DAY,     # cool to 1 year
    4 * YEAR,      # cool to 5 years
    5 * YEAR,      # cool to 10 years
    40 * YEAR,     # cool to 50 years
    50 * YEAR,     # cool to 100 years
    100 * YEAR,    # cool to 200 years
    100 * YEAR,    # cool to 300 years
    200 * YEAR,    # cool to 500 years
    500 * YEAR,    # cool to 1000 years
    1000 * YEAR,   # cool to 2000 years
    3000 * YEAR,   # cool to 5000 years
    5000 * YEAR,   # cool to 10000 years
]

SOURCE_RATES = [FLUX] + [0.0] * 16

COOLING_LABELS = [
    "shutdown", "1 hour", "1 day", "1 week",
    "1 month", "1 year", "5 years", "10 years", "50 years", "100 years",
    "200 years", "300 years", "500 years", "1000 years", "2000 years",
    "5000 years", "10000 years",
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

    Returns (flux, energy_group_structure_name) both in ascending energy order
    as OpenMC expects.
    """
    with open(path) as f:
        lines = f.read().strip().split("\n")

    values = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            values.append(float(line))
        except ValueError:
            break

    n_boundaries = (len(values) + 1) // 2
    n_groups = n_boundaries - 1

    energies = np.array(values[:n_boundaries])
    flux = np.array(values[n_boundaries:n_boundaries + n_groups])

    energies = energies[::-1]
    flux = flux[::-1]

    peak_idx = np.argmax(flux)
    peak_mid = 0.5 * (energies[peak_idx] + energies[peak_idx + 1]) / 1e6
    print(f"Loaded spectrum: {path}")
    print(f"  {n_groups} groups, {energies[0]:.2e} – {energies[-1]:.2e} eV")
    print(f"  Peak at {peak_mid:.2f} MeV (group {peak_idx})")

    return flux, "LLNL-616"


# ============================================================================
# Depletion driver
# ============================================================================

def deplete_material(name, mat_def, spectrum, energy_groups, chain_file=None):
    """Deplete a material and extract activation metrics at each cooling time.

    Returns dict with waste disposal rating (Fetter limits), contact dose
    (Sv/hr), and dominant contributing nuclides at each cooling-time snapshot.
    """
    mat = openmc.Material()
    percent_type = mat_def["percent_type"]
    for element, fraction in mat_def["elements"].items():
        mat.add_element(element, fraction, percent_type=percent_type)
    mat.set_density("g/cm3", mat_def["density"])
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

    depleted = results_list[1:]

    wdr_values = []
    dose_rates = []
    activity_values = []
    dominant_wdr = []
    dominant_dose = []
    dominant_activity = []

    for mat_step in depleted:
        wdr = mat_step.waste_disposal_rating(limits="StrlSchV_metal_recycling")
        wdr_nuc = mat_step.waste_disposal_rating(limits="StrlSchV_metal_recycling", by_nuclide=True)
        wdr_values.append(wdr)
        if wdr_nuc:
            dominant_wdr.append(max(wdr_nuc, key=wdr_nuc.get))
        else:
            dominant_wdr.append(None)

        dose = mat_step.get_photon_contact_dose_rate(dose_quantity="effective")
        dose_nuc = mat_step.get_photon_contact_dose_rate(
            dose_quantity="effective", by_nuclide=True
        )
        dose_rates.append(dose)
        if isinstance(dose_nuc, dict) and dose_nuc:
            dominant_dose.append(max(dose_nuc, key=dose_nuc.get))
        else:
            dominant_dose.append(None)

        activity = mat_step.get_activity(units="Bq/kg")
        activity_nuc = mat_step.get_activity(units="Bq/kg", by_nuclide=True)
        activity_values.append(activity)
        if isinstance(activity_nuc, dict) and activity_nuc:
            dominant_activity.append(max(activity_nuc, key=activity_nuc.get))
        else:
            dominant_activity.append(None)

    return {
        "name": name,
        "density_g_cm3": mat_def["density"],
        "percent_type": percent_type,
        "elements": mat_def["elements"],
        "waste_disposal_rating": wdr_values,
        "contact_dose_Sv_per_hr": dose_rates,
        "activity_Bq_per_kg": activity_values,
        "dominant_nuclide_wdr": dominant_wdr,
        "dominant_nuclide_dose": dominant_dose,
        "dominant_nuclide_activity": dominant_activity,
    }


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Activation screening of structural materials under DT fusion irradiation"
    )
    parser.add_argument(
        "--test", action="store_true",
        help="Run only 3 test materials (eurofer, SS_316L_N_IG, tungsten)",
    )
    parser.add_argument("--materials", nargs="+", default=None,
                        help="Specific materials to run (keys from MATERIALS dict)")
    parser.add_argument("--chain", default=None,
                        help="Path to depletion chain XML file")
    parser.add_argument("--cross_sections", default=None,
                        help="Path to cross_sections.xml file")
    parser.add_argument("--spectrum", default="hcpb_fw_616.txt",
                        help="Path to FISPACT-style spectrum file (default: hcpb_fw_616.txt)")
    parser.add_argument("--output", default="results/material_data.json",
                        help="Output JSON path")
    args = parser.parse_args()

    if args.cross_sections:
        openmc.config["cross_sections"] = args.cross_sections
    if args.chain:
        openmc.config["chain_file"] = args.chain

    if args.materials:
        materials = args.materials
        for m in materials:
            if m not in MATERIALS:
                print(f"Unknown material: {m}")
                print(f"Available: {', '.join(ALL_MATERIALS)}")
                return
    elif args.test:
        materials = TEST_MATERIALS
    else:
        materials = ALL_MATERIALS

    print(f"Activation screening: {len(materials)} materials")
    print(f"Irradiation: {IRRADIATION_YEARS} FPY @ {WALL_LOADING_MW_M2} MW/m²")
    print(f"Flux: {FLUX:.2e} n/cm²/s")
    print(f"Cooling times: {', '.join(COOLING_LABELS)}\n")

    spectrum, energy_groups = load_spectrum_file(args.spectrum)
    print()

    full_chain = Chain.from_xml(args.chain) if args.chain else Chain.from_xml(
        openmc.config["chain_file"]
    )
    print(f"Chain: {len(full_chain)} nuclides (will reduce per material)\n")

    results = {}
    failed = []

    for i, name in enumerate(materials):
        mat_def = MATERIALS[name]
        t0 = time.time()
        print(
            f"[{i+1}/{len(materials)}] {name} (density={mat_def['density']} g/cm³)...",
            end=" ", flush=True,
        )
        try:
            # Reduce chain to nuclides reachable within 5 transmutation
            # steps from this material's constituent isotopes
            mat_tmp = openmc.Material()
            for element, fraction in mat_def["elements"].items():
                mat_tmp.add_element(element, fraction,
                                    percent_type=mat_def["percent_type"])
            initial_nuclides = list(mat_tmp.get_nuclides())
            reduced_chain = full_chain.reduce(initial_nuclides, level=5)

            data = deplete_material(name, mat_def, spectrum, energy_groups,
                                    chain_file=reduced_chain)
            results[name] = data
            dt = time.time() - t0
            wdr0 = data["waste_disposal_rating"][0]
            dose0 = data["contact_dose_Sv_per_hr"][0]
            print(f"OK ({dt:.1f}s)  shutdown: WDR={wdr0:.2e}  {dose0:.2e} Sv/hr")
        except Exception as exc:
            dt = time.time() - t0
            print(f"FAILED ({dt:.1f}s): {exc}")
            traceback.print_exc()
            failed.append(name)

    # ---- save ----
    output = {
        "metadata": {
            "description": (
                "Systematic activation screening of structural materials "
                "under DT fusion first-wall irradiation"
            ),
            "spectrum": "Approximate DT first-wall, VITAMIN-J-42",
            "wall_loading_MW_m2": WALL_LOADING_MW_M2,
            "availability": AVAILABILITY,
            "total_flux_n_cm2_s": FLUX,
            "irradiation_fpy": IRRADIATION_YEARS,
            "cooling_times": COOLING_LABELS,
            "volume_cm3": 1.0,
            "waste_limits": "StrlSchV_metal_recycling",
            "notes": (
                "Each material uses its real density. Waste disposal rating uses "
                "German StrlSchV 2018 Anlage 4 metal scrap recycling limits "
                "(sum-of-fractions; <1 = meets clearance for metal recycling). "
                "Contact dose in Sv/hr."
            ),
        },
        "materials": results,
    }

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\nSaved: {args.output}")
    if failed:
        print(f"Failed materials: {', '.join(failed)}")
    print(f"Successful: {len(results)}/{len(materials)}")


if __name__ == "__main__":
    main()
