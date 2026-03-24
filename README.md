# Fusion Activation Screening

Systematic activation screening of elements (Z=1–83) under DT fusion
first-wall irradiation using OpenMC's transport-free material depletion.

Each element is irradiated for 40 full-power years at 1 MW/m² wall loading
with the HCPB first-wall neutron spectrum, then cooled from shutdown to
100 years. Two metrics are computed:

- **Waste disposal rating** (Fetter limits) — determines whether irradiated
  material qualifies as low-level waste (rating < 1.0). Plotted at
  waste-relevant timescales: 1 year, 10 years, 100 years.
- **Contact dose rate** (Sv/hr) — determines when maintenance personnel can
  access irradiated components. Plotted at maintenance-relevant timescales:
  1 day, 1 week, 1 month.

Results are visualized as color-coded periodic tables.

See the full write-up in the [report (PDF)](report.pdf).

## Prerequisites

- Python 3.10+
- OpenMC >= 0.15.3 (with shared library support)
- numpy, matplotlib
- ENDF/B-VIII.0 cross sections and depletion chain
- HCPB-FW spectrum file (`hcpb_fw_616.txt`, included)

## Quick start

```bash
# 1. Run 10 test elements (Fe, Cr, W, Ni, Co, Mo, Nb, V, Ti, Ta)
python generate_data.py --test \
  --chain /path/to/chain-endf-b8.0.xml \
  --cross_sections /path/to/endf-b8.0-hdf5/cross_sections.xml

# 2. Plot the input spectrum
python plot_spectrum.py

# 3. Generate periodic table heatmaps and summary plots
python plot_periodic_table.py
```

## Full run (all 81 elements)

```bash
python generate_data.py \
  --chain /path/to/chain-endf-b8.0.xml \
  --cross_sections /path/to/endf-b8.0-hdf5/cross_sections.xml
```

## Scripts

### `generate_data.py`

Depletes each element for 40 full-power years under the HCPB first-wall
spectrum, then tracks waste disposal rating and contact dose rate through
10 cooling times (shutdown to 100 years). Outputs `results/element_data.json`.

```
usage: generate_data.py [-h] [--test] [--elements ELEMENTS [ELEMENTS ...]]
                        [--chain CHAIN] [--cross_sections CROSS_SECTIONS]
                        [--spectrum SPECTRUM] [--output OUTPUT]

  --test                Run only 10 test elements
  --elements            Specific elements to run (e.g. --elements Fe Cr W)
  --chain               Path to depletion chain XML file
  --cross_sections      Path to cross_sections.xml file
  --spectrum            Path to FISPACT-style spectrum file (default: hcpb_fw_616.txt)
  --output              Output JSON path (default: results/element_data.json)
```

### `plot_periodic_table.py`

Reads `element_data.json` and generates:
- Waste disposal rating heatmaps at 1 year, 10 years, 100 years
- Contact dose rate heatmaps at 1 day, 1 week, 1 month
- 6 summary plots identifying problematic elements at each timescale
- 1 combined summary showing elements problematic for waste, dose, or both

```
  --input               Input JSON (default: results/element_data.json)
  --output-dir          Output directory (default: results)
  --metric              wdr, dose, or both (default: both)
  --top-n               Number of worst elements to flag as problematic (default: 10)
```

### `plot_spectrum.py`

Plots the input neutron spectrum.

```
  --spectrum            Path to spectrum file (default: hcpb_fw_616.txt)
  --output              Output image path (default: results/hcpb_spectrum.png)
```

## Output

All results go to `results/`:

```
results/
  element_data.json                  # Raw data for all elements
  hcpb_spectrum.png                  # Input spectrum plot
  periodic_table_wdr_*.png           # Waste disposal rating heatmaps (1yr, 10yr, 100yr)
  periodic_table_dose_*.png          # Contact dose rate heatmaps (1d, 1wk, 1mo)
  top10_waste_short.png              # Problematic elements: waste at 1 year
  top10_waste_medium.png             # Problematic elements: waste at 10 years
  top10_waste_long.png               # Problematic elements: waste at 100 years
  top10_dose_short.png               # Problematic elements: dose at 1 day
  top10_dose_medium.png              # Problematic elements: dose at 1 week
  top10_dose_long.png                # Problematic elements: dose at 1 month
  top10_combined.png                 # Combined: problematic for waste, dose, or both
```
