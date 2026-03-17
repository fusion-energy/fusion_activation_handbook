# Fusion Activation Screening

Systematic activation screening of elements (Z=1–83) under DT fusion
first-wall irradiation using OpenMC's transport-free material depletion.

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

# 3. Generate periodic table heatmaps and top-N summary plots
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
spectrum, then tracks activity and dose through 10 cooling times
(shutdown to 100 years). Outputs `results/element_data.json`.

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

```
python generate_data.py --chain /home/jon/nuclear_data/chain-endf-b8.0.xml --cross_sections /home/jon/nuclear_data/endf-b8.0-hdf5/cross_sections.xml --spectrum hcpb_fw_616.txt
```

### `plot_periodic_table.py`

Reads `element_data.json` and generates:
- Individual periodic table heatmaps for each metric at each cooling time
- Top-N summary plots (elements ever in worst N at any cooling time)

```
  --input               Input JSON (default: results/element_data.json)
  --output-dir          Output directory (default: results)
  --metric              activity, dose, or both (default: both)
  --top-n               Number of worst elements to flag (default: 10)
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
  periodic_table_activity_*.png      # Activity heatmaps per cooling time
  periodic_table_dose_*.png          # Dose heatmaps per cooling time
  top10_activity.png                 # Summary: worst elements by activity
  top10_dose.png                     # Summary: worst elements by dose
  top10_combined.png                 # Summary: worst elements by either
```
