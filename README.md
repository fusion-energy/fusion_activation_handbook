# WORK IN PROGRESS

# fusion_activation_handbook

Code for automated production of a fusion materials activation handbook using OpenMC.

The handbook contains the following simulated values on a wide range of elements and materials exposed to neutron irradiation.
- Simulated activity
- Atom inventory
- Dose
- Decay heat

Download the latest PDF report [here](link to latest release)

# To reproduce the handbook locally

## 1. Install OpenMC and other dependencies

From a conda environment you can install the OpenMC
```bash
conda install -c conda-forge openmc==0.13.3
pip install openmc_depletion_plotter
```

## 2. Install the nuclear data
```bash
pip install openmc_data
download_nndc -r b8.0
download_nndc_chain -r b8.0
```

## 3. Run the simulation scripts and produce plots
```bash
python simulate_irradiation.py
```

At this point the png files for each figure have been created in the figs folder

## 5. Produce the report in MarkDown

```bash
python produce_markdown_report.py -nuc_data=ENDF/B8.0
```
This takes a report template and fills in details for the version of OpenMC used, nuclear data used in MarkDown format.

## 6. Produce the report in PDF
To create the report in PDF this can be made by converting the MarkDown format but we need to install Latex and Pandoc.
```bash
apt-get install -y texlive-latex-base texlive-fonts-recommended
apt-get install -y pandoc
pandoc activation_handbook.md activation_handbook.pdf
```
