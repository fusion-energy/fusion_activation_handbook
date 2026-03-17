#!/usr/bin/env python3
"""Plot a FISPACT-style neutron spectrum file.

Usage:
    python plot_spectrum.py
    python plot_spectrum.py --spectrum hcpb_fw_616.txt
"""

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from generate_data import load_spectrum_file


def main():
    parser = argparse.ArgumentParser(description="Plot neutron spectrum")
    parser.add_argument("--spectrum", default="hcpb_fw_616.txt",
                        help="Path to FISPACT-style spectrum file")
    parser.add_argument("--output", default="results/hcpb_spectrum.png",
                        help="Output image path")
    args = parser.parse_args()

    flux, energy_group_name = load_spectrum_file(args.spectrum)

    import openmc.mgxs
    boundaries = np.array(openmc.mgxs.GROUP_STRUCTURES[energy_group_name])
    midpoints = np.sqrt(boundaries[:-1] * boundaries[1:])

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.loglog(midpoints / 1e6, flux, linewidth=0.8)
    ax.set_xlabel("Energy [MeV]")
    ax.set_ylabel("Flux [n/cm²/s]")
    ax.set_title("HCPB First-Wall Neutron Spectrum (616 groups, ascending energy)")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1e-11, 25)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    plt.close()
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
