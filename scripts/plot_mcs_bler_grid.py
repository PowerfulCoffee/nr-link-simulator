#!/usr/bin/env python3
"""
Plot C++ (LS-CE + noise est) vs C++ (Ideal/Perfect CE) vs Sionna BLER waterfall curves for MCS 3-27.
Each MCS gets its own subplot in a 5x5 grid.
"""

import csv
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from pathlib import Path

matplotlib.use('Agg')

BASE_DIR = Path("/workspace/nr-link-simulator")
CPP_CSV = BASE_DIR / "results/bler_batch/cpp_bler_mcs3-27.csv"
CPP_PERFECT_CSV = BASE_DIR / "results/bler_batch/cpp_bler_perfect_mcs3-27.csv"
SIONNA_CSV = BASE_DIR / "results/bler_batch/sionna_bler_mcs2-27.csv"
OUT_DIR = BASE_DIR / "results/bler_batch/plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_bler_data(csv_path):
    """Load BLER data keyed by mcs -> sorted list of (snr_db, bler, n_blocks)"""
    data = {}
    if not csv_path.exists():
        print(f"Warning: {csv_path} not found, skipping")
        return data
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mcs = int(row['mcs'])
            snr = float(row['snr_db'])
            bler = float(row['bler'])
            n_blocks = int(row['n_blocks'])
            if mcs not in data:
                data[mcs] = {}
            data[mcs][snr] = (bler, n_blocks)
    return data


def prepare_plot_points(snr_dict):
    """Convert {snr: (bler, n_blocks)} to sorted (snrs, blers) with BLER=0 replaced by floor."""
    snrs = sorted(snr_dict.keys())
    blers = []
    for s in snrs:
        bler, n_blk = snr_dict[s]
        if bler <= 0.0:
            bler = 3.0 / max(n_blk, 30)
        blers.append(max(bler, 1e-4))
    return snrs, blers


def main():
    cpp_data = load_bler_data(CPP_CSV)
    cpp_perfect_data = load_bler_data(CPP_PERFECT_CSV)
    sionna_data = load_bler_data(SIONNA_CSV)

    mcs_range = list(range(3, 28))
    n_mcs = len(mcs_range)
    n_cols = 5
    n_rows = (n_mcs + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 18), sharex=False, sharey=False)
    fig.suptitle('BLER Performance Comparison (AWGN, SISO, 25 PRBs)\n'
                 'C++ LS-CE (DMRS) vs C++ Ideal CE (Perfect CSI) vs Sionna LS-CE\n'
                 'MCS 3-27, PDSCH Table 2 (256QAM)',
                 fontsize=14, fontweight='bold', y=0.995)

    qam_labels = {2: 'QPSK', 4: '16QAM', 6: '64QAM', 8: '256QAM'}

    for idx, mcs in enumerate(mcs_range):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        if mcs <= 9:
            qm = 2
        elif mcs <= 16:
            qm = 4
        elif mcs <= 26:
            qm = 6
        else:
            qm = 8

        all_snrs = []

        if mcs in cpp_perfect_data:
            snrs, blers = prepare_plot_points(cpp_perfect_data[mcs])
            ax.semilogy(snrs, blers, '^-', color='#2ca02c', linewidth=1.8, markersize=5,
                        label='C++ Ideal CE', alpha=0.9, zorder=3)
            all_snrs.extend(snrs)

        if mcs in cpp_data:
            snrs, blers = prepare_plot_points(cpp_data[mcs])
            ax.semilogy(snrs, blers, 'o-', color='#1f77b4', linewidth=1.8, markersize=5,
                        label='C++ LS-CE', alpha=0.85, zorder=2)
            all_snrs.extend(snrs)

        if mcs in sionna_data:
            snrs, blers = prepare_plot_points(sionna_data[mcs])
            ax.semilogy(snrs, blers, 's--', color='#d62728', linewidth=1.5, markersize=4,
                        label='Sionna', alpha=0.85, zorder=1)
            all_snrs.extend(snrs)

        if all_snrs:
            snr_min = min(all_snrs)
            snr_max = max(all_snrs)
            margin = max(1.0, (snr_max - snr_min) * 0.08)
            ax.set_xlim(snr_min - margin, snr_max + margin)

        ax.set_title(f'MCS {mcs} ({qam_labels.get(qm, "QAM")})', fontsize=10, fontweight='bold')
        ax.set_xlabel('SNR (dB)', fontsize=8)
        ax.set_ylabel('BLER', fontsize=8)
        ax.grid(True, which='both', linestyle='--', alpha=0.5)
        ax.set_ylim([5e-4, 3.0])
        ax.tick_params(axis='both', labelsize=7)
        ax.legend(fontsize=6.5, loc='upper right')

        ax.axhline(y=0.1, color='gray', linestyle=':', linewidth=0.7, alpha=0.6)
        ax.axhline(y=0.01, color='gray', linestyle=':', linewidth=0.5, alpha=0.4)

    for idx in range(n_mcs, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.982])
    out_path = OUT_DIR / "bler_mcs3_27_three_way_grid.png"
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved grid plot to: {out_path}")

    out_pdf = OUT_DIR / "bler_mcs3_27_three_way_grid.pdf"
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f"Saved PDF to: {out_pdf}")

    plt.close()


if __name__ == "__main__":
    main()
