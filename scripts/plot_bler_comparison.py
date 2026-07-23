#!/usr/bin/env python3
"""
Plot SNR vs BLER comparison for 4 MCS levels in a 2x2 grid
Methods: Ideal CE (Perfect), C++ LS-Doppler, Sionna LS (no Doppler)
"""
import csv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 9,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
})

RESULTS_DIR = Path("/workspace/nr-link-simulator/scripts/results/tdl_doppler")

MCS_INFO = {
    3:  {"name": "QPSK (MCS 3)",  "mod": "QPSK",  "R": 0.27},
    10: {"name": "16QAM (MCS 10)", "mod": "16QAM", "R": 0.50},
    17: {"name": "64QAM (MCS 17)", "mod": "64QAM", "R": 0.62},
    27: {"name": "64QAM high-rate (MCS 27)", "mod": "64QAM (high R)", "R": 0.89},
}

METHODS = [
    ("cpp_perfect",        "C++ Ideal CE",        "-",  "o", "#1f77b4", 1.5, 6),
    ("cpp_ls_doppler",     "C++ LS + Doppler",    "--", "s", "#2ca02c", 1.5, 6),
    ("sionna_perfect",     "Sionna Ideal CE",     ":",  "^", "#1f77b4", 2.0, 6),
    ("sionna_ls_nodoppler","Sionna LS (no Doppler)", "-.", "v", "#d62728", 1.5, 6),
]

def load_csv(filepath):
    snr, bler = [], []
    if not filepath.exists():
        return None, None
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                s = float(row['SINR_dB'])
                b = float(row['BLER'])
                n = int(row['Blocks'])
                if n >= 30 and b <= 1.0:
                    snr.append(s)
                    bler.append(b)
            except:
                pass
    return np.array(snr), np.array(bler)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

mcs_list = [3, 10, 17, 27]

plot_data = {}
for mcs in mcs_list:
    plot_data[mcs] = {}
    for method_key, _, _, _, _, _, _ in METHODS:
        f = RESULTS_DIR / f"{method_key}_mcs{mcs}.csv"
        snr, bler = load_csv(f)
        if snr is not None and len(snr) > 0:
            plot_data[mcs][method_key] = (snr, bler)

for idx, mcs in enumerate(mcs_list):
    ax = axes[idx]
    info = MCS_INFO[mcs]
    
    for method_key, label, ls, marker, color, lw, ms in METHODS:
        if method_key in plot_data[mcs]:
            snr, bler = plot_data[mcs][method_key]
            if len(snr) > 0:
                ax.semilogy(snr, bler, label=label, linestyle=ls, marker=marker,
                           color=color, linewidth=lw, markersize=ms, markevery=max(1, len(snr)//10))
    
    ax.set_xlabel('SNR (dB)')
    ax.set_ylabel('BLER')
    ax.set_title(f"{info['name']}, R≈{info['R']}", fontweight='bold')
    ax.grid(True, which='both', alpha=0.3, linestyle='-')
    ax.set_ylim([0.008, 1.0])
    ax.axhline(y=0.1, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.text(ax.get_xlim()[0]+0.5 if ax.get_xlim()[0] < ax.get_xlim()[1] else 0, 0.12, 'BLER=0.1', 
            fontsize=8, color='gray', alpha=0.7)
    ax.legend(loc='lower left', framealpha=0.9)

fig.suptitle(
    'TDL-A Channel Performance Comparison\n'
    'DS=100ns, $f_d$=70Hz, 3PRB SISO, DMRS Type1 add-pos=1 (dual DMRS)\n'
    'Ideal CE vs C++ LS with Doppler estimation/compensation vs Sionna LS (linear interp only)',
    fontsize=13, fontweight='bold', y=1.02
)

plt.tight_layout()
out_path = "/workspace/nr-link-simulator/scripts/results/tdl_doppler/bler_comparison_grid.png"
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"Plot saved to {out_path}")

print("\n" + "="*70)
print("PERFORMANCE SUMMARY (SNR required for BLER = 0.1)")
print("="*70)
for mcs in mcs_list:
    info = MCS_INFO[mcs]
    print(f"\n{info['name']} (R≈{info['R']}):")
    for method_key, label, _, _, _, _, _ in METHODS:
        if method_key in plot_data[mcs]:
            snr, bler = plot_data[mcs][method_key]
            if len(snr) > 1:
                bler = np.array(bler)
                snr = np.array(snr)
                sort_idx = np.argsort(snr)
                snr_s = snr[sort_idx]
                bler_s = bler[sort_idx]
                if bler_s[0] > 0.1 and bler_s[-1] < 0.1:
                    snr_at_01 = np.interp(0.1, bler_s[::-1], snr_s[::-1])
                    print(f"  {label:30s}: {snr_at_01:5.1f} dB")
                elif bler_s[-1] >= 0.1:
                    above_idx = np.where(bler_s < 0.15)[0]
                    if len(above_idx) > 0:
                        snr_est = snr_s[above_idx[0]]
                        print(f"  {label:30s}: ~{snr_est:4.1f} dB (BLER≈{bler_s[above_idx[0]]:.2f} at highest SNR)")
                    else:
                        print(f"  {label:30s}: >{snr_s[-1]:.1f} dB (BLER={bler_s[-1]:.3f})")
