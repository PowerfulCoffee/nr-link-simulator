#!/usr/bin/env python3
"""
Run C++ vs Sionna 4x4 MIMO TDL-A BLER comparison
Config: 30PRB, 4Rx, rank4, DMRS Type1 add-pos=1 (symbols 2&7), TDL-A DS=100ns
MCS Table 2, comparing QPSK/16QAM/64QAM/256QAM with Ideal CE and LS CE
"""
import os
import sys
import subprocess
import time
from pathlib import Path
import csv

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent.parent
BUILD_DIR = ROOT / "build"
SCRIPTS_DIR = ROOT / "scripts"
RESULTS_DIR = ROOT / "results" / "tdl_crosscheck_4x4"
CPP_BIN = BUILD_DIR / "examples" / "pdsch_bler_tdl_mimo"
SIONNA_SCRIPT = SCRIPTS_DIR / "sionna_tdl_crosscheck.py"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MCS_CONFIGS = [
    {"mcs": 4,  "mod": "QPSK",   "qam": 2, "snr_start": -4, "snr_end": 16, "snr_step": 2},
    {"mcs": 10, "mod": "16QAM",  "qam": 4, "snr_start": 4,  "snr_end": 24, "snr_step": 2},
    {"mcs": 16, "mod": "64QAM",  "qam": 6, "snr_start": 10, "snr_end": 30, "snr_step": 2},
    {"mcs": 24, "mod": "256QAM", "qam": 8, "snr_start": 16, "snr_end": 36, "snr_step": 2},
]

def run_cpp(mcs, snr_start, snr_end, snr_step, perfect_csi, max_blocks=200, target_errors=50, out_file=None):
    if out_file is None:
        ce_type = "ideal" if perfect_csi else "ls"
        out_file = RESULTS_DIR / f"cpp_{ce_type}_mcs{mcs}.csv"
    if out_file.exists() and out_file.stat().st_size > 50:
        print(f"  [C++ {ce_type}] Already exists: {out_file.name}")
        return out_file

    args = [
        str(CPP_BIN),
        str(mcs),
        str(snr_start),
        str(snr_end),
        str(snr_step),
        "1" if perfect_csi else "0",
        str(max_blocks),
        str(target_errors),
        "30",
        "4",
        "4",
        "1",
        "100e-9",
        "5.0",
        "1",
        "0",
        "2",
        str(out_file),
    ]
    print(f"  [C++ {'ideal' if perfect_csi else 'ls'}] Running MCS{mcs}...")
    t0 = time.time()
    result = subprocess.run(args, cwd=str(BUILD_DIR), capture_output=True, text=True, timeout=3600)
    dt = time.time() - t0
    if result.returncode != 0:
        print(f"  [C++ ERROR] {result.stderr[-500:]}")
        return None
    print(f"  [C++ {'ideal' if perfect_csi else 'ls'}] Done in {dt:.1f}s -> {out_file.name}")
    return out_file

def run_sionna(mcs, perfect_csi, max_blocks=100, out_file=None):
    if out_file is None:
        ce_type = "ideal" if perfect_csi else "ls"
        out_file = RESULTS_DIR / f"sionna_{ce_type}_mcs{mcs}.csv"
    if out_file.exists() and out_file.stat().st_size > 50:
        print(f"  [Sionna {ce_type}] Already exists: {out_file.name}")
        return out_file

    args = [sys.executable, str(SIONNA_SCRIPT), str(mcs), "1" if perfect_csi else "0", str(out_file)]
    print(f"  [Sionna {'ideal' if perfect_csi else 'ls'}] Running MCS{mcs}...")
    t0 = time.time()
    result = subprocess.run(args, cwd=str(ROOT), capture_output=True, text=True, timeout=7200)
    dt = time.time() - t0
    if result.returncode != 0:
        print(f"  [Sionna ERROR] {result.stderr[-1000:]}")
        return None
    print(f"  [Sionna {'ideal' if perfect_csi else 'ls'}] Done in {dt:.1f}s -> {out_file.name}")
    return out_file

def read_csv(filepath):
    if not filepath or not Path(filepath).exists():
        return [], [], []
    snrs, blers, blocks = [], [], []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                snrs.append(float(row['SINR_dB']))
                blers.append(float(row['BLER']))
                blocks.append(int(row['Blocks']))
            except (ValueError, KeyError):
                continue
    return np.array(snrs), np.array(blers), np.array(blocks)

def plot_all_results():
    fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharey=True)
    if not isinstance(axes, np.ndarray):
        axes = [axes]

    markers = {'cpp_ideal': 'o-', 'cpp_ls': 's-', 'sionna_ideal': '^--', 'sionna_ls': 'v--'}
    colors = {'cpp_ideal': '#1f77b4', 'cpp_ls': '#ff7f0e', 'sionna_ideal': '#2ca02c', 'sionna_ls': '#d62728'}
    labels = {'cpp_ideal': 'C++ Ideal CE', 'cpp_ls': 'C++ LS CE', 'sionna_ideal': 'Sionna Ideal CE', 'sionna_ls': 'Sionna LS CE'}

    for ax, cfg in zip(axes, MCS_CONFIGS):
        mcs = cfg['mcs']
        mod = cfg['mod']
        for key in ['cpp_ideal', 'cpp_ls', 'sionna_ideal', 'sionna_ls']:
            ce_type = key.split('_')[1]
            is_perfect = ce_type == 'ideal'
            prefix = key.split('_')[0]
            filepath = RESULTS_DIR / f"{prefix}_{ce_type}_mcs{mcs}.csv"
            snrs, blers, _ = read_csv(filepath)
            if len(snrs) > 0:
                ax.semilogy(snrs, blers, markers[key], color=colors[key], label=labels[key],
                           markersize=5, linewidth=1.5, alpha=0.8)

        ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5)
        ax.set_xlabel('SNR (dB)')
        ax.set_title(f'{mod} (MCS{mcs}, Table 2)\n30PRB 4x4 Rank4, TDL-A DS=100ns')
        ax.grid(True, which='both', alpha=0.3)
        ax.set_ylim([1e-3, 1.05])
        ax.legend(fontsize=7, loc='lower left')

    axes[0].set_ylabel('BLER')
    plt.tight_layout()

    plot_file = RESULTS_DIR / "bler_comparison_4x4_tdl.png"
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.savefig(RESULTS_DIR / "bler_comparison_4x4_tdl.pdf", bbox_inches='tight')
    print(f"\nPlot saved to {plot_file}")
    return plot_file

def main():
    print("=" * 60)
    print("NR PDSCH BLER: C++ vs Sionna 4x4 MIMO TDL-A Comparison")
    print("=" * 60)

    if not CPP_BIN.exists():
        print(f"ERROR: C++ binary not found at {CPP_BIN}")
        print("Please build first: cd build && cmake .. && make -j")
        sys.exit(1)

    all_tasks = []
    for cfg in MCS_CONFIGS:
        all_tasks.append(('cpp', cfg, True))
        all_tasks.append(('cpp', cfg, False))
        all_tasks.append(('sionna', cfg, True))
        all_tasks.append(('sionna', cfg, False))

    results = {}
    for i, (engine, cfg, perfect) in enumerate(all_tasks):
        mcs = cfg['mcs']
        mod = cfg['mod']
        ce_type = "ideal" if perfect else "ls"
        print(f"\n[{i+1}/{len(all_tasks)}] {engine.upper()} {mod} MCS{mcs} {ce_type} CE")

        if engine == 'cpp':
            f = run_cpp(mcs, cfg['snr_start'], cfg['snr_end'], cfg['snr_step'], perfect,
                       max_blocks=200, target_errors=50)
        else:
            f = run_sionna(mcs, perfect, max_blocks=100)
        results[(engine, ce_type, mcs)] = f

    print("\n" + "=" * 60)
    print("Generating comparison plot...")
    print("=" * 60)
    plot_file = plot_all_results()

    print("\nResults summary:")
    for cfg in MCS_CONFIGS:
        mcs = cfg['mcs']
        mod = cfg['mod']
        print(f"\n--- {mod} (MCS{mcs}) ---")
        for key in ['cpp_ideal', 'cpp_ls', 'sionna_ideal', 'sionna_ls']:
            prefix, ce_type = key.split('_')
            filepath = RESULTS_DIR / f"{prefix}_{ce_type}_mcs{mcs}.csv"
            snrs, blers, _ = read_csv(filepath)
            if len(snrs) > 0:
                bler10_idx = np.argmin(np.abs(blers - 0.1)) if np.any(blers < 0.1) else None
                if bler10_idx is not None and blers[bler10_idx] < 0.5:
                    print(f"  {key:15s}: BLER=0.1 @ ~{snrs[bler10_idx]:.1f} dB")
                else:
                    min_bler = blers.min() if len(blers) > 0 else 1.0
                    print(f"  {key:15s}: min BLER={min_bler:.4f} at {snrs[-1]:.1f} dB")
            else:
                print(f"  {key:15s}: NO DATA")

if __name__ == "__main__":
    main()
