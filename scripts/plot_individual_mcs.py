#!/usr/bin/env python3
"""Plot individual SNR vs BLER comparison for each MCS (C++ vs Sionna)."""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MCS_TABLE_1 = [
    (2, 120), (2, 157), (2, 193), (2, 251), (2, 308), (2, 379), (2, 449), (2, 526), (2, 602), (2, 679),
    (4, 340), (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658), (6, 438), (6, 466), (6, 517),
    (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873), (6, 910), (6, 948)
]

def mod_label(qm):
    return {2: 'QPSK', 4: '16QAM', 6: '64QAM', 8: '256QAM'}.get(qm, f'{qm}QAM')

def load_csv(path):
    data = {}
    if not os.path.exists(path):
        print(f"Warning: {path} not found")
        return data
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mcs = int(row['mcs'])
            snr = float(row['snr_db'])
            bler = float(row['bler'])
            n_blk = int(row['n_blocks'])
            n_err = int(row['n_errors'])
            if mcs not in data:
                data[mcs] = {'snr': [], 'bler': [], 'n_blk': [], 'n_err': []}
            data[mcs]['snr'].append(snr)
            data[mcs]['bler'].append(bler)
            data[mcs]['n_blk'].append(n_blk)
            data[mcs]['n_err'].append(n_err)
    for mcs in data:
        for k in data[mcs]:
            data[mcs][k] = np.array(data[mcs][k])
        sort_idx = np.argsort(data[mcs]['snr'])
        for k in data[mcs]:
            data[mcs][k] = data[mcs][k][sort_idx]
    return data

def find_snr_at_bler(snr_arr, bler_arr, target=0.1):
    if len(snr_arr) < 2:
        return None
    bler_work = np.where(bler_arr < 1e-10, 1e-10, bler_arr)
    bler_work = np.clip(bler_work, 1e-10, 1.0)
    if bler_work[0] < target or bler_work[-1] > target:
        return None
    for i in range(len(snr_arr) - 1):
        b1, b2 = bler_work[i], bler_work[i+1]
        if (b1 >= target >= b2) or (b1 <= target <= b2):
            if b1 == b2:
                return snr_arr[i]
            t = (np.log10(target) - np.log10(b1)) / (np.log10(b2) - np.log10(b1))
            return snr_arr[i] + t * (snr_arr[i+1] - snr_arr[i])
    return None

def plot_individual_mcs(cpp_data, sionna_data, mcs, out_dir):
    qm, R_x1024 = MCS_TABLE_1[mcs]
    R = R_x1024 / 1024.0
    mod = mod_label(qm)
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    has_cpp = mcs in cpp_data
    has_sionna = mcs in sionna_data
    
    snr_all = []
    bler_min_plot = 1e-4
    
    if has_cpp:
        snr_c = cpp_data[mcs]['snr']
        bler_c = cpp_data[mcs]['bler'].copy()
        n_blk_c = cpp_data[mcs]['n_blk']
        snr_all.extend(snr_c.tolist())
        bler_c_plot = np.where(bler_c < 1e-10, bler_min_plot, bler_c)
        ax.semilogy(snr_c, bler_c_plot, 'b-o', linewidth=2, markersize=5,
                    markerfacecolor='none', markeredgewidth=1.5, label='C++', zorder=3)
        
        for i, (s, b, nb, ne) in enumerate(zip(snr_c, bler_c, n_blk_c, cpp_data[mcs]['n_err'])):
            if b == 0 or b < 1e-10:
                ax.plot(s, bler_min_plot, 'bv', markersize=6, markerfacecolor='blue', zorder=4)
    
    if has_sionna:
        snr_s = sionna_data[mcs]['snr']
        bler_s = sionna_data[mcs]['bler'].copy()
        n_blk_s = sionna_data[mcs]['n_blk']
        snr_all.extend(snr_s.tolist())
        bler_s_plot = np.where(bler_s < 1e-10, bler_min_plot, bler_s)
        ax.semilogy(snr_s, bler_s_plot, 'r--s', linewidth=2, markersize=5,
                    markerfacecolor='none', markeredgewidth=1.5, label='Sionna', zorder=3)
        
        for i, (s, b, nb, ne) in enumerate(zip(snr_s, bler_s, n_blk_s, sionna_data[mcs]['n_err'])):
            if b == 0 or b < 1e-10:
                ax.plot(s, bler_min_plot, 'rv', markersize=6, markerfacecolor='red', zorder=4)
    
    ax.axhline(y=0.1, color='gray', linestyle=':', linewidth=1.0, alpha=0.7, label='BLER=0.1')
    
    ax.set_xlabel('SNR (dB)', fontsize=13)
    ax.set_ylabel('BLER', fontsize=13)
    ax.set_title(f'MCS {mcs} ({mod}, R={R:.3f}) — AWGN, LS-CE + LMMSE, 25 PRB', fontsize=14)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, which='both', alpha=0.3)
    
    if snr_all:
        snr_min = min(snr_all) - 1
        snr_max = max(snr_all) + 1
        ax.set_xlim([snr_min, snr_max])
    ax.set_ylim([bler_min_plot * 0.5, 1.5])
    
    cpp_snr10 = find_snr_at_bler(cpp_data[mcs]['snr'], cpp_data[mcs]['bler'], 0.1) if has_cpp else None
    sio_snr10 = find_snr_at_bler(sionna_data[mcs]['snr'], sionna_data[mcs]['bler'], 0.1) if has_sionna else None
    
    info_text = []
    if cpp_snr10 is not None:
        info_text.append(f'C++  @ BLER=0.1: SNR={cpp_snr10:.2f} dB')
    if sio_snr10 is not None:
        info_text.append(f'Sionna @ BLER=0.1: SNR={sio_snr10:.2f} dB')
    if cpp_snr10 is not None and sio_snr10 is not None:
        diff = sio_snr10 - cpp_snr10
        info_text.append(f'Δ = {diff:+.2f} dB')
    
    if info_text:
        ax.text(0.02, 0.02, '\n'.join(info_text), transform=ax.transAxes,
                fontsize=10, verticalalignment='bottom', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    out_path = os.path.join(out_dir, f'mcs{mcs:02d}_bler_comparison.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")
    return cpp_snr10, sio_snr10

def main():
    out_dir = "/workspace/nr-link-simulator/results/bler_batch/plots"
    os.makedirs(out_dir, exist_ok=True)
    
    cpp_csv = "/workspace/nr-link-simulator/results/bler_batch/cpp_bler_mcs2-27.csv"
    sio_csv = "/workspace/nr-link-simulator/results/bler_batch/sionna_bler_mcs2-27.csv"
    
    print("Loading C++ data...")
    cpp_data = load_csv(cpp_csv)
    print(f"  Loaded {len(cpp_data)} MCS values from C++")
    
    print("Loading Sionna data...")
    sio_data = load_csv(sio_csv)
    print(f"  Loaded {len(sio_data)} MCS values from Sionna")
    
    mcs_list = sorted(set(list(cpp_data.keys()) + list(sio_data.keys())))
    print(f"\nPlotting {len(mcs_list)} MCS comparison charts...")
    
    summary = []
    for mcs in mcs_list:
        qm, R_x1024 = MCS_TABLE_1[mcs]
        R = R_x1024 / 1024.0
        mod = mod_label(qm)
        cpp_s, sio_s = plot_individual_mcs(cpp_data, sio_data, mcs, out_dir)
        summary.append((mcs, mod, R, cpp_s, sio_s))
    
    print("\n" + "="*70)
    print(f"{'MCS':>4} {'Mod':>6} {'Rate':>6} {'C++ SNR@0.1':>12} {'Sionna SNR@0.1':>14} {'Δ(dB)':>8}")
    print("-"*70)
    for mcs, mod, R, cs, ss in summary:
        cs_str = f"{cs:.2f}" if cs is not None else "N/A"
        ss_str = f"{ss:.2f}" if ss is not None else "N/A"
        if cs is not None and ss is not None:
            d = ss - cs
            d_str = f"{d:+.2f}"
        else:
            d_str = "N/A"
        print(f"{mcs:4d} {mod:>6} {R:6.3f} {cs_str:>12} {ss_str:>14} {d_str:>8}")
    print("="*70)
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 16))
    groups = [
        ('QPSK (MCS 2-9)', [m for m in range(2, 10)]),
        ('16QAM (MCS 10-16)', [m for m in range(10, 17)]),
        ('64QAM (MCS 17-27)', [m for m in range(17, 28)]),
    ]
    cmap = plt.cm.tab20
    
    for ax, (title, mcs_sub) in zip(axes, groups):
        for idx, mcs in enumerate(mcs_sub):
            color = cmap(idx / max(len(mcs_sub)-1, 1))
            qm, Rx = MCS_TABLE_1[mcs]
            R = Rx/1024.0
            if mcs in cpp_data:
                s = cpp_data[mcs]['snr']
                b = np.where(cpp_data[mcs]['bler'] < 1e-10, 1e-4, cpp_data[mcs]['bler'])
                ax.semilogy(s, b, '-', color=color, linewidth=1.8, label=f'MCS{mcs} C++')
            if mcs in sio_data:
                s = sio_data[mcs]['snr']
                b = np.where(sio_data[mcs]['bler'] < 1e-10, 1e-4, sio_data[mcs]['bler'])
                ax.semilogy(s, b, '--', color=color, linewidth=1.8, label=f'MCS{mcs} Sio' if idx==0 else None)
        ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5)
        ax.set_xlabel('SNR (dB)')
        ax.set_ylabel('BLER')
        ax.set_title(title)
        ax.legend(fontsize=7, ncol=3, loc='lower left')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([1e-4, 2])
    
    plt.tight_layout()
    overview_path = os.path.join(out_dir, 'overview_all_mcs.png')
    fig.savefig(overview_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nOverview saved to {overview_path}")
    print(f"Individual plots saved to {out_dir}/")

if __name__ == "__main__":
    main()
