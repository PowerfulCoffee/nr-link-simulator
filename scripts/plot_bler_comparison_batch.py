#!/usr/bin/env python3
"""Plot BLER vs SNR comparison: C++ vs Sionna for MCS 2-27."""

import csv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

MCS_TABLE_1 = [
    (2, 120), (2, 157), (2, 193), (2, 251), (2, 308), (2, 379), (2, 449), (2, 526), (2, 602), (2, 679),
    (4, 340), (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658), (6, 438), (6, 466), (6, 517),
    (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873), (6, 910), (6, 948)
]

RESULTS_DIR = "/workspace/nr-link-simulator/results/bler_batch"
OUT_DIR = "/workspace/nr-link-simulator/results/bler_batch"

def load_csv(filepath):
    data = {}
    if not os.path.exists(filepath):
        return data
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mcs = int(row['mcs'])
            snr = float(row['snr_db'])
            bler = float(row['bler'])
            n_blk = int(row['n_blocks'])
            n_err = int(row['n_errors'])
            if mcs not in data:
                data[mcs] = {'snr': [], 'bler': [], 'n_blocks': [], 'n_errors': []}
            data[mcs]['snr'].append(snr)
            data[mcs]['bler'].append(bler)
            data[mcs]['n_blocks'].append(n_blk)
            data[mcs]['n_errors'].append(n_err)
    for mcs in data:
        for k in data[mcs]:
            data[mcs][k] = np.array(data[mcs][k])
    return data

def mod_label(qm):
    if qm == 2: return 'QPSK'
    elif qm == 4: return '16QAM'
    elif qm == 6: return '64QAM'
    elif qm == 8: return '256QAM'
    return f'{qm}QAM'

def plot_all_mcs(cpp_data, sionna_data, output_path):
    mcs_list = list(range(2, 28))
    
    qpsk_mcs = [m for m in mcs_list if MCS_TABLE_1[m][0] == 2]
    qam16_mcs = [m for m in mcs_list if MCS_TABLE_1[m][0] == 4]
    qam64_mcs = [m for m in mcs_list if MCS_TABLE_1[m][0] == 6]
    
    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.25)
    
    groups = [
        ('QPSK (MCS 2-9)', qpsk_mcs, 0, 0),
        ('16QAM (MCS 10-16)', qam16_mcs, 0, 1),
        ('64QAM (MCS 17-27)', qam64_mcs, 0, 2),
        ('QPSK (MCS 2-9) Zoom', qpsk_mcs, 1, 0),
        ('16QAM (MCS 10-16) Zoom', qam16_mcs, 1, 1),
        ('64QAM (MCS 17-27) Zoom', qam64_mcs, 1, 2),
    ]
    
    cmap = plt.cm.viridis
    
    for title, mcs_subset, row, col in groups:
        ax = fig.add_subplot(gs[row, col])
        n_mcs = len(mcs_subset)
        
        for idx, mcs in enumerate(mcs_subset):
            qm, R_x1024 = MCS_TABLE_1[mcs]
            R = R_x1024 / 1024.0
            color = cmap(idx / max(n_mcs - 1, 1))
            label = f'MCS {mcs} ({mod_label(qm)}, R={R:.2f})'
            
            if mcs in cpp_data:
                snr = cpp_data[mcs]['snr']
                bler = cpp_data[mcs]['bler']
                valid = bler > 0
                if np.any(bler < 1e-10):
                    bler_plot = np.where(bler < 1e-10, 1e-5, bler)
                else:
                    bler_plot = bler
                ax.semilogy(snr, bler_plot, '-', color=color, linewidth=1.5, label=f'C++ {label}' if idx == 0 else None)
                ax.semilogy(snr, bler_plot, 'o', color=color, markersize=3, markerfacecolor='none')
            
            if mcs in sionna_data:
                snr_s = sionna_data[mcs]['snr']
                bler_s = sionna_data[mcs]['bler']
                if np.any(bler_s < 1e-10):
                    bler_s_plot = np.where(bler_s < 1e-10, 1e-5, bler_s)
                else:
                    bler_s_plot = bler_s
                ax.semilogy(snr_s, bler_s_plot, '--', color=color, linewidth=1.5, label=f'Sionna {label}' if idx == 0 else None)
                ax.semilogy(snr_s, bler_s_plot, 's', color=color, markersize=3, markerfacecolor='none')
        
        ax.axhline(y=0.1, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
        ax.set_xlabel('SNR (dB)')
        ax.set_ylabel('BLER')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.set_ylim([1e-4, 1.0])
        
        if 'Zoom' in title:
            ax.set_ylim([1e-3, 1.0])
    
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color='black', linestyle='-', marker='o', markersize=5, label='C++ (LS CE + LMMSE)'),
        Line2D([0], [0], color='black', linestyle='--', marker='s', markersize=5, label='Sionna (LS CE + LMMSE)'),
    ]
    fig.legend(handles=custom_lines, loc='upper center', bbox_to_anchor=(0.5, 0.98), ncol=2, fontsize=12)
    
    fig.suptitle('NR PDSCH BLER vs SNR: C++ nr-link-simulator vs Sionna (AWGN, 25 PRB, 1x1 SISO, 1dB step)', 
                 fontsize=14, y=1.0)
    
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Overview plot saved to {output_path}")

def plot_per_mcs(cpp_data, sionna_data, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    for mcs in range(2, 28):
        qm, R_x1024 = MCS_TABLE_1[mcs]
        R = R_x1024 / 1024.0
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        has_data = False
        
        if mcs in cpp_data:
            snr = cpp_data[mcs]['snr']
            bler = cpp_data[mcs]['bler'].copy()
            n_blk = cpp_data[mcs]['n_blocks']
            n_err = cpp_data[mcs]['n_errors']
            bler[bler < 1e-10] = 1e-5
            ax.semilogy(snr, bler, 'b-o', linewidth=2, markersize=5, label='C++ (LS CE + LMMSE)')
            for i, (s, b, nb, ne) in enumerate(zip(snr, bler, n_blk, n_err)):
                if ne > 0 and nb > 0:
                    ax.annotate(f'{ne}/{nb}', (s, b), textcoords="offset points", 
                               xytext=(5, 5), fontsize=7, color='blue')
            has_data = True
        
        if mcs in sionna_data:
            snr_s = sionna_data[mcs]['snr']
            bler_s = sionna_data[mcs]['bler'].copy()
            n_blk_s = sionna_data[mcs]['n_blocks']
            n_err_s = sionna_data[mcs]['n_errors']
            bler_s[bler_s < 1e-10] = 1e-5
            ax.semilogy(snr_s, bler_s, 'r--s', linewidth=2, markersize=5, label='Sionna (LS CE + LMMSE)')
            for i, (s, b, nb, ne) in enumerate(zip(snr_s, bler_s, n_blk_s, n_err_s)):
                if ne > 0 and nb > 0:
                    ax.annotate(f'{ne}/{nb}', (s, b), textcoords="offset points",
                               xytext=(5, -10), fontsize=7, color='red')
            has_data = True
        
        if not has_data:
            plt.close(fig)
            continue
        
        ax.axhline(y=0.1, color='gray', linestyle=':', linewidth=1, alpha=0.7, label='BLER=0.1')
        ax.set_xlabel('SNR (dB)', fontsize=12)
        ax.set_ylabel('BLER', fontsize=12)
        ax.set_title(f'MCS {mcs}: {mod_label(qm)}, R={R:.3f} (AWGN, 25 PRB, 1x1 SISO)', fontsize=13)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_ylim([1e-4, 1.0])
        
        out_path = os.path.join(output_dir, f'bler_mcs{mcs:02d}.png')
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    
    print(f"Per-MCS plots saved to {output_dir}")

def compute_gap_summary(cpp_data, sionna_data):
    print("\n=== BLER=0.1 SNR Gap Summary (C++ - Sionna) ===")
    print(f"{'MCS':>4} {'Mod':>6} {'R':>6} {'C++ SNR':>8} {'Sionna SNR':>10} {'Gap(dB)':>8}")
    print("-" * 50)
    
    gaps = []
    for mcs in range(2, 28):
        qm, R_x1024 = MCS_TABLE_1[mcs]
        R = R_x1024 / 1024.0
        
        def find_snr_at_bler(data, mcs, target=0.1):
            if mcs not in data:
                return None
            snr = data[mcs]['snr']
            bler = data[mcs]['bler']
            idx = np.argsort(snr)
            snr_sorted = snr[idx]
            bler_sorted = bler[idx]
            for i in range(len(bler_sorted) - 1):
                b1, b2 = bler_sorted[i], bler_sorted[i+1]
                s1, s2 = snr_sorted[i], snr_sorted[i+1]
                if (b1 >= target and b2 <= target) or (b1 <= target and b2 >= target):
                    if abs(b1 - b2) < 1e-12:
                        return s1
                    t = (np.log10(target) - np.log10(max(b1, 1e-10))) / (np.log10(max(b2, 1e-10)) - np.log10(max(b1, 1e-10)))
                    return s1 + t * (s2 - s1)
            return None
        
        snr_cpp = find_snr_at_bler(cpp_data, mcs)
        snr_sio = find_snr_at_bler(sionna_data, mcs)
        
        if snr_cpp is not None and snr_sio is not None:
            gap = snr_cpp - snr_sio
            gaps.append(gap)
            print(f"{mcs:4d} {mod_label(qm):>6} {R:6.3f} {snr_cpp:8.2f} {snr_sio:10.2f} {gap:+8.2f}")
        else:
            print(f"{mcs:4d} {mod_label(qm):>6} {R:6.3f} {'N/A':>8} {'N/A':>10} {'N/A':>8}")
    
    if gaps:
        print("-" * 50)
        print(f"Mean gap: {np.mean(gaps):+.2f} dB, Max gap: {np.max(np.abs(gaps)):.2f} dB")
    return gaps

def main():
    cpp_csv = os.path.join(RESULTS_DIR, "cpp_bler_mcs2-27.csv")
    sionna_csv = os.path.join(RESULTS_DIR, "sionna_bler_mcs2-27.csv")
    
    print("Loading C++ results...")
    cpp_data = load_csv(cpp_csv)
    print(f"  Loaded {len(cpp_data)} MCS values")
    
    print("Loading Sionna results...")
    sionna_data = load_csv(sionna_csv)
    print(f"  Loaded {len(sionna_data)} MCS values")
    
    if not cpp_data and not sionna_data:
        print("ERROR: No data found! Run simulations first.")
        return
    
    out_overview = os.path.join(OUT_DIR, "bler_comparison_overview.png")
    plot_all_mcs(cpp_data, sionna_data, out_overview)
    
    out_per_mcs = os.path.join(OUT_DIR, "per_mcs")
    plot_per_mcs(cpp_data, sionna_data, out_per_mcs)
    
    compute_gap_summary(cpp_data, sionna_data)

if __name__ == "__main__":
    main()
