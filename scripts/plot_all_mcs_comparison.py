#!/usr/bin/env python3
"""
Plot MCS 0-27 BLER comparison: nr-link-simulator (C++) vs Sionna, AWGN SISO.
Generates two figures:
  1. A grid of subplots (4x7) showing all MCS BLER curves
  2. A summary plot showing BLER=0.1 SNR points for all MCS
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import csv

def load_csv(path):
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mcs = int(row['MCS'])
            if mcs not in data:
                data[mcs] = {'snr': [], 'bler': [], 'blocks': [], 'errors': [],
                            'mod': row['Modulation'], 'R': float(row['CodeRate']), 'qm': int(row['Qm'])}
            data[mcs]['snr'].append(float(row['SINR_dB']))
            data[mcs]['bler'].append(float(row['BLER']))
            data[mcs]['blocks'].append(int(row['Blocks']))
            data[mcs]['errors'].append(int(row['Errors']))
    for mcs in data:
        for k in ['snr','bler','blocks','errors']:
            data[mcs][k] = np.array(data[mcs][k])
    return data

def find_snr_at_bler(snrs, blers, target=0.1):
    if len(blers) < 2:
        return None
    snrs = np.asarray(snrs)
    blers = np.asarray(blers)
    for i in range(len(blers)-1):
        b1, b2 = blers[i], blers[i+1]
        s1, s2 = snrs[i], snrs[i+1]
        if b1 >= target and b2 <= target:
            if b1 == b2:
                return (s1+s2)/2
            t = (b1 - target) / (b1 - b2)
            return s1 + t*(s2-s1)
        if b1 == target:
            return s1
    return None

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    cpp_path = os.path.join(base_dir, 'build', 'examples', 'results', 'cpp_awgn_all_mcs.csv')
    sio_path = os.path.join(script_dir, 'results', 'sionna_awgn_all_mcs.csv')
    out_dir = os.path.join(script_dir, 'results')
    os.makedirs(out_dir, exist_ok=True)

    cpp_data = load_csv(cpp_path)
    sio_data = load_csv(sio_path) if os.path.exists(sio_path) else {}

    print(f"Loaded C++ data: MCS {min(cpp_data.keys())}-{max(cpp_data.keys())}")
    if sio_data:
        print(f"Loaded Sionna data: MCS {min(sio_data.keys())}-{max(sio_data.keys())}")

    mcs_list = sorted(cpp_data.keys())
    n_mcs = len(mcs_list)
    ncols = 4
    nrows = int(np.ceil(n_mcs / ncols))

    fig = plt.figure(figsize=(20, 24))
    gs = GridSpec(nrows, ncols, figure=fig, hspace=0.35, wspace=0.3)
    fig.suptitle('AWGN SISO BLER: nr-link-simulator (C++) vs Sionna, MCS 0-27\n'
                 '3 PRB, 15kHz SCS, LDPC 20 iter, perfect CSI, Min-Sum vs BP',
                 fontsize=16, fontweight='bold', y=0.98)

    colors = {'cpp': '#1f77b4', 'sionna': '#d62728'}
    markers = {'cpp': 'o', 'sionna': 's'}

    summary = []

    for idx, mcs in enumerate(mcs_list):
        ax = fig.add_subplot(gs[idx // ncols, idx % ncols])
        cd = cpp_data[mcs]
        mod = cd['mod']; R = cd['R']; qm = cd['qm']

        snr_cpp = cd['snr']; bler_cpp = cd['bler']
        err_cpp = cd['errors']; n_cpp = cd['blocks']
        ax.semilogy(snr_cpp, bler_cpp, color=colors['cpp'], marker=markers['cpp'],
                    markersize=4, linewidth=1.5, label='nr-link-sim')

        for i, (s, b, e, n) in enumerate(zip(snr_cpp, bler_cpp, err_cpp, n_cpp)):
            if e > 0 and n > 0:
                ci = 1.96 * np.sqrt(b * (1-b) / n)
                ax.errorbar(s, b, yerr=min(ci, b*0.8), color=colors['cpp'], capsize=2, elinewidth=0.8)

        if mcs in sio_data:
            sd = sio_data[mcs]
            snr_sio = sd['snr']; bler_sio = sd['bler']
            err_sio = sd['errors']; n_sio = sd['blocks']
            ax.semilogy(snr_sio, bler_sio, color=colors['sionna'], marker=markers['sionna'],
                        markersize=4, linewidth=1.5, label='Sionna')
            for s, b, e, n in zip(snr_sio, bler_sio, err_sio, n_sio):
                if e > 0 and n > 0:
                    ci = 1.96 * np.sqrt(b * (1-b) / n)
                    ax.errorbar(s, b, yerr=min(ci, b*0.8), color=colors['sionna'], capsize=2, elinewidth=0.8)

        snr10_cpp = find_snr_at_bler(snr_cpp, bler_cpp, 0.1)
        snr10_sio = find_snr_at_bler(sio_data[mcs]['snr'], sio_data[mcs]['bler'], 0.1) if mcs in sio_data else None
        gap = (snr10_cpp - snr10_sio) if (snr10_cpp is not None and snr10_sio is not None) else None
        summary.append((mcs, mod, R, qm, snr10_cpp, snr10_sio, gap))

        ax.axhline(y=0.1, color='gray', linestyle='--', linewidth=0.7, alpha=0.6)
        ax.set_xlabel('SNR (dB)', fontsize=8)
        ax.set_ylabel('BLER', fontsize=8)
        title = f'MCS{mcs} {mod} R={R:.3f}'
        if gap is not None:
            title += f' Δ={gap:+.1f}dB'
        ax.set_title(title, fontsize=9, fontweight='bold')
        ax.set_ylim(1e-4, 1.5)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc='upper right')
        ax.tick_params(labelsize=7)

    out_path = os.path.join(out_dir, 'bler_all_mcs_comparison.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved grid plot to {out_path}")

    fig2, axes = plt.subplots(1, 2, figsize=(16, 6))

    ax1 = axes[0]
    mcs_arr = np.array([s[0] for s in summary])
    snr10_cpp_arr = np.array([s[4] if s[4] is not None else np.nan for s in summary])
    snr10_sio_arr = np.array([s[5] if s[5] is not None else np.nan for s in summary])
    gap_arr = np.array([s[6] if s[6] is not None else np.nan for s in summary])

    ax1.plot(mcs_arr, snr10_cpp_arr, color=colors['cpp'], marker='o', linewidth=2, label='nr-link-sim')
    mask = ~np.isnan(snr10_sio_arr)
    ax1.plot(mcs_arr[mask], snr10_sio_arr[mask], color=colors['sionna'], marker='s', linewidth=2, label='Sionna')
    ax1.set_xlabel('MCS Index', fontsize=12)
    ax1.set_ylabel('SNR @ BLER=0.1 (dB)', fontsize=12)
    ax1.set_title('Required SNR for BLER=0.1 across all MCS', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(mcs_arr)

    ax2 = axes[1]
    ax2.bar(mcs_arr[mask], gap_arr[mask], color=['#2ca02c' if g < 1.0 else '#ff7f0e' if g < 2.0 else '#d62728' for g in gap_arr[mask]],
            edgecolor='black', linewidth=0.5)
    ax2.axhline(y=0, color='black', linewidth=0.5)
    ax2.axhline(y=1.0, color='green', linestyle='--', linewidth=1, alpha=0.5, label='1 dB gap')
    ax2.set_xlabel('MCS Index', fontsize=12)
    ax2.set_ylabel('SNR Gap (dB) [C++ - Sionna]', fontsize=12)
    ax2.set_title('Performance Gap (C++ relative to Sionna)', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xticks(mcs_arr)

    fig2.suptitle('MCS 0-27 Performance Summary (AWGN SISO)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    out_path2 = os.path.join(out_dir, 'bler_all_mcs_summary.png')
    plt.savefig(out_path2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved summary plot to {out_path2}")

    print("\n=== Summary: BLER=0.1 SNR (dB) ===")
    print(f"{'MCS':>4} {'Mod':>6} {'R':>6} {'Qm':>3} {'C++':>8} {'Sionna':>8} {'Gap':>8}")
    print("-" * 50)
    for mcs, mod, R, qm, s10c, s10s, gap in summary:
        sc = f'{s10c:+.2f}' if s10c is not None else '  N/A '
        ss = f'{s10s:+.2f}' if s10s is not None else '  N/A '
        sg = f'{gap:+.2f}' if gap is not None else '  N/A '
        print(f"{mcs:4d} {mod:>6} {R:6.3f} {qm:3d} {sc:>8} {ss:>8} {sg:>8}")

    csv_path = os.path.join(out_dir, 'bler_all_mcs_summary.csv')
    with open(csv_path, 'w') as f:
        f.write('MCS,Modulation,CodeRate,Qm,SNR_BLER10_Cpp,SNR_BLER10_Sionna,Gap_dB\n')
        for mcs, mod, R, qm, s10c, s10s, gap in summary:
            f.write(f'{mcs},{mod},{R:.6f},{qm},')
            f.write(f'{s10c:.4f},' if s10c is not None else ',')
            f.write(f'{s10s:.4f},' if s10s is not None else ',')
            f.write(f'{gap:.4f}\n' if gap is not None else '\n')
    print(f"\nSaved summary CSV to {csv_path}")

if __name__ == '__main__':
    main()
