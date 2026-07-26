#!/usr/bin/env python3
"""Plot MCS 3-27 SNR vs BLER comparison (C++ vs Sionna) in a grid layout,
   matching the style of output/bler_per_mcs_comparison.png."""

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

MOD_COLOR = {2: '#1f77b4', 4: '#ff7f0e', 6: '#2ca02c', 8: '#d62728'}

def mod_label(qm):
    return {2: 'QPSK', 4: '16QAM', 6: '64QAM', 8: '256QAM'}.get(qm, f'{qm}QAM')

def load_csv(path, has_tbs=False):
    data = {}
    tbs_map = {}
    if not os.path.exists(path):
        return data, tbs_map
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mcs = int(row['mcs'])
            snr = float(row['snr_db'])
            bler = float(row['bler'])
            n_blk = int(row['n_blocks'])
            n_err = int(row['n_errors'])
            if has_tbs and 'tbs' in row:
                tbs_map[mcs] = int(row['tbs'])
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
    return data, tbs_map

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

def get_snr_range(mcs):
    qm, R_x1024 = MCS_TABLE_1[mcs]
    R = R_x1024 / 1024.0
    if qm == 2:
        if R < 0.25:
            return (-4, 6)
        elif R < 0.45:
            return (-2, 10)
        else:
            return (1, 13)
    elif qm == 4:
        if R < 0.45:
            return (3, 15)
        elif R < 0.60:
            return (6, 18)
        else:
            return (8, 20)
    elif qm == 6:
        if R < 0.50:
            return (8, 20)
        elif R < 0.65:
            return (10, 23)
        elif R < 0.80:
            return (13, 26)
        else:
            return (16, 30)
    return (0, 30)

def main():
    out_dir = "/workspace/nr-link-simulator/results/bler_batch/plots"
    os.makedirs(out_dir, exist_ok=True)

    cpp_csv = "/workspace/nr-link-simulator/results/bler_batch/cpp_bler_mcs2-27.csv"
    sio_csv = "/workspace/nr-link-simulator/results/bler_batch/sionna_bler_mcs2-27.csv"

    cpp_data, _ = load_csv(cpp_csv, has_tbs=False)
    sio_data, sio_tbs = load_csv(sio_csv, has_tbs=True)

    mcs_list = list(range(3, 28))
    n_mcs = len(mcs_list)
    ncols = 4
    nrows = int(np.ceil(n_mcs / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4 * nrows))
    fig.suptitle(
        'NR PUSCH BLER vs SNR: C++ nr-link-simulator vs Sionna (AWGN SISO, 25PRB, LS-CE + DMRS-noise-est LMMSE)\n'
        'Solid= C++, Dashed= Sionna',
        fontsize=13, fontweight='bold', y=0.995
    )

    axes_flat = axes.flatten() if nrows > 1 or ncols > 1 else [axes]

    bler_min = 1e-4

    for idx, mcs in enumerate(mcs_list):
        ax = axes_flat[idx]
        qm, R_x1024 = MCS_TABLE_1[mcs]
        R = R_x1024 / 1024.0
        mod = mod_label(qm)
        color = MOD_COLOR.get(qm, '#333333')
        snr_lo, snr_hi = get_snr_range(mcs)

        has_cpp = mcs in cpp_data
        has_sio = mcs in sio_data

        if has_cpp:
            snr_c = cpp_data[mcs]['snr']
            bler_c = cpp_data[mcs]['bler'].copy()
            bler_c_plot = np.where(bler_c < 1e-10, bler_min, bler_c)
            ax.semilogy(snr_c, bler_c_plot, '-o', color=color, linewidth=1.5,
                        markersize=4, markerfacecolor='none', markeredgewidth=1.2,
                        markeredgecolor=color, label='C++', zorder=3)

        if has_sio:
            snr_s = sio_data[mcs]['snr']
            bler_s = sio_data[mcs]['bler'].copy()
            bler_s_plot = np.where(bler_s < 1e-10, bler_min, bler_s)
            ax.semilogy(snr_s, bler_s_plot, '--s', color=color, linewidth=1.5,
                        markersize=4, markerfacecolor='none', markeredgewidth=1.2,
                        markeredgecolor=color, label='Sionna', zorder=3)

        ax.axhline(y=0.1, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)
        ax.axvline(x=0, color='lightgray', linestyle='-', linewidth=0.4, alpha=0.3)

        ax.set_ylim([bler_min * 0.5, 2.0])
        ax.set_xlim([snr_lo - 1, snr_hi + 1])
        ax.grid(True, which='both', alpha=0.3, linewidth=0.5)
        ax.tick_params(labelsize=8)
        ax.set_xlabel('SNR (dB)', fontsize=8)
        ax.set_ylabel('BLER', fontsize=8)

        cpp_snr10 = find_snr_at_bler(cpp_data[mcs]['snr'], cpp_data[mcs]['bler'], 0.1) if has_cpp else None
        sio_snr10 = find_snr_at_bler(sio_data[mcs]['snr'], sio_data[mcs]['bler'], 0.1) if has_sio else None

        gap_str = ''
        if cpp_snr10 is not None and sio_snr10 is not None:
            gap = sio_snr10 - cpp_snr10
            gap_str = f' gap: {gap:+.2f} dB'
        elif cpp_snr10 is not None and not has_sio:
            gap_str = ' Sionna: N/A'
        elif sio_snr10 is not None and not has_cpp:
            gap_str = ' C++: N/A'

        tbs = sio_tbs.get(mcs, None)
        tbs_str = f' TBS={tbs}' if tbs is not None else ''

        title = f'MCS {mcs} {mod} R={R:.3f}{tbs_str}\nBLER=0.1{gap_str}'
        ax.set_title(title, fontsize=8, fontweight='bold')

        ax.legend(loc='lower left', fontsize=7, framealpha=0.8)

    for idx in range(n_mcs, nrows * ncols):
        axes_flat[idx].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.975])
    out_path = os.path.join(out_dir, 'bler_per_mcs_grid_mcs3-27.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved grid plot to {out_path}")

if __name__ == "__main__":
    main()
