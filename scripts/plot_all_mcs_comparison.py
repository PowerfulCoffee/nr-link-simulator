import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'output')
os.makedirs(output_dir, exist_ok=True)
cpp_csv = os.path.join(results_dir, 'cpp_all_mcs_parsed.csv')
sionna_csv = os.path.join(results_dir, 'sionna_all_mcs_awgn_fast.csv')

def load_csv(path):
    data = {}
    if not os.path.exists(path):
        return data
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mcs = int(row['mcs'])
            snr = float(row['snr_db'])
            bler = float(row['bler'])
            if mcs not in data:
                data[mcs] = []
            data[mcs].append((snr, bler))
    for mcs in data:
        data[mcs].sort()
    return data

def calc_bler10(snr_bler_list):
    if not snr_bler_list:
        return None
    sbl = sorted(snr_bler_list)
    for i in range(len(sbl)-1):
        s1, b1 = sbl[i]
        s2, b2 = sbl[i+1]
        if b1 >= 0.1 and b2 <= 0.1:
            return s1 + (0.1 - b1) * (s2 - s1) / (b2 - b1)
    for s, b in sbl:
        if b < 0.1:
            return s
    return None

MCS_TABLE = [
    (2, 120), (2, 193), (2, 308), (2, 449), (2, 602),
    (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658),
    (6, 466), (6, 517), (6, 567), (6, 616), (6, 666),
    (6, 719), (6, 772), (6, 822), (6, 873),
    (8, 682.5), (8, 711), (8, 754), (8, 797), (8, 841), (8, 885), (8, 916.5), (8, 948),
]
mod_names = {2:'QPSK', 4:'16QAM', 6:'64QAM', 8:'256QAM'}
mod_colors = {2: '#1f77b4', 4: '#ff7f0e', 6: '#2ca02c', 8: '#d62728'}

cpp_data = load_csv(cpp_csv)
sionna_data = load_csv(sionna_csv)

print(f"C++ data: {len(cpp_data)} MCS loaded")
print(f"Sionna data: {len(sionna_data)} MCS loaded")

fig, axes = plt.subplots(2, 2, figsize=(18, 14))
fig.suptitle('NR PDSCH BLER vs Es/N0 - C++ nr-link-simulator vs Sionna (AWGN SISO, 3PRB)', fontsize=14, fontweight='bold')

mod_groups = {
    'QPSK (MCS 0-4)': list(range(0, 5)),
    '16QAM (MCS 5-10)': list(range(5, 11)),
    '64QAM (MCS 11-19)': list(range(11, 20)),
    '256QAM (MCS 20-27)': list(range(20, 28)),
}

for ax_idx, (title, mcs_list) in enumerate(mod_groups.items()):
    ax = axes[ax_idx // 2][ax_idx % 2]
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('Es/N0 (dB)', fontsize=10)
    ax.set_ylabel('BLER', fontsize=10)
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 1.0)
    ax.grid(True, alpha=0.3, which='both')
    ax.axhline(y=0.1, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    for mcs in mcs_list:
        qm, r1024 = MCS_TABLE[mcs]
        R = r1024 / 1024.0
        color = mod_colors[qm]
        alpha_base = 0.4 + 0.6 * (mcs - mcs_list[0]) / max(1, len(mcs_list)-1)

        if mcs in cpp_data:
            snrs = [x[0] for x in cpp_data[mcs]]
            blers = [x[1] for x in cpp_data[mcs]]
            ax.plot(snrs, blers, '-', color=color, alpha=alpha_base, linewidth=1.5,
                    label=f'MCS{mcs} C++ (R={R:.2f})')
            ax.plot(snrs, blers, 'o', color=color, alpha=alpha_base, markersize=3)

        if mcs in sionna_data:
            snrs = [x[0] for x in sionna_data[mcs]]
            blers = [x[1] for x in sionna_data[mcs]]
            ax.plot(snrs, blers, '--', color=color, alpha=alpha_base, linewidth=1.5,
                    label=f'MCS{mcs} Sionna (R={R:.2f})')
            ax.plot(snrs, blers, 's', color=color, alpha=alpha_base, markersize=3)

    ax.legend(fontsize=7, loc='lower left', ncol=2)

plt.tight_layout()
out_path1 = os.path.join(output_dir, 'bler_all_mcs_comparison.png')
plt.savefig(out_path1, dpi=150, bbox_inches='tight')
print(f"Saved BLER curves to {out_path1}")
plt.close()

bler10_cpp = {}
bler10_sionna = {}
for mcs in range(28):
    b10 = calc_bler10(cpp_data.get(mcs, []))
    if b10 is not None:
        bler10_cpp[mcs] = b10
    b10s = calc_bler10(sionna_data.get(mcs, []))
    if b10s is not None:
        bler10_sionna[mcs] = b10s

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

mcs_list_cpp = sorted(bler10_cpp.keys())
snrs_cpp = [bler10_cpp[m] for m in mcs_list_cpp]
rates = [MCS_TABLE[m][1]/1024.0 for m in mcs_list_cpp]
qms = [MCS_TABLE[m][0] for m in mcs_list_cpp]

for qm in [2, 4, 6, 8]:
    mask = [i for i, m in enumerate(mcs_list_cpp) if qms[i] == qm]
    if mask:
        ax1.scatter([rates[i] for i in mask], [snrs_cpp[i] for i in mask],
                   c=mod_colors[qm], label=f'{mod_names[qm]} C++', marker='o', s=60, alpha=0.8)

mcs_list_s = sorted(bler10_sionna.keys())
if mcs_list_s:
    snrs_s = [bler10_sionna[m] for m in mcs_list_s]
    rates_s = [MCS_TABLE[m][1]/1024.0 for m in mcs_list_s]
    qms_s = [MCS_TABLE[m][0] for m in mcs_list_s]
    for qm in [2, 4, 6, 8]:
        mask = [i for i, m in enumerate(mcs_list_s) if qms_s[i] == qm]
        if mask:
            ax1.scatter([rates_s[i] for i in mask], [snrs_s[i] for i in mask],
                       c=mod_colors[qm], label=f'{mod_names[qm]} Sionna', marker='s', s=60, alpha=0.5)

ax1.set_xlabel('Code Rate R', fontsize=11)
ax1.set_ylabel('Es/N0 @ BLER=0.1 (dB)', fontsize=11)
ax1.set_title('Required SNR for BLER=0.1 vs Code Rate', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=9)

diffs = []
mcs_diffs = []
for mcs in sorted(set(bler10_cpp.keys()) & set(bler10_sionna.keys())):
    d = bler10_cpp[mcs] - bler10_sionna[mcs]
    diffs.append(d)
    mcs_diffs.append(mcs)
    qm = MCS_TABLE[mcs][0]
    ax2.bar(mcs, d, color=mod_colors[qm], alpha=0.7, width=0.7)

if diffs:
    avg_diff = np.mean(diffs)
    ax2.axhline(y=avg_diff, color='black', linestyle='--', linewidth=1.5, label=f'Mean gap = {avg_diff:.2f} dB')
    ax2.set_title(f'C++ vs Sionna SNR Gap @ BLER=0.1 (positive = C++ worse)', fontsize=12, fontweight='bold')
    ax2.set_xlabel('MCS Index', fontsize=11)
    ax2.set_ylabel('SNR Gap (dB)', fontsize=11)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.legend(fontsize=10)
    ax2.set_xticks(range(0, 28, 2))
    print(f"\nC++ vs Sionna SNR gaps at BLER=0.1:")
    for m, d in zip(mcs_diffs, diffs):
        qm, r1024 = MCS_TABLE[m]
        print(f"  MCS{m:2d} ({mod_names[qm]:6s} R={r1024/1024:.3f}): gap = {d:+.2f} dB (C++ at {bler10_cpp[m]:.2f}, Sionna at {bler10_sionna[m]:.2f})")
    print(f"  Mean gap: {avg_diff:.2f} dB")
else:
    ax2.set_title('C++ vs Sionna SNR Gap (waiting for Sionna results...)', fontsize=12, fontweight='bold')
    ax2.text(0.5, 0.5, 'Sionna simulation running...\nGap plot will update when complete',
             transform=ax2.transAxes, ha='center', va='center', fontsize=14)

plt.tight_layout()
out_path2 = os.path.join(output_dir, 'bler_all_mcs_summary.png')
plt.savefig(out_path2, dpi=150, bbox_inches='tight')
print(f"Saved summary plot to {out_path2}")
plt.close()
