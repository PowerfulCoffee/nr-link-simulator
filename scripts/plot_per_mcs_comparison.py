import os, csv, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results')
OUTPUT_DIR = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUTPUT_DIR, exist_ok=True)

MCS_TABLE = [
    (2, 120),   # 0
    (2, 193),   # 1
    (2, 308),   # 2
    (2, 449),   # 3
    (2, 602),   # 4
    (4, 378),   # 5
    (4, 434),   # 6
    (4, 490),   # 7
    (4, 553),   # 8
    (4, 616),   # 9
    (4, 658),   # 10
    (6, 466),   # 11
    (6, 517),   # 12
    (6, 567),   # 13
    (6, 616),   # 14
    (6, 666),   # 15
    (6, 719),   # 16
    (6, 772),   # 17
    (6, 822),   # 18
    (6, 873),   # 19
    (8, 682.5), # 20
    (8, 711),   # 21
    (8, 754),   # 22
    (8, 797),   # 23
    (8, 841),   # 24
    (8, 885),   # 25
    (8, 916.5), # 26
    (8, 948),   # 27
]
TBS_MAP = {
    0:104, 1:176, 2:288, 3:408, 4:552,
    5:704, 6:808, 7:888, 8:1032, 9:1128, 10:1224,
    11:1288, 12:1416, 13:1608, 14:1736, 15:1864,
    16:2024, 17:2152, 18:2280, 19:2408,
    20:2472, 21:2600, 22:2792, 23:2976, 24:3104, 25:3240, 26:3368, 27:3496
}
MOD_NAMES = {2: 'QPSK', 4: '16QAM', 6: '64QAM', 8: '256QAM'}
MOD_COLORS = {2: '#1f77b4', 4: '#ff7f0e', 6: '#2ca02c', 8: '#d62728'}

def load_csv(path):
    data = {}
    if not os.path.exists(path):
        return data
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            mcs = int(row['mcs'])
            snr = float(row['snr_db'])
            bler = float(row['bler'])
            if mcs not in data:
                data[mcs] = []
            data[mcs].append((snr, bler))
    for mcs in data:
        data[mcs].sort()
    return data

def bler10(sbl):
    sbl = sorted(sbl)
    for i in range(len(sbl) - 1):
        s1, b1 = sbl[i]
        s2, b2 = sbl[i + 1]
        if b1 >= 0.1 and b2 <= 0.1:
            return s1 + (0.1 - b1) * (s2 - s1) / (b2 - b1)
    return None

cpp = load_csv(os.path.join(RESULTS_DIR, 'cpp_all_mcs_parsed.csv'))
sionna_fast = load_csv(os.path.join(RESULTS_DIR, 'sionna_all_mcs_awgn_fast.csv'))
refine_path = os.path.join(RESULTS_DIR, 'sionna_refined_key_mcs.csv')
sionna_refine = load_csv(refine_path) if os.path.exists(refine_path) else {}

sionna = dict(sionna_fast)
for mcs, pts in sionna_refine.items():
    rs = set(round(s, 1) for s, b in pts)
    existing = [(s, b) for s, b in sionna.get(mcs, []) if round(s, 1) not in rs]
    sionna[mcs] = sorted(existing + pts)

comparable_mcs = sorted(set(cpp.keys()) & set(sionna.keys()))
print(f"Comparable MCS (both C++ and Sionna): {comparable_mcs}")

n = len(comparable_mcs)
ncols = 4
nrows = (n + ncols - 1) // ncols

fig = plt.figure(figsize=(5 * ncols, 4.2 * nrows))
fig.suptitle('NR PDSCH BLER vs Es/N0: C++ nr-link-simulator vs Sionna (AWGN SISO, 3PRB)\n'
             'Solid= C++ (Offset Min-Sum LLR=Max-Log), Dashed= Sionna (BP LLR=log-MAP)',
             fontsize=14, fontweight='bold', y=1.01)

gs = GridSpec(nrows, ncols, figure=fig, hspace=0.45, wspace=0.35)

for idx, mcs in enumerate(comparable_mcs):
    ax = fig.add_subplot(gs[idx])
    qm, r1024 = MCS_TABLE[mcs]
    R = r1024 / 1024.0
    tbs = TBS_MAP[mcs]
    color = MOD_COLORS[qm]

    cpp_data = cpp[mcs]
    s_s_cpp = [x[0] for x in cpp_data]
    b_s_cpp = [x[1] for x in cpp_data]
    ax.semilogy(s_s_cpp, b_s_cpp, '-', color=color, linewidth=2.0, label='C++', zorder=3)
    ax.plot(s_s_cpp, b_s_cpp, 'o', color=color, markersize=4, zorder=4)

    s_data = sionna[mcs]
    s_s_s = [x[0] for x in s_data]
    b_s_s = [x[1] for x in s_data]
    ax.semilogy(s_s_s, b_s_s, '--', color=color, linewidth=2.0, label='Sionna', zorder=3)
    ax.plot(s_s_s, b_s_s, 's', color=color, markersize=4, markerfacecolor='white', markeredgewidth=1.5, zorder=4)

    ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.6, linewidth=1.0)

    b10_c = bler10(cpp_data)
    b10_s = bler10(s_data)
    gap = b10_c - b10_s if (b10_c is not None and b10_s is not None) else None

    title = f'MCS {mcs}  {MOD_NAMES[qm]}  R={R:.3f}  TBS={tbs}'
    if gap is not None:
        title += f'\nBLER=0.1 gap: {gap:+.2f} dB'
    ax.set_title(title, fontsize=9, fontweight='bold')

    ax.set_xlabel('Es/N0 (dB)', fontsize=8)
    ax.set_ylabel('BLER', fontsize=8)
    ax.set_ylim(5e-4, 1.2)
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=8, loc='lower left')

    if b10_c is not None:
        ax.axvline(x=b10_c, color=color, linestyle='-', alpha=0.2, linewidth=0.8)
    if b10_s is not None:
        ax.axvline(x=b10_s, color=color, linestyle='--', alpha=0.2, linewidth=0.8)

for idx in range(n, nrows * ncols):
    ax = fig.add_subplot(gs[idx])
    ax.set_visible(False)

output_path = os.path.join(OUTPUT_DIR, 'bler_per_mcs_comparison.png')
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved per-MCS comparison to: {output_path}")

gaps = []
for mcs in comparable_mcs:
    qm, r1024 = MCS_TABLE[mcs]
    R = r1024 / 1024.0
    bc = bler10(cpp[mcs])
    bs = bler10(sionna[mcs])
    if bc is not None and bs is not None:
        g = bc - bs
        gaps.append((mcs, qm, R, bc, bs, g))
        print(f"  MCS{mcs:2d} ({MOD_NAMES[qm]:6s} R={R:.3f}): C++={bc:6.2f}dB  Sionna={bs:6.2f}dB  gap={g:+.2f}dB")

if gaps:
    gap_vals = [g[5] for g in gaps]
    print(f"\nMean gap: {np.mean(gap_vals):+.2f} dB, Std: {np.std(gap_vals):.2f} dB")

summary_fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
b10c_d = {}; b10s_d = {}
for mcs in comparable_mcs:
    v = bler10(cpp[mcs])
    if v is not None: b10c_d[mcs] = v
    v = bler10(sionna[mcs])
    if v is not None: b10s_d[mcs] = v

for qm in [2, 4, 6, 8]:
    mc = sorted(b10c_d.keys())
    idx = [i for i, m in enumerate(mc) if MCS_TABLE[m][0] == qm]
    if idx:
        ax1.scatter([MCS_TABLE[mc[i]][1]/1024 for i in idx], [b10c_d[mc[i]] for i in idx],
                    c=MOD_COLORS[qm], label=f'{MOD_NAMES[qm]} C++', marker='o', s=70, zorder=5, edgecolors='k', linewidth=0.5)
    ms = sorted(b10s_d.keys())
    idx = [i for i, m in enumerate(ms) if MCS_TABLE[m][0] == qm]
    if idx:
        ax1.scatter([MCS_TABLE[ms[i]][1]/1024 for i in idx], [b10s_d[ms[i]] for i in idx],
                    c=MOD_COLORS[qm], label=f'{MOD_NAMES[qm]} Sionna', marker='s', s=55, zorder=4, alpha=0.6)
ax1.set_xlabel('Code Rate R', fontsize=11)
ax1.set_ylabel('Es/N0 @ BLER=0.1 (dB)', fontsize=11)
ax1.set_title('Required SNR for BLER=0.1 vs Code Rate', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=8, ncol=2, loc='upper left')

common = sorted(set(b10c_d.keys()) & set(b10s_d.keys()))
diffs = [b10c_d[m] - b10s_d[m] for m in common]
for m in common:
    qm = MCS_TABLE[m][0]
    ax2.bar(m, b10c_d[m] - b10s_d[m], color=MOD_COLORS[qm], alpha=0.8, width=0.7, edgecolor='k', linewidth=0.3)
md = np.mean(diffs) if diffs else 0
ax2.axhline(md, color='black', ls='--', lw=1.5, label=f'Mean = {md:+.2f} dB')
ax2.axhspan(-0.3, 0.3, alpha=0.1, color='green')
ax2.set_title('C++ - Sionna SNR Gap @ BLER=0.1', fontsize=12, fontweight='bold')
ax2.set_xlabel('MCS Index', fontsize=11)
ax2.set_ylabel('SNR Gap (dB)', fontsize=11)
ax2.grid(True, alpha=0.3, axis='y')
ax2.legend(fontsize=10)
ax2.set_xticks(range(0, 28, 2))
plt.tight_layout()
summary_path = os.path.join(OUTPUT_DIR, 'bler_summary_gap.png')
plt.savefig(summary_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved summary gap plot to: {summary_path}")
