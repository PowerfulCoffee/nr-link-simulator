import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
import os
from scipy import stats

def load_csv(path):
    sinr = []
    bler = []
    nblocks = []
    nerrors = []
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sinr.append(float(row['SINR_dB']))
            bler.append(float(row['BLER']))
            nblocks.append(int(row['Blocks']))
            nerrors.append(int(row['Errors']))
    return np.array(sinr), np.array(bler), np.array(nblocks), np.array(nerrors)

def binomial_ci(k, n, confidence=0.95):
    if n == 0:
        return 0.0, 0.0
    p_hat = k / n
    if k == 0:
        lower = 0.0
        upper = 1 - (1 - confidence) ** (1 / n)
    elif k == n:
        lower = (1 - confidence) ** (1 / n)
        upper = 1.0
    else:
        alpha = 1 - confidence
        z = stats.norm.ppf(1 - alpha / 2)
        se = np.sqrt(p_hat * (1 - p_hat) / n)
        lower = max(0.0, p_hat - z * se)
        upper = min(1.0, p_hat + z * se)
    return p_hat - lower, upper - p_hat

script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, 'results')

cpp_tdla_path = os.path.join(results_dir, 'bler_mcs27_tdla_perfect_csi.csv')
sionna_tdla_path = os.path.join(results_dir, 'sionna_mcs27_tdla_perfect_csi.csv')

fig, ax = plt.subplots(1, 1, figsize=(11, 8))

if os.path.exists(cpp_tdla_path):
    sinr_cpp, bler_cpp, n_cpp, k_cpp = load_csv(cpp_tdla_path)
    err_low = np.zeros_like(bler_cpp)
    err_high = np.zeros_like(bler_cpp)
    for i in range(len(bler_cpp)):
        el, eh = binomial_ci(k_cpp[i], n_cpp[i])
        err_low[i] = el
        err_high[i] = eh
    ax.errorbar(sinr_cpp, bler_cpp, yerr=[err_low, err_high],
                fmt='bo-', linewidth=2, markersize=7, capsize=4,
                label='C++ nr-link-simulator (TDL-A, perfect CSI)',
                markerfacecolor='blue', markeredgecolor='blue', ecolor='blue', elinewidth=1)
    for i in range(len(bler_cpp)):
        ax.annotate(f'{k_cpp[i]}/{n_cpp[i]}',
                    (sinr_cpp[i], bler_cpp[i]),
                    textcoords="offset points", xytext=(8, 8),
                    fontsize=7, color='blue', alpha=0.7)

if os.path.exists(sionna_tdla_path):
    sinr_sio, bler_sio, n_sio, k_sio = load_csv(sionna_tdla_path)
    ax.semilogy(sinr_sio, bler_sio, 'rs--', linewidth=2, markersize=8,
                label='Sionna 2.x (TDL-A, perfect CSI)')
    for i in range(len(bler_sio)):
        ax.annotate(f'{k_sio[i]}/{n_sio[i]}',
                    (sinr_sio[i], bler_sio[i]),
                    textcoords="offset points", xytext=(8, -12),
                    fontsize=7, color='red', alpha=0.7)

awgn_path = os.path.join(results_dir, 'bler_mcs27_awgn.csv')
if os.path.exists(awgn_path):
    sinr_awgn, bler_awgn, n_awgn, k_awgn = load_csv(awgn_path)
    ax.semilogy(sinr_awgn, bler_awgn, 'g^-', linewidth=2, markersize=7,
                label='C++ AWGN SISO (reference)')

ax.set_xlabel('Es/N0 (dB)', fontsize=13)
ax.set_ylabel('BLER (log scale)', fontsize=13)
ax.set_title('NR PDSCH MCS27 BLER vs Es/N0\n'
             '(3 PRB, 15kHz SCS, 64QAM R=910/1024, TBS=2472, LDPC 20 iter, Offset Min-Sum)',
             fontsize=12)
ax.set_yscale('log')
ax.set_ylim([5e-4, 1.0])
ax.set_xlim([18, 42])
ax.grid(True, which='both', linestyle='--', alpha=0.7)
ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5)
ax.text(41.5, 0.11, 'BLER=0.1', ha='right', fontsize=9, color='gray')

ax.text(0.02, 0.02,
        'Bars: 95% binomial confidence interval\n'
        'Labels: errors/total blocks\n'
        'Note: At BLER ≈ 0.001-0.002 with N=1000 blocks,\n'
        'only 1-2 errors are observed; statistical\n'
        'uncertainty is large. The apparent uptick at\n'
        '38-40 dB is within noise.',
        transform=ax.transAxes, fontsize=8, verticalalignment='bottom',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

ax.legend(fontsize=10, loc='lower left')

plt.tight_layout()
out_path = os.path.join(script_dir, 'bler_mcs27_waterfall.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"Plot saved to {out_path}")

print("\n=== Summary ===")
if os.path.exists(cpp_tdla_path):
    print("C++ TDL-A results (errors/blocks = BLER, 95% CI):")
    for i, (s, b) in enumerate(zip(sinr_cpp, bler_cpp)):
        el, eh = binomial_ci(k_cpp[i], n_cpp[i])
        print(f"  Es/N0={s:5.1f} dB  {k_cpp[i]:3d}/{n_cpp[i]:5d}  BLER={b:.4f}  [{b-el:.4f}, {b+eh:.4f}]")
if os.path.exists(sionna_tdla_path):
    print("\nSionna TDL-A results:")
    for s, b in zip(sinr_sio, bler_sio):
        print(f"  Es/N0={s:5.1f} dB  BLER={b:.4f}")
