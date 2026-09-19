import csv
import matplotlib.pyplot as plt
import numpy as np
import os

results_dir = '/workspace/nr-link-simulator/results/ce_compare'

MCS_CONFIGS = [
    {"mcs": 4,  "mod": "QPSK",   "color": "#1f77b4", "mcs_label": "MCS4 (R=0.588)"},
    {"mcs": 8,  "mod": "16QAM",  "color": "#d62728", "mcs_label": "MCS8 (R=0.540)"},
    {"mcs": 15, "mod": "64QAM",  "color": "#2ca02c", "mcs_label": "MCS15 (R=0.650)"},
]

CURVE_STYLES = [
    {"label": "C++ Ideal CE",   "ls": "-",  "marker": "o", "mf": "full"},
    {"label": "C++ LS CE",      "ls": "-",  "marker": "^", "mf": "full"},
    {"label": "Sionna Ideal CE","ls": "--", "marker": "s", "mf": "none"},
    {"label": "Sionna LS CE",   "ls": "--", "marker": "v", "mf": "none"},
]

FILE_PATTERNS = [
    "cpp_ideal_mcs{}.csv",
    "cpp_ls_mcs{}.csv",
    "sionna_ideal_mcs{}.csv",
    "sionna_ls_mcs{}.csv",
]

def read_csv(filepath):
    snr_list = []
    bler_list = []
    if not os.path.exists(filepath):
        print(f"Warning: file not found: {filepath}")
        return snr_list, bler_list
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            snr_list.append(float(row['SINR_dB']))
            bler_list.append(float(row['BLER']))
    return snr_list, bler_list

plt.figure(figsize=(18, 6))

plt_idx = 1
for cfg in MCS_CONFIGS:
    mcs = cfg['mcs']
    mod = cfg['mod']
    color = cfg['color']
    
    ax = plt.subplot(1, 3, plt_idx)
    
    for style_idx, (style, fpat) in enumerate(zip(CURVE_STYLES, FILE_PATTERNS)):
        filepath = os.path.join(results_dir, fpat.format(mcs))
        snr, bler = read_csv(filepath)
        if snr:
            ax.semilogy(snr, bler, color=color, linestyle=style["ls"], 
                        marker=style["marker"], linewidth=2, markersize=7,
                        markerfacecolor=color if style["mf"] == "full" else "none",
                        markeredgecolor=color,
                        label=style["label"])
    
    ax.set_xlabel('SNR (dB)', fontsize=11)
    ax.set_ylabel('BLER', fontsize=11)
    ax.set_title(f'{mod} ({cfg["mcs_label"]})', fontsize=12)
    ax.grid(True, which='both', linestyle='--', alpha=0.7)
    ax.legend(fontsize=9)
    ax.set_ylim(1e-3, 1.5)
    
    if mod == "QPSK":
        ax.set_xlim(-3, 20)
    elif mod == "16QAM":
        ax.set_xlim(5, 30)
    else:
        ax.set_xlim(13, 40)
    
    plt_idx += 1

plt.suptitle('C++ vs Sionna BLER Comparison: 4x4 MIMO TDL-A (DS=15ns, fd=5Hz)\n30 PRB, Rank 4, DMRS Type1 (sym 2&7), Ideal vs LS Channel Estimation', fontsize=13)
plt.tight_layout()

output_dir = '/workspace/nr-link-simulator/results/tdl_crosscheck_4x4'
os.makedirs(output_dir, exist_ok=True)
subplot_file = os.path.join(output_dir, 'bler_comparison_subplots.png')
plt.savefig(subplot_file, dpi=150, bbox_inches='tight')
print(f"Subplot figure saved to {subplot_file}")

plt.figure(figsize=(12, 8))
for cfg in MCS_CONFIGS:
    mcs = cfg['mcs']
    mod = cfg['mod']
    color = cfg['color']
    for style_idx, (style, fpat) in enumerate(zip(CURVE_STYLES, FILE_PATTERNS)):
        filepath = os.path.join(results_dir, fpat.format(mcs))
        snr, bler = read_csv(filepath)
        if snr:
            if style_idx == 0:
                plt.semilogy(snr, bler, color=color, linestyle=style["ls"],
                          marker=style["marker"], linewidth=2, markersize=7,
                          markerfacecolor=color, label=f'{mod} {style["label"]}')
            elif style_idx == 3:
                plt.semilogy(snr, bler, color=color, linestyle=style["ls"],
                          marker=style["marker"], linewidth=2, markersize=7,
                          markerfacecolor="none", markeredgecolor=color, label=f'{mod} {style["label"]}')
            else:
                plt.semilogy(snr, bler, color=color, linestyle=style["ls"],
                          marker=style["marker"], linewidth=1.5, markersize=6,
                          markerfacecolor=color if style["mf"] == "full" else "none",
                          markeredgecolor=color, alpha=0.8)

plt.xlabel('SNR (dB)', fontsize=12)
plt.ylabel('BLER', fontsize=12)
plt.title('BLER vs SNR: 30 PRB, 4x4 Rank 4, TDL-A (DS=15ns, fd=5Hz)\nDMRS Type1 Add-Pos=1, Ideal vs LS Channel Estimation', fontsize=12)
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(fontsize=9, bbox_to_anchor=(1.02, 1), loc='upper left')
plt.xlim(-3, 40)
plt.ylim(1e-3, 1.5)
plt.tight_layout()
allinone_file = os.path.join(output_dir, 'bler_comparison_all.png')
plt.savefig(allinone_file, dpi=150, bbox_inches='tight')
print(f"All-in-one figure saved to {allinone_file}")

print("\n=== Summary of Results ===")
for cfg in MCS_CONFIGS:
    mcs = cfg['mcs']
    mod = cfg['mod']
    print(f"\n{mod} ({cfg['mcs_label']}):")
    for style, fpat in zip(CURVE_STYLES, FILE_PATTERNS):
        filepath = os.path.join(results_dir, fpat.format(mcs))
        snr, bler = read_csv(filepath)
        if snr:
            bler_arr = np.array(bler)
            snr_arr = np.array(snr)
            below1 = np.where(bler_arr <= 0.1)[0]
            below01 = np.where(bler_arr <= 0.01)[0]
            snr1 = f"{snr_arr[below1[0]]:.0f}dB" if len(below1) > 0 else "N/A"
            snr01 = f"{snr_arr[below01[0]]:.0f}dB" if len(below01) > 0 else "N/A"
            print(f"  {style['label']:20s}: SNR range [{snr[0]:.0f}, {snr[-1]:.0f}]dB, BLER<0.1 @ {snr1}, BLER<0.01 @ {snr01}")
