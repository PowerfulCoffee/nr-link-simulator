import csv
import matplotlib.pyplot as plt
import numpy as np
import os

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
output_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
os.makedirs(output_dir, exist_ok=True)

MCS_CONFIGS = [
    {"mcs": 4,  "mod": "QPSK",   "color": "#1f77b4", "mcs_label": "MCS4 (QPSK, R=0.588)",  "x_max": 20},
    {"mcs": 8,  "mod": "16QAM",  "color": "#d62728", "mcs_label": "MCS8 (16QAM, R=0.540)", "x_max": 30},
    {"mcs": 15, "mod": "64QAM",  "color": "#2ca02c", "mcs_label": "MCS15 (64QAM, R=0.650)","x_max": 38},
    {"mcs": 27, "mod": "256QAM", "color": "#9467bd", "mcs_label": "MCS27 (256QAM, R=0.926)","x_max": 46},
]

CURVE_STYLES = [
    {"label": "C++ Ideal CE",   "ls": "-",  "marker": "o", "mf": "full",  "lw": 2.0, "fname": "cpp_ideal_mcs{}.csv"},
    {"label": "C++ LS CE",      "ls": "-",  "marker": "^", "mf": "full",  "lw": 2.0, "fname": "cpp_ls_mcs{}.csv"},
    {"label": "Sionna Ideal CE","ls": "--", "marker": "s", "mf": "none",  "lw": 1.5, "fname": "sionna_ideal_mcs{}.csv"},
    {"label": "Sionna LS CE",   "ls": "--", "marker": "v", "mf": "none",  "lw": 1.5, "fname": "sionna_ls_mcs{}.csv"},
]

def read_csv(filepath, snr_min=-5, snr_max=50):
    snr_list = []
    bler_list = []
    if not os.path.exists(filepath):
        return snr_list, bler_list
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                s = float(row['SINR_dB'])
                b = float(row['BLER'])
                if snr_min <= s <= snr_max:
                    snr_list.append(s)
                    if b <= 0:
                        b = 5e-4
                    bler_list.append(b)
            except:
                pass
    # Sort by SNR
    if snr_list:
        pairs = sorted(zip(snr_list, bler_list))
        snr_list = [p[0] for p in pairs]
        bler_list = [p[1] for p in pairs]
    return snr_list, bler_list

plt.rcParams.update({'font.size': 10})
fig, axes = plt.subplots(2, 2, figsize=(14, 11))
axes = axes.flatten()

for plt_idx, cfg in enumerate(MCS_CONFIGS):
    ax = axes[plt_idx]
    mcs = cfg['mcs']
    color = cfg['color']
    x_max = cfg['x_max']
    
    has_label = plt_idx == 0
    
    for style_idx, style in enumerate(CURVE_STYLES):
        filepath = os.path.join(results_dir, style["fname"].format(mcs))
        snr, bler = read_csv(filepath, snr_max=x_max)
        if snr:
            label = style["label"] if has_label else None
            ax.semilogy(snr, bler, color=color, linestyle=style["ls"], 
                        marker=style["marker"], linewidth=style["lw"], markersize=5,
                        markerfacecolor=color if style["mf"] == "full" else "none",
                        markeredgecolor=color, markeredgewidth=1.5,
                        label=label)
    
    ax.set_xlabel('SNR (dB)', fontsize=11)
    ax.set_ylabel('BLER', fontsize=11)
    ax.set_title(f'{cfg["mcs_label"]}', fontsize=12, fontweight='bold')
    ax.grid(True, which='both', linestyle='--', alpha=0.6)
    ax.set_ylim(5e-4, 2.0)
    ax.set_xlim(-5, x_max)
    ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5, linewidth=1)
    ax.set_xticks(range(-5, x_max + 1, 5))

fig.suptitle('C++ vs Sionna BLER: 4x4 MIMO TDL-A (DS=100ns, $f_d$=10Hz)\n30 PRB, Rank 4, DMRS Type1 Add-Pos=1 (sym 2&7), Ideal vs LS CE', fontsize=13, y=1.01)

handles = []
labels = []
for style in CURVE_STYLES:
    h = plt.Line2D([0], [0], color='black', linestyle=style["ls"], 
                   marker=style["marker"], linewidth=style["lw"], markersize=6,
                   markerfacecolor='black' if style["mf"] == "full" else "none",
                   markeredgecolor='black')
    handles.append(h)
    labels.append(style["label"])
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.01), ncol=4, fontsize=11)

plt.tight_layout()

subplot_file = os.path.join(output_dir, 'bler_comparison_subplots.png')
plt.savefig(subplot_file, dpi=150, bbox_inches='tight')
print(f"Figure saved to {subplot_file}")
plt.close()

# Print summary
print("\n" + "="*70)
print("SUMMARY: BLER vs SNR for all MCS")
print("="*70)
for cfg in MCS_CONFIGS:
    mcs = cfg['mcs']
    print(f"\n{cfg['mcs_label']}:")
    for style in CURVE_STYLES:
        filepath = os.path.join(results_dir, style["fname"].format(mcs))
        snr, bler = read_csv(filepath, snr_max=cfg['x_max'])
        if snr:
            bler_arr = np.array(bler)
            snr_arr = np.array(snr)
            below1 = np.where(bler_arr <= 0.1)[0]
            if len(below1) > 0:
                snr1 = f"{snr_arr[below1[0]]:.0f}dB"
            else:
                min_bler = np.min(bler_arr)
                min_snr = snr_arr[np.argmin(bler_arr)]
                snr1 = f">={cfg['x_max']}dB (min BLER={min_bler:.3f} @ {min_snr:.0f}dB)"
            
            # Also find BLER=1.0 to 0.1 transition
            bler1_idx = np.where(bler_arr < 1.0)[0]
            start = f"{snr_arr[bler1_idx[0]]:.0f}dB" if len(bler1_idx) > 0 else "N/A"
            
            print(f"  {style['label']:20s}: waterfall starts ~{start}, BLER<0.1 @ {snr1}")
        else:
            print(f"  {style['label']:20s}: NO DATA")
