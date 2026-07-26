import csv
import matplotlib.pyplot as plt
import numpy as np
import os

results_dir = '/workspace/nr-link-simulator/results/tdl_crosscheck_4x4'

MCS_CONFIGS = [
    {"mcs": 4,  "mod": "QPSK",   "color": "blue",    "marker_cpp": "o", "marker_sionna": "s"},
    {"mcs": 10, "mod": "16QAM",  "color": "red",     "marker_cpp": "o", "marker_sionna": "s"},
    {"mcs": 16, "mod": "64QAM",  "color": "green",   "marker_cpp": "o", "marker_sionna": "s"},
    {"mcs": 24, "mod": "256QAM", "color": "purple",  "marker_cpp": "o", "marker_sionna": "s"},
]

def read_csv(filepath):
    snr_list = []
    bler_list = []
    if not os.path.exists(filepath):
        return snr_list, bler_list
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            snr_list.append(float(row['SINR_dB']))
            bler_list.append(float(row['BLER']))
    return snr_list, bler_list

plt.figure(figsize=(10, 7))

for cfg in MCS_CONFIGS:
    mcs = cfg['mcs']
    mod = cfg['mod']
    color = cfg['color']
    
    cpp_file = os.path.join(results_dir, f'cpp_ideal_mcs{mcs}.csv')
    sionna_file = os.path.join(results_dir, f'sionna_ideal_mcs{mcs}.csv')
    
    cpp_snr, cpp_bler = read_csv(cpp_file)
    sionna_snr, sionna_bler = read_csv(sionna_file)
    
    if cpp_snr:
        plt.semilogy(cpp_snr, cpp_bler, color=color, marker=cfg['marker_cpp'], 
                     linestyle='-', linewidth=2, markersize=8, label=f'C++ Ideal CE {mod}')
    
    if sionna_snr:
        plt.semilogy(sionna_snr, sionna_bler, color=color, marker=cfg['marker_sionna'],
                     linestyle='--', linewidth=2, markersize=8, markerfacecolor='none',
                     label=f'Sionna Ideal CE {mod}')

plt.xlabel('SNR (dB)', fontsize=12)
plt.ylabel('BLER', fontsize=12)
plt.title('BLER vs SNR: 30 PRB, 4x4 Rank 4, TDL-A (DS=100ns), DMRS Type1 Add-Pos=1\nIdeal Channel Estimation', fontsize=12)
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(fontsize=10, bbox_to_anchor=(1.02, 1), loc='upper left')
plt.xlim(-6, 32)
plt.ylim(1e-3, 1.5)
plt.tight_layout()

output_file = os.path.join(results_dir, 'bler_comparison_ideal_ce.png')
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"Plot saved to {output_file}")

plt.figure(figsize=(12, 10))

plt_idx = 1
for cfg in MCS_CONFIGS:
    mcs = cfg['mcs']
    mod = cfg['mod']
    
    plt.subplot(2, 2, plt_idx)
    
    cpp_file = os.path.join(results_dir, f'cpp_ideal_mcs{mcs}.csv')
    sionna_file = os.path.join(results_dir, f'sionna_ideal_mcs{mcs}.csv')
    
    cpp_snr, cpp_bler = read_csv(cpp_file)
    sionna_snr, sionna_bler = read_csv(sionna_file)
    
    if cpp_snr:
        plt.semilogy(cpp_snr, cpp_bler, 'b-o', linewidth=2, markersize=7, label='C++ Ideal CE')
    if sionna_snr:
        plt.semilogy(sionna_snr, sionna_bler, 'r--s', linewidth=2, markersize=7, 
                     markerfacecolor='none', label='Sionna Ideal CE')
    
    plt.xlabel('SNR (dB)', fontsize=10)
    plt.ylabel('BLER', fontsize=10)
    plt.title(f'{mod} (MCS{mcs}, Table 2)', fontsize=11)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.legend(fontsize=9)
    plt.ylim(1e-3, 1.5)
    
    plt_idx += 1

plt.suptitle('C++ vs Sionna BLER Comparison: 4x4 MIMO TDL-A (DS=100ns)\n30 PRB, Rank 4, DMRS Type1 (sym 2&7), Ideal CE', fontsize=13)
plt.tight_layout()

subplot_file = os.path.join(results_dir, 'bler_comparison_subplots.png')
plt.savefig(subplot_file, dpi=150, bbox_inches='tight')
print(f"Subplot figure saved to {subplot_file}")

print("\n=== Summary of Results ===")
for cfg in MCS_CONFIGS:
    mcs = cfg['mcs']
    mod = cfg['mod']
    cpp_file = os.path.join(results_dir, f'cpp_ideal_mcs{mcs}.csv')
    sionna_file = os.path.join(results_dir, f'sionna_ideal_mcs{mcs}.csv')
    cpp_snr, cpp_bler = read_csv(cpp_file)
    sionna_snr, sionna_bler = read_csv(sionna_file)
    print(f"\n{mod} (MCS{mcs}):")
    if cpp_snr:
        print(f"  C++:    SNR range [{cpp_snr[0]}, {cpp_snr[-1]}] dB, points={len(cpp_snr)}")
    if sionna_snr:
        print(f"  Sionna: SNR range [{sionna_snr[0]}, {sionna_snr[-1]}] dB, points={len(sionna_snr)}")
