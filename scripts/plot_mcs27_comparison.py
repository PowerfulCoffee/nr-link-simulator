import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(script_dir)
output_dir = os.path.join(base_dir, 'output')
os.makedirs(output_dir, exist_ok=True)

cpp_csv = os.path.join(base_dir, 'results', 'bler_mcs27_awgn_fine.csv')
sionna_txt = os.path.join(base_dir, 'results', 'sionna_mcs27_awgn.txt')

cpp_snr = []
cpp_bler = []
with open(cpp_csv, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        cpp_snr.append(float(row['snr_db']))
        cpp_bler.append(float(row['bler']))

sionna_snr = []
sionna_bler = []
with open(sionna_txt, 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split(',')
        if len(parts) >= 4:
            sionna_snr.append(float(parts[0]))
            sionna_bler.append(float(parts[3]))

plt.figure(figsize=(10, 6))
plt.semilogy(cpp_snr, cpp_bler, 'bo-', linewidth=2, markersize=8, label='nr-link-simulator (C++, Min-Sum)')
plt.semilogy(sionna_snr, sionna_bler, 'rs-', linewidth=2, markersize=8, label='NVIDIA Sionna (Python, BP)')

plt.xlabel('SNR (dB)', fontsize=12)
plt.ylabel('BLER', fontsize=12)
plt.title('MCS 27 (64QAM, R≈0.889) BLER vs SNR - AWGN Channel, Perfect CSI\n24 PRB, 1x1 SISO, LDPC 30 iterations', fontsize=13)
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(fontsize=11)
plt.ylim(1e-4, 2)
plt.xlim(15, 22)

plt.axhline(y=0.1, color='g', linestyle=':', alpha=0.5, label='BLER=0.1')

plt.tight_layout()
out_png = os.path.join(output_dir, 'mcs27_bler_awgn_comparison.png')
plt.savefig(out_png, dpi=150)
print(f"Plot saved to: {out_png}")

print("\n=== MCS27 AWGN BLER Comparison ===")
print(f"{'SNR(dB)':<10} {'C++ BLER':<12} {'Sionna BLER':<12}")
print("-" * 35)
for i, snr in enumerate(cpp_snr):
    c_bler = cpp_bler[i]
    s_match = [s for s, b in zip(sionna_snr, sionna_bler) if abs(s - snr) < 0.3]
    s_bler_str = "---"
    if s_match:
        idx = sionna_snr.index(s_match[0])
        s_bler_str = f"{sionna_bler[idx]:.4f}"
    print(f"{snr:<10.1f} {c_bler:<12.4f} {s_bler_str:<12}")

print("\n=== Sionna Reference Points ===")
for snr, bler in zip(sionna_snr, sionna_bler):
    print(f"  {snr:.1f} dB: BLER = {bler:.4f}")

def find_bler_at(snr_list, bler_list, target=0.1):
    for i in range(len(snr_list)-1):
        if (bler_list[i] - target) * (bler_list[i+1] - target) <= 0:
            x1, x2 = snr_list[i], snr_list[i+1]
            y1, y2 = bler_list[i], bler_list[i+1]
            if y1 != y2:
                return x1 + (target - y1) * (x2 - x1) / (y2 - y1)
    return None

cpp_bler10 = find_bler_at(cpp_snr, cpp_bler, 0.1)
sionna_bler10 = find_bler_at(sionna_snr, sionna_bler, 0.1)
print(f"\nBLER=0.1 threshold:")
if cpp_bler10:
    print(f"  nr-link-simulator: {cpp_bler10:.2f} dB")
if sionna_bler10:
    print(f"  Sionna:           {sionna_bler10:.2f} dB")
if cpp_bler10 and sionna_bler10:
    print(f"  Gap:              {cpp_bler10 - sionna_bler10:.2f} dB")
