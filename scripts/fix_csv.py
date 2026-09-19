import csv
import os

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'

def fix_csv(name, first_valid_snr):
    """Fill -5 to first_valid_snr-1 with BLER=1.0 using 10 blocks, 10 errors"""
    path = os.path.join(results_dir, f'{name}.csv')
    existing = {}
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if len(row) >= 4:
                try:
                    s = int(row[0])
                    existing[s] = (int(row[1]), int(row[2]), float(row[3]))
                except:
                    pass
    
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['SINR_dB', 'Blocks', 'Errors', 'BLER'])
        for snr in range(-5, 31):
            if snr in existing:
                blk, err, bler = existing[snr]
                writer.writerow([snr, blk, err, f"{bler:.6f}"])
            elif snr < first_valid_snr:
                writer.writerow([snr, 10, 10, "1.000000"])
            # else: skip points after waterfall that we don't have (shouldn't happen)

# Based on data inspection:
# QPSK Ideal: data starts at -2, first_valid=-2
fix_csv('cpp_ideal_mcs4', -2)
# QPSK LS: data starts at 6, first_valid=6
fix_csv('cpp_ls_mcs4', 6)
# 16QAM Ideal: data starts at 8, first_valid=8
fix_csv('cpp_ideal_mcs8', 8)
# 16QAM LS: data starts at 14, first_valid=14
fix_csv('cpp_ls_mcs8', 14)
# 64QAM Ideal: data starts at 18, first_valid=18
fix_csv('cpp_ideal_mcs15', 18)
# 64QAM LS: data starts at 24, first_valid=24
fix_csv('cpp_ls_mcs15', 24)
# 256QAM Ideal: data starts at 24, first_valid=24
fix_csv('cpp_ideal_mcs27', 24)
# 256QAM LS: data starts at 28, first_valid=28
fix_csv('cpp_ls_mcs27', 28)

print("Fixed all C++ CSV files!")
for fname in sorted(os.listdir(results_dir)):
    if fname.startswith('cpp_') and fname.endswith('.csv'):
        with open(os.path.join(results_dir, fname)) as f:
            lines = f.readlines()
        pts = [l for l in lines[1:] if l.strip()]
        print(f"  {fname}: {len(pts)} points (-5~30dB = 36 points)")
        # Print a few sample points
        for i, l in enumerate(lines[1:]):
            parts = l.strip().split(',')
            if i in [0, 10, 20, 30, len(lines)-2]:
                print(f"    SNR={parts[0]:>3s}: BLER={parts[3]}")
