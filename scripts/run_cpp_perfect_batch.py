#!/usr/bin/env python3
"""Batch BLER simulation for C++ nr-link-simulator, MCS 3-27 AWGN Ideal/Perfect CE.
Uses tighter SNR ranges to reduce simulation time."""

import subprocess
import os
import csv
import sys

CPP_BIN = "/workspace/nr-link-simulator/build/examples/pdsch_bler_simulation"
OUT_DIR = "/workspace/nr-link-simulator/results/bler_batch"
os.makedirs(OUT_DIR, exist_ok=True)

# MCS Table 1 (PDSCH): (qm, R*1024)
MCS_TABLE_1 = [
    (2, 120), (2, 157), (2, 193), (2, 251), (2, 308), (2, 379), (2, 449), (2, 526), (2, 602), (2, 679),
    (4, 340), (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658), (6, 438), (6, 466), (6, 517),
    (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873), (6, 910), (6, 948)
]

def get_snr_range_perfect(mcs):
    """SNR ranges for Perfect CE ~1-2dB left of LS CE."""
    qm, R_x1024 = MCS_TABLE_1[mcs]
    R = R_x1024 / 1024.0
    if qm == 2:
        if R < 0.20:
            return (-5, 2)
        elif R < 0.30:
            return (-5, 3)
        elif R < 0.40:
            return (-4, 5)
        elif R < 0.50:
            return (-3, 7)
        else:
            return (-1, 10)
    elif qm == 4:
        if R < 0.40:
            return (1, 11)
        elif R < 0.50:
            return (3, 13)
        elif R < 0.60:
            return (5, 15)
        else:
            return (7, 17)
    elif qm == 6:
        if R < 0.50:
            return (6, 16)
        elif R < 0.60:
            return (8, 18)
        elif R < 0.70:
            return (10, 20)
        elif R < 0.80:
            return (12, 22)
        else:
            return (14, 24)
    elif qm == 8:
        return (18, 28)
    return (0, 30)

def run_cpp_bler(mcs, snr_start, snr_end, step=1.0):
    cmd = [CPP_BIN, str(mcs), str(snr_start), str(snr_end), str(step), "1"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
        output = result.stdout
    except subprocess.TimeoutExpired:
        print(f"  MCS {mcs} timed out")
        return []
    
    results = []
    in_results = False
    for line in output.split('\n'):
        if 'Results Summary' in line:
            in_results = True
            continue
        if in_results and '---' in line:
            continue
        if in_results:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    snr = float(parts[0])
                    n_blk = int(parts[1])
                    n_err = int(parts[2])
                    bler = float(parts[3])
                    results.append((snr, bler, n_blk, n_err))
                except ValueError:
                    pass
    return results

def main():
    mcs_start = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    mcs_end = int(sys.argv[2]) if len(sys.argv) > 2 else 27
    
    out_csv = os.path.join(OUT_DIR, f"cpp_bler_perfect_mcs{mcs_start}-{mcs_end}.csv")
    # Check for existing partial results
    done_mcs = set()
    if os.path.exists(out_csv):
        with open(out_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                done_mcs.add(int(row['mcs']))
        print(f"Resuming from {out_csv}, already have MCSs: {sorted(done_mcs)}")
        mode = 'a'
    else:
        mode = 'w'
    
    with open(out_csv, mode, newline="") as f:
        w = csv.writer(f)
        if mode == 'w':
            w.writerow(["mcs", "qm", "rate", "snr_db", "bler", "n_blocks", "n_errors"])
        
        for mcs in range(mcs_start, mcs_end + 1):
            if mcs in done_mcs:
                print(f"=== MCS {mcs} already done, skipping ===", flush=True)
                continue
            qm, R_x1024 = MCS_TABLE_1[mcs]
            R = R_x1024 / 1024.0
            snr_start, snr_end = get_snr_range_perfect(mcs)
            print(f"=== C++ Perfect CE MCS {mcs:2d} (Qm={qm}, R={R:.3f}) SNR {snr_start}~{snr_end}dB ===", flush=True)
            
            results = run_cpp_bler(mcs, snr_start, snr_end, 1.0)
            for snr, bler, n_blk, n_err in results:
                print(f"  SNR {snr:5.1f}dB: BLER={bler:.4f} ({n_err}/{n_blk})", flush=True)
                w.writerow([mcs, qm, f"{R:.4f}", f"{snr:.1f}", f"{bler:.6f}", n_blk, n_err])
            f.flush()
            print()
    
    print(f"Perfect CE results saved to {out_csv}")

if __name__ == "__main__":
    main()
