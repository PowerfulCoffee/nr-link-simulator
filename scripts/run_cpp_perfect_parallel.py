#!/usr/bin/env python3
"""Parallel batch BLER simulation for C++ Perfect CE, MCS 3-27.
Incremental CSV writing, 2 workers to avoid CPU thrashing."""

import subprocess
import os
import csv
import sys
from multiprocessing import Pool

CPP_BIN = "/workspace/nr-link-simulator/build/examples/pdsch_bler_simulation"
OUT_DIR = "/workspace/nr-link-simulator/results/bler_batch"
OUT_CSV = os.path.join(OUT_DIR, "cpp_bler_perfect_mcs3-27.csv")
MAX_BLOCKS = 100
TARGET_ERRORS = 20
TIMEOUT = 1800

MCS_TABLE_1 = [
    (2, 120), (2, 157), (2, 193), (2, 251), (2, 308), (2, 379), (2, 449), (2, 526), (2, 602), (2, 679),
    (4, 340), (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658), (6, 438), (6, 466), (6, 517),
    (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873), (6, 910), (6, 948)
]

def get_snr_range_perfect(mcs):
    qm, R_x1024 = MCS_TABLE_1[mcs]
    R = R_x1024 / 1024.0
    if qm == 2:
        if R < 0.20: return (-5, 1)
        elif R < 0.30: return (-5, 2)
        elif R < 0.40: return (-4, 4)
        elif R < 0.50: return (-3, 6)
        else: return (-1, 9)
    elif qm == 4:
        if R < 0.40: return (1, 10)
        elif R < 0.50: return (3, 12)
        elif R < 0.60: return (5, 14)
        else: return (7, 16)
    elif qm == 6:
        if R < 0.50: return (6, 15)
        elif R < 0.60: return (8, 17)
        elif R < 0.70: return (10, 19)
        elif R < 0.80: return (12, 21)
        else: return (14, 24)
    elif qm == 8:
        return (18, 27)
    return (0, 30)

def run_one_mcs(mcs):
    qm, R_x1024 = MCS_TABLE_1[mcs]
    R = R_x1024 / 1024.0
    snr_start, snr_end = get_snr_range_perfect(mcs)
    
    cmd = [CPP_BIN, str(mcs), str(snr_start), str(snr_end), "1.0", "1", str(MAX_BLOCKS), str(TARGET_ERRORS)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        output = result.stdout
    except subprocess.TimeoutExpired:
        print(f"[MCS {mcs:2d}] TIMEOUT after {TIMEOUT}s", flush=True)
        return (mcs, [])
    
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
                    print(f"[MCS {mcs:2d}] SNR {snr:5.1f}dB: BLER={bler:.4f} ({n_err}/{n_blk})", flush=True)
                except ValueError:
                    pass
    print(f"[MCS {mcs:2d}] Done.", flush=True)
    return (mcs, results)

def main():
    mcs_start = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    mcs_end = int(sys.argv[2]) if len(sys.argv) > 2 else 27
    n_workers = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    
    done_mcs = set()
    if os.path.exists(OUT_CSV):
        with open(OUT_CSV, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                done_mcs.add(int(row['mcs']))
    
    mcs_list = [m for m in range(mcs_start, mcs_end + 1) if m not in done_mcs]
    print(f"Running MCS {mcs_list} with {n_workers} workers (max_blocks={MAX_BLOCKS}, target_err={TARGET_ERRORS})", flush=True)
    
    os.makedirs(OUT_DIR, exist_ok=True)
    write_header = not os.path.exists(OUT_CSV)
    
    csv_lock = None  # We write incrementally in main process, so use imap_unordered
    
    with open(OUT_CSV, 'a', newline='') as f:
        w = csv.writer(f)
        if write_header:
            w.writerow(["mcs", "qm", "rate", "snr_db", "bler", "n_blocks", "n_errors"])
        
        with Pool(n_workers) as pool:
            for (mcs, results) in pool.imap_unordered(run_one_mcs, mcs_list):
                qm, R_x1024 = MCS_TABLE_1[mcs]
                R = R_x1024 / 1024.0
                for (snr, bler, n_blk, n_err) in results:
                    w.writerow([mcs, qm, f"{R:.4f}", f"{snr:.1f}", f"{bler:.6f}", n_blk, n_err])
                f.flush()
    
    print(f"\nAll results saved to {OUT_CSV}")

if __name__ == "__main__":
    main()
