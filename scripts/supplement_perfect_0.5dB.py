#!/usr/bin/env python3
"""Targeted supplement of Perfect CE data: add 0.5dB points only around the
BLER cliff region where BLER transitions from ~1 to ~0.
Uses iterative refinement: add points, check if more are needed."""

import subprocess
import os
import csv
from multiprocessing import Pool

CPP_BIN = "/workspace/nr-link-simulator/build/examples/pdsch_bler_simulation"
OUT_CSV = "/workspace/nr-link-simulator/results/bler_batch/cpp_bler_perfect_mcs3-27.csv"
MAX_BLOCKS = 200
TARGET_ERRORS = 30
TIMEOUT = 300
N_WORKERS = 2

MCS_TABLE_1 = [
    (2, 120), (2, 157), (2, 193), (2, 251), (2, 308), (2, 379), (2, 449), (2, 526), (2, 602), (2, 679),
    (4, 340), (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658), (6, 438), (6, 466), (6, 517),
    (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873), (6, 910), (6, 948)
]

def load_data(csv_path):
    data = {}
    existing = set()
    if os.path.exists(csv_path):
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mcs = int(row['mcs'])
                snr = float(row['snr_db'])
                bler = float(row['bler'])
                if mcs not in data:
                    data[mcs] = {}
                data[mcs][snr] = bler
                existing.add((mcs, snr))
    return data, existing

def run_point(args):
    mcs, snr = args
    cmd = [CPP_BIN, str(mcs), f"{snr:.1f}", f"{snr:.1f}", "0.5", "1", str(MAX_BLOCKS), str(TARGET_ERRORS)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        output = result.stdout
    except subprocess.TimeoutExpired:
        print(f"[MCS {mcs:2d} SNR {snr:.1f}] TIMEOUT", flush=True)
        return None

    for line in output.split('\n'):
        parts = line.split()
        if len(parts) >= 4:
            try:
                s = float(parts[0])
                n_blk = int(parts[1])
                n_err = int(parts[2])
                bler = float(parts[3])
                if abs(s - snr) < 0.1:
                    print(f"[MCS {mcs:2d} SNR {snr:4.1f}dB] BLER={bler:.4f} ({n_err}/{n_blk})", flush=True)
                    qm, R_x1024 = MCS_TABLE_1[mcs]
                    R = R_x1024 / 1024.0
                    return [mcs, qm, f"{R:.4f}", f"{s:.1f}", f"{bler:.6f}", n_blk, n_err]
            except ValueError:
                pass
    return None

def find_cliff_gaps(mcs, data):
    """Find 0.5dB midpoints between consecutive points where BLER crosses the transition region."""
    if mcs not in data:
        return []
    snrs = sorted(data[mcs].keys())
    needed = []

    for i in range(len(snrs) - 1):
        s_lo = snrs[i]
        s_hi = snrs[i+1]
        gap = s_hi - s_lo
        if gap > 1.01 or gap < 0.49:
            continue
        b_lo = data[mcs][s_lo]
        b_hi = data[mcs][s_hi]

        crosses_transition = (b_lo > 0.1 and b_hi < 0.9) or (b_lo > 0.9 and b_hi < 0.1)
        if crosses_transition:
            mid = round(s_lo + 0.5, 1)
            needed.append(mid)

    return needed

def main():
    data, existing = load_data(OUT_CSV)
    print(f"Loaded {len(existing)} existing points", flush=True)

    to_run = []
    for mcs in range(3, 28):
        needed = find_cliff_gaps(mcs, data)
        for snr in needed:
            if (mcs, snr) not in existing:
                to_run.append((mcs, snr))

    print(f"Need to run {len(to_run)} transition midpoints:", flush=True)
    for (m, s) in to_run:
        print(f"  MCS {m:2d} @ {s:.1f}dB", flush=True)

    if not to_run:
        print("No points needed, done.")
        return

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    results = []
    with Pool(N_WORKERS) as pool:
        for res in pool.imap_unordered(run_point, to_run):
            if res is not None:
                results.append(res)

    with open(OUT_CSV, 'a', newline='') as f:
        w = csv.writer(f)
        for row in results:
            w.writerow(row)
        f.flush()

    print(f"\nWrote {len(results)} new points to {OUT_CSV}", flush=True)

if __name__ == "__main__":
    main()
