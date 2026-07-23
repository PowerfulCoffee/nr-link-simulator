#!/usr/bin/env python3
"""
Run C++ and Sionna BLER simulations for TDL-A with Doppler
Compares: C++ LS-Doppler, C++ Ideal CE, Sionna LS (no Doppler), Sionna Perfect
"""
import subprocess
import os
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

CPP_BIN = "/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo"
SIONNA_SCRIPT = "/workspace/nr-link-simulator/scripts/sionna_tdl_doppler.py"
RESULTS_DIR = Path("/workspace/nr-link-simulator/scripts/results/tdl_doppler")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MCS_LIST = [3, 10, 17, 24]
MOD_NAMES = {3: "QPSK", 10: "16QAM", 17: "64QAM", 24: "256QAM"}

def get_sinr_range(mcs):
    if mcs == 3:
        return (-6, 10, 1.0)
    elif mcs == 10:
        return (2, 18, 1.0)
    elif mcs == 17:
        return (8, 24, 1.0)
    elif mcs == 24:
        return (14, 30, 1.0)
    return (0, 20, 1.0)

def run_cpp(mcs, perfect):
    sinr_start, sinr_end, sinr_step = get_sinr_range(mcs)
    max_blocks = 300
    target_errors = 50
    n_rb = 25
    n_rx = 1
    n_layers = 1
    ch_type = 1
    delay_spread = 100e-9
    max_doppler = 70.0
    dmrs_add = 1
    
    method = "cpp_perfect" if perfect else "cpp_ls_doppler"
    out_file = RESULTS_DIR / f"{method}_mcs{mcs}.csv"
    
    if out_file.exists():
        existing = out_file.read_text().strip().split('\n')
        if len(existing) > 2:
            last_line = existing[-1].split(',')
            if len(last_line) >= 4 and float(last_line[3]) < 0.02:
                print(f"[SKIP] {out_file.name} already complete")
                return str(out_file)
    
    cmd = [
        CPP_BIN,
        str(mcs), str(sinr_start), str(sinr_end), str(sinr_step),
        "1" if perfect else "0",
        str(max_blocks), str(target_errors),
        str(n_rb), str(n_rx), str(n_layers),
        str(ch_type), str(delay_spread), str(max_doppler), str(dmrs_add),
        str(out_file)
    ]
    
    label = f"C++ {'Perfect' if perfect else 'LS-Doppler'} MCS{mcs} ({MOD_NAMES[mcs]})"
    print(f"[RUN] {label}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            print(f"[FAIL] {label}: {result.stderr[-500:]}")
        else:
            print(f"[DONE] {label}")
    except subprocess.TimeoutExpired:
        print(f"[TIMEOUT] {label}")
    return str(out_file)

def run_sionna(mcs, perfect):
    method = "sionna_perfect" if perfect else "sionna_ls_nodoppler"
    out_file = RESULTS_DIR / f"{method}_mcs{mcs}.csv"
    
    if out_file.exists():
        existing = out_file.read_text().strip().split('\n')
        if len(existing) > 2:
            last_line = existing[-1].split(',')
            if len(last_line) >= 4 and float(last_line[3]) < 0.02:
                print(f"[SKIP] {out_file.name} already complete")
                return str(out_file)
    
    cmd = [
        "python3", SIONNA_SCRIPT,
        str(mcs), "1" if perfect else "0", str(out_file)
    ]
    
    label = f"Sionna {'Perfect' if perfect else 'LS-noDoppler'} MCS{mcs} ({MOD_NAMES[mcs]})"
    print(f"[RUN] {label}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, cwd="/workspace/nr-link-simulator/scripts")
        if result.returncode != 0:
            print(f"[FAIL] {label}: {result.stderr[-1000:]}")
        else:
            print(f"[DONE] {label}")
    except subprocess.TimeoutExpired:
        print(f"[TIMEOUT] {label}")
    return str(out_file)

def main():
    print("="*60)
    print("TDL-A Doppler BLER Simulation Suite")
    print(f"Channel: 25PRB SISO, TDL-A 100ns, fd=70Hz, DMRS add-pos=1")
    print("="*60)
    
    tasks = []
    for mcs in MCS_LIST:
        tasks.append(('cpp', mcs, False))
        tasks.append(('cpp', mcs, True))
        tasks.append(('sionna', mcs, False))
        tasks.append(('sionna', mcs, True))
    
    n_workers = 2
    print(f"\nRunning {len(tasks)} tasks with {n_workers} workers...\n")
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {}
        for task_type, mcs, perfect in tasks:
            if task_type == 'cpp':
                fut = executor.submit(run_cpp, mcs, perfect)
            else:
                fut = executor.submit(run_sionna, mcs, perfect)
            futures[fut] = (task_type, mcs, perfect)
        
        for fut in as_completed(futures):
            pass
    
    print("\n" + "="*60)
    print("All simulations complete!")
    print(f"Results in: {RESULTS_DIR}")
    print("="*60)

if __name__ == "__main__":
    main()
