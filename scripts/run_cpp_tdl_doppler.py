#!/usr/bin/env python3
"""
Run C++ BLER simulations for TDL-A with Doppler (3PRB SISO)
"""
import subprocess
import os
from pathlib import Path

CPP_BIN = "/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo"
RESULTS_DIR = Path("/workspace/nr-link-simulator/scripts/results/tdl_doppler")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MCS_LIST = [
    (3,  -6, 8,  1.0, "QPSK"),
    (10,  2, 18, 1.0, "16QAM"),
    (17,  8, 24, 1.0, "64QAM"),
    (27, 16, 32, 1.0, "64QAM-high"),
]

def run_cpp(mcs, sinr_start, sinr_end, sinr_step, perfect, label):
    max_blocks = 500
    target_errors = 100
    n_rb = 3
    n_rx = 1
    n_layers = 1
    ch_type = 1
    delay_spread = 100e-9
    max_doppler = 70.0
    dmrs_add = 1
    
    method = "cpp_perfect" if perfect else "cpp_ls_doppler"
    out_file = RESULTS_DIR / f"{method}_mcs{mcs}.csv"
    
    if out_file.exists():
        lines = out_file.read_text().strip().split('\n')
        if len(lines) > 3:
            last = lines[-1].split(',')
            if len(last) >= 4 and float(last[3]) < 0.02:
                print(f"[SKIP] {label} already done")
                return
    
    cmd = [
        CPP_BIN,
        str(mcs), str(sinr_start), str(sinr_end), str(sinr_step),
        "1" if perfect else "0",
        str(max_blocks), str(target_errors),
        str(n_rb), str(n_rx), str(n_layers),
        str(ch_type), str(delay_spread), str(max_doppler), str(dmrs_add),
        str(out_file)
    ]
    
    print(f"[RUN] {label}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    print(f"[DONE] {label}")
    if result.returncode != 0:
        print(f"  stderr: {result.stderr[-300:]}")

def main():
    print("Running C++ TDL-A Doppler simulations (3PRB SISO)...")
    for mcs, s_start, s_end, s_step, mod in MCS_LIST:
        run_cpp(mcs, s_start, s_end, s_step, False, f"C++ LS-Doppler MCS{mcs} ({mod})")
        run_cpp(mcs, s_start, s_end, s_step, True, f"C++ Ideal CE MCS{mcs} ({mod})")
    print("C++ simulations complete!")

if __name__ == "__main__":
    main()
