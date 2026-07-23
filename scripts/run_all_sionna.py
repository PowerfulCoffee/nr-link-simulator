#!/usr/bin/env python3
"""Run all Sionna simulations (Perfect CE + LS) for 4 MCS levels"""
import subprocess
import os
from pathlib import Path

SIONNA_SCRIPT = "/workspace/nr-link-simulator/scripts/sionna_tdl_v3.py"
RESULTS_DIR = Path("/workspace/nr-link-simulator/scripts/results/tdl_doppler")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MCS_LIST = [3, 10, 17, 27]

def run_sionna(mcs, perfect):
    method = "sionna_perfect" if perfect else "sionna_ls_nodoppler"
    out_file = RESULTS_DIR / f"{method}_mcs{mcs}.csv"
    
    cmd = ["python3", "-u", SIONNA_SCRIPT, str(mcs), "1" if perfect else "0", str(out_file)]
    label = f"Sionna {'Perfect' if perfect else 'LS-noDoppler'} MCS{mcs}"
    print(f"[RUN] {label}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800, cwd="/workspace/nr-link-simulator/scripts")
    print(f"[DONE] {label}")
    if result.returncode != 0:
        print(f"  stderr: {result.stderr[-500:]}")
    return out_file

if __name__ == "__main__":
    import sys
    only_perfect = "--perfect" in sys.argv
    only_ls = "--ls" in sys.argv
    
    for mcs in MCS_LIST:
        if not only_ls:
            run_sionna(mcs, True)
        if not only_perfect:
            run_sionna(mcs, False)
    print("All Sionna simulations complete!")
