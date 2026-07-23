#!/usr/bin/env python3
"""Run remaining Sionna simulations with fast settings"""
import subprocess
import os
from pathlib import Path

SIONNA_SCRIPT = "/workspace/nr-link-simulator/scripts/sionna_tdl_v3.py"
RESULTS_DIR = Path("/workspace/nr-link-simulator/scripts/results/tdl_doppler")

tasks = [
    (10, False, "MCS10 LS"),
    (17, True,  "MCS17 Perfect"),
    (17, False, "MCS17 LS"),
    (27, True,  "MCS27 Perfect"),
    (27, False, "MCS27 LS"),
]

for mcs, perfect, label in tasks:
    method = "sionna_perfect" if perfect else "sionna_ls_nodoppler"
    out_file = RESULTS_DIR / f"{method}_mcs{mcs}.csv"
    if out_file.exists():
        lines = out_file.read_text().strip().split('\n')
        if len(lines) > 5:
            print(f"[SKIP] {label} already has {len(lines)-1} points")
            continue
    cmd = ["python3", "-u", SIONNA_SCRIPT, str(mcs), "1" if perfect else "0", str(out_file)]
    print(f"[RUN] {label}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd="/workspace/nr-link-simulator/scripts")
    print(f"[DONE] {label}")
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-300:]}")

print("All simulations complete!")
