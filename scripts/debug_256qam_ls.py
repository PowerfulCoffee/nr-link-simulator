#!/usr/bin/env python3
"""Debug LS CE for 256QAM - run at 45dB with 100 blocks and print debug info."""
import subprocess, os, sys

# First, let's also check Sionna LS CE at high SNR to see if it has the same issue
# But first run C++ with 30 blocks at 40, 45, 50dB to confirm error floor
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'

for snr in [40, 45, 50]:
    tmp_out = f'/tmp/debug_ls_{snr}.csv'
    print(f"\nRunning C++ MCS27 LS CE at SNR={snr}dB, 30 blocks...", flush=True)
    cmd = [
        cpp_bin, '27', str(snr), str(snr), '1',
        '0', '30', '10',
        '30', '4', '4', '1', '100', '10',
        '1', '0', '2', tmp_out
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                               cwd='/workspace/nr-link-simulator', timeout=1800)
        if os.path.exists(tmp_out):
            with open(tmp_out) as f:
                lines = f.readlines()
                for line in lines:
                    print(f"  {line.strip()}", flush=True)
        if result.returncode != 0:
            print(f"  STDERR: {result.stderr[-500:]}", flush=True)
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT", flush=True)
        if os.path.exists(tmp_out):
            with open(tmp_out) as f:
                print(f"  Partial: {f.read()}", flush=True)
