#!/usr/bin/env python3
"""Verify C++ LS CE error floor at high SNR with more blocks."""
import subprocess, os

cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'
tmp_out = '/tmp/verify_ls_floor.csv'

# Run at 45dB with 50 blocks to confirm error floor
print("Running C++ MCS27 LS CE at 45dB with 50 blocks...", flush=True)
cmd = [
    cpp_bin, '27', '45', '45', '1',
    '0', '50', '20',
    '30', '4', '4', '1', '100', '10',
    '1', '0', '2', tmp_out
]
result = subprocess.run(cmd, capture_output=True, text=True, cwd='/workspace/nr-link-simulator', timeout=1200)
print("STDOUT:", result.stdout[-500:] if result.stdout else "")
print("STDERR:", result.stderr[-500:] if result.stderr else "")
print("Return code:", result.returncode)

if os.path.exists(tmp_out):
    with open(tmp_out) as f:
        print("\nResult:", f.read())

# Also run Ideal CE at same point for comparison
print("\n\nRunning C++ MCS27 Ideal CE at 45dB with 30 blocks...", flush=True)
tmp_out2 = '/tmp/verify_ideal.csv'
cmd2 = [
    cpp_bin, '27', '45', '45', '1',
    '1', '30', '1',
    '30', '4', '4', '1', '100', '10',
    '1', '0', '2', tmp_out2
]
result2 = subprocess.run(cmd2, capture_output=True, text=True, cwd='/workspace/nr-link-simulator', timeout=600)
if os.path.exists(tmp_out2):
    with open(tmp_out2) as f:
        print("Result:", f.read())
