#!/usr/bin/env python3
"""Extend C++ LS CE 256QAM to 45dB and run Sionna 256QAM to 40dB."""
import subprocess, os, sys, csv

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'
sionna_script = '/workspace/nr-link-simulator/scripts/sionna_tdl_crosscheck.py'

def read_csv(path):
    data = {}
    if not os.path.exists(path):
        return data
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                snr = int(round(float(row['SINR_dB'])))
                data[snr] = {
                    'blocks': int(row['Blocks']),
                    'errors': int(row['Errors']),
                    'bler': float(row['BLER'])
                }
            except:
                pass
    return data

def write_csv(path, data):
    with open(path, 'w') as f:
        f.write('SINR_dB,Blocks,Errors,BLER\n')
        for snr in sorted(data.keys()):
            d = data[snr]
            f.write(f"{snr},{d['blocks']},{d['errors']},{d['bler']:.6f}\n")

def run_cpp_point(mcs, ce, snr, max_blocks, target_errors):
    perfect_csi = 1 if ce == 'ideal' else 0
    tmp_out = f'/tmp/cpp_s_{mcs}_{ce}_{snr}.csv'
    cmd = [
        cpp_bin, str(mcs), str(snr), str(snr), '1',
        str(perfect_csi), str(max_blocks), str(target_errors),
        '30', '4', '4', '1', '100', '10',
        '1', '0', '2', tmp_out
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True,
                       cwd='/workspace/nr-link-simulator', timeout=600)
        data = read_csv(tmp_out)
        return data.get(snr, None)
    except subprocess.TimeoutExpired:
        if os.path.exists(tmp_out):
            data = read_csv(tmp_out)
            return data.get(snr, None)
        return None
    except Exception as e:
        print(f"    ERROR: {e}", flush=True)
        return None

# Phase 1: C++ 256QAM LS CE extend to 45dB (15 blocks for speed)
print("="*60, flush=True)
print("Phase 1: C++ 256QAM LS CE 41-45dB (15 blocks)", flush=True)
print("="*60, flush=True)

ls_file = os.path.join(results_dir, 'cpp_ls_mcs27.csv')
ls_data = read_csv(ls_file)

for snr in range(41, 46):
    result = run_cpp_point(27, 'ls', snr, 15, 5)
    if result:
        ls_data[snr] = result
        write_csv(ls_file, ls_data)
        print(f"  SNR={snr}dB: BLER={result['bler']:.4f} ({result['errors']}/{result['blocks']})", flush=True)
    else:
        print(f"  SNR={snr}dB: FAILED", flush=True)

# Also extend Ideal CE to 45dB for completeness
print("\n" + "="*60, flush=True)
print("Phase 2: C++ 256QAM Ideal CE 41-45dB (10 blocks, confirm 0)", flush=True)
print("="*60, flush=True)

ideal_file = os.path.join(results_dir, 'cpp_ideal_mcs27.csv')
ideal_data = read_csv(ideal_file)

for snr in range(41, 46):
    result = run_cpp_point(27, 'ideal', snr, 10, 1)
    if result:
        ideal_data[snr] = result
        write_csv(ideal_file, ideal_data)
        print(f"  SNR={snr}dB: BLER={result['bler']:.4f} ({result['errors']}/{result['blocks']})", flush=True)

# Phase 3: Sionna 256QAM Ideal CE 31-40dB
print("\n" + "="*60, flush=True)
print("Phase 3: Sionna 256QAM Ideal CE 31-40dB", flush=True)
print("="*60, flush=True)

for seg_start in range(31, 41, 5):
    seg_end = min(seg_start + 4, 40)
    out_file = os.path.join(results_dir, 'sionna_ideal_mcs27.csv')
    cmd = [
        sys.executable, '-u', sionna_script,
        '--mcs', '27', '--snr-start', str(seg_start), '--snr-end', str(seg_end),
        '--snr-step', '1', '--ce-type', 'ideal',
        '--delay-spread', '100', '--doppler', '10',
        '--max-blocks', '12', '--target-errors', '4', '--batch-size', '2',
        '--output', out_file
    ]
    print(f"  Running Sionna Ideal: {seg_start}-{seg_end}dB...", flush=True)
    try:
        proc = subprocess.Popen(cmd, cwd='/workspace/nr-link-simulator',
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, bufsize=1)
        for line in proc.stdout:
            line = line.strip()
            if 'SNR=' in line and 'BLER=' in line:
                print(f"    {line}", flush=True)
        proc.wait()
    except Exception as e:
        print(f"    ERROR: {e}", flush=True)

# Phase 4: Sionna 256QAM LS CE 31-40dB
print("\n" + "="*60, flush=True)
print("Phase 4: Sionna 256QAM LS CE 31-40dB", flush=True)
print("="*60, flush=True)

for seg_start in range(31, 41, 5):
    seg_end = min(seg_start + 4, 40)
    out_file = os.path.join(results_dir, 'sionna_ls_mcs27.csv')
    cmd = [
        sys.executable, '-u', sionna_script,
        '--mcs', '27', '--snr-start', str(seg_start), '--snr-end', str(seg_end),
        '--snr-step', '1', '--ce-type', 'ls',
        '--delay-spread', '100', '--doppler', '10',
        '--max-blocks', '12', '--target-errors', '4', '--batch-size', '2',
        '--output', out_file
    ]
    print(f"  Running Sionna LS: {seg_start}-{seg_end}dB...", flush=True)
    try:
        proc = subprocess.Popen(cmd, cwd='/workspace/nr-link-simulator',
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, bufsize=1)
        for line in proc.stdout:
            line = line.strip()
            if 'SNR=' in line and 'BLER=' in line:
                print(f"    {line}", flush=True)
        proc.wait()
    except Exception as e:
        print(f"    ERROR: {e}", flush=True)

# Print final summary
print("\n" + "="*60, flush=True)
print("Final Summary", flush=True)
print("="*60, flush=True)
for name in ['cpp_ideal_mcs27', 'cpp_ls_mcs27', 'sionna_ideal_mcs27', 'sionna_ls_mcs27']:
    path = os.path.join(results_dir, f'{name}.csv')
    data = read_csv(path)
    print(f"\n{name}:")
    for snr in sorted(data.keys()):
        if snr >= 28:
            d = data[snr]
            print(f"  SNR={snr:2d}dB: BLER={d['bler']:.4f} ({d['errors']:3d}/{d['blocks']:3d})")
