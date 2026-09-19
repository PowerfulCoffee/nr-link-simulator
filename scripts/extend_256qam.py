#!/usr/bin/env python3
"""Extend 256QAM simulations to 40dB and add high-SNR points for other MCSs."""
import subprocess, os, sys, csv

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'
sionna_script = '/workspace/nr-link-simulator/scripts/sionna_tdl_crosscheck.py'
os.makedirs(results_dir, exist_ok=True)

DELAY_SPREAD = 100
MAX_DOPPLER = 10

def read_csv(path):
    data = {}
    if not os.path.exists(path):
        return data
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            snr = int(round(float(row['SINR_dB'])))
            data[snr] = {
                'blocks': int(row['Blocks']),
                'errors': int(row['Errors']),
                'bler': float(row['BLER'])
            }
    return data

def write_csv(path, data):
    with open(path, 'w') as f:
        f.write('SINR_dB,Blocks,Errors,BLER\n')
        for snr in sorted(data.keys()):
            d = data[snr]
            f.write(f"{snr},{d['blocks']},{d['errors']},{d['bler']:.6f}\n")

def run_cpp(mcs, ce, snr_start, snr_end, max_blocks, target_errors, name):
    out_file = os.path.join(results_dir, f'{name}.csv')
    perfect_csi = 1 if ce == 'ideal' else 0
    data = read_csv(out_file)
    
    needed = []
    for snr in range(snr_start, snr_end + 1):
        if snr not in data or data[snr]['blocks'] < max_blocks // 2:
            needed.append(snr)
    
    if not needed:
        print(f"  [SKIP] {name} already has {len(data)} points up to {max(data.keys())}dB", flush=True)
        return data
    
    # Run in segments of 5dB to avoid timeout
    seg_size = 5
    for seg_start in range(min(needed), max(needed) + 1, seg_size):
        seg_end = min(seg_start + seg_size - 1, max(needed))
        tmp_out = f'/tmp/cpp_{name}_{seg_start}_{seg_end}.csv'
        cmd = [
            cpp_bin, str(mcs), str(seg_start), str(seg_end), '1',
            str(perfect_csi), str(max_blocks), str(target_errors),
            '30', '4', '4', '1', str(DELAY_SPREAD), str(MAX_DOPPLER),
            '1', '0', '2', tmp_out
        ]
        print(f"  [RUN] C++ {name}: {seg_start}-{seg_end}dB, max_blocks={max_blocks}", flush=True)
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True,
                                   cwd='/workspace/nr-link-simulator', timeout=1800)
            seg_data = read_csv(tmp_out)
            for snr, d in seg_data.items():
                data[snr] = d
            write_csv(out_file, data)
            # Print results
            for snr in sorted(seg_data.keys()):
                print(f"    SNR={snr}dB: BLER={seg_data[snr]['bler']:.4f} ({seg_data[snr]['errors']}/{seg_data[snr]['blocks']})", flush=True)
        except subprocess.TimeoutExpired:
            print(f"    [TIMEOUT] at segment {seg_start}-{seg_end}", flush=True)
            if os.path.exists(tmp_out):
                seg_data = read_csv(tmp_out)
                for snr, d in seg_data.items():
                    data[snr] = d
                write_csv(out_file, data)
        except Exception as e:
            print(f"    [ERROR] {e}", flush=True)
    
    return data

def run_sionna(mcs, ce, snr_start, snr_end, max_blocks, target_errors, name):
    out_file = os.path.join(results_dir, f'{name}.csv')
    data = read_csv(out_file)
    
    needed = []
    for snr in range(snr_start, snr_end + 1):
        if snr not in data or data[snr]['blocks'] < max_blocks // 2:
            needed.append(snr)
    
    if not needed:
        print(f"  [SKIP] {name} already complete", flush=True)
        return data
    
    # Run in segments of 4dB for Sionna (slower)
    seg_size = 4
    for seg_start in range(min(needed), max(needed) + 1, seg_size):
        seg_end = min(seg_start + seg_size - 1, max(needed))
        cmd = [
            sys.executable, '-u', sionna_script,
            '--mcs', str(mcs),
            '--snr-start', str(seg_start),
            '--snr-end', str(seg_end),
            '--snr-step', '1',
            '--ce-type', ce,
            '--delay-spread', str(DELAY_SPREAD),
            '--doppler', str(MAX_DOPPLER),
            '--max-blocks', str(max_blocks),
            '--target-errors', str(target_errors),
            '--batch-size', '2',
            '--output', out_file
        ]
        print(f"  [RUN] Sionna {name}: {seg_start}-{seg_end}dB, max_blocks={max_blocks}", flush=True)
        try:
            proc = subprocess.Popen(cmd, cwd='/workspace/nr-link-simulator',
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, bufsize=1)
            for line in proc.stdout:
                line = line.strip()
                if 'SNR=' in line and 'BLER=' in line:
                    print(f"    {line}", flush=True)
            proc.wait()
            data = read_csv(out_file)
        except Exception as e:
            print(f"    [ERROR] {e}", flush=True)
    
    return data

print("="*60, flush=True)
print("Phase 1: C++ 256QAM MCS27 Ideal CE 24-40dB", flush=True)
print("="*60, flush=True)
# Use more blocks for statistical accuracy in waterfall
run_cpp(27, 'ideal', 24, 40, 50, 15, 'cpp_ideal_mcs27')

print("\n" + "="*60, flush=True)
print("Phase 2: C++ 256QAM MCS27 LS CE 28-40dB", flush=True)
print("="*60, flush=True)
run_cpp(27, 'ls', 28, 40, 50, 15, 'cpp_ls_mcs27')

print("\n" + "="*60, flush=True)
print("Phase 3: Sionna 256QAM MCS27 Ideal CE 31-40dB", flush=True)
print("="*60, flush=True)
run_sionna(27, 'ideal', 31, 40, 30, 10, 'sionna_ideal_mcs27')

print("\n" + "="*60, flush=True)
print("Phase 4: Sionna 256QAM MCS27 LS CE 31-40dB", flush=True)
print("="*60, flush=True)
run_sionna(27, 'ls', 31, 40, 30, 10, 'sionna_ls_mcs27')

# Also extend 64QAM LS to make sure it fully converges
print("\n" + "="*60, flush=True)
print("Phase 5: C++ 64QAM MCS15 LS CE 30-36dB (confirm convergence)", flush=True)
print("="*60, flush=True)
run_cpp(15, 'ls', 30, 36, 30, 10, 'cpp_ls_mcs15')

print("\n" + "="*60, flush=True)
print("All extended simulations complete!", flush=True)
print("="*60, flush=True)

# Print summary
for name in ['cpp_ideal_mcs27', 'cpp_ls_mcs27', 'sionna_ideal_mcs27', 'sionna_ls_mcs27']:
    path = os.path.join(results_dir, f'{name}.csv')
    data = read_csv(path)
    print(f"\n{name}:")
    for snr in sorted(data.keys()):
        if snr >= 24:
            d = data[snr]
            print(f"  SNR={snr:2d}dB: BLER={d['bler']:.4f} ({d['errors']:3d}/{d['blocks']:3d})")
