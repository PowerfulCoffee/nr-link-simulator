#!/usr/bin/env python3
"""Quick C++ scan to find 256QAM waterfall region, then Sionna follow-up."""
import subprocess, os, sys, csv

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'

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
    """Run a single SNR point."""
    perfect_csi = 1 if ce == 'ideal' else 0
    tmp_out = f'/tmp/cpp_single_{mcs}_{ce}_{snr}.csv'
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
        # Read whatever we have
        if os.path.exists(tmp_out):
            data = read_csv(tmp_out)
            return data.get(snr, None)
        return None
    except Exception as e:
        print(f"    ERROR at SNR={snr}: {e}", flush=True)
        return None

# Phase 1: Quick scan of C++ 256QAM Ideal CE 28-40dB with 20 blocks
print("="*60, flush=True)
print("Phase 1: C++ 256QAM Ideal CE quick scan (20 blocks, 5 errors)", flush=True)
print("="*60, flush=True)

ideal_file = os.path.join(results_dir, 'cpp_ideal_mcs27.csv')
ideal_data = read_csv(ideal_file)

for snr in range(28, 41):
    if snr in ideal_data and ideal_data[snr]['blocks'] >= 15 and ideal_data[snr]['bler'] < 0.1:
        print(f"  SNR={snr}dB: BLER={ideal_data[snr]['bler']:.4f} (existing)", flush=True)
        continue
    if snr in ideal_data and ideal_data[snr]['blocks'] >= 20:
        print(f"  SNR={snr}dB: BLER={ideal_data[snr]['bler']:.4f} (existing)", flush=True)
        continue
    
    result = run_cpp_point(27, 'ideal', snr, 20, 5)
    if result:
        ideal_data[snr] = result
        write_csv(ideal_file, ideal_data)
        print(f"  SNR={snr}dB: BLER={result['bler']:.4f} ({result['errors']}/{result['blocks']})", flush=True)
    else:
        print(f"  SNR={snr}dB: FAILED", flush=True)

# Phase 2: C++ 256QAM LS CE 30-40dB
print("\n" + "="*60, flush=True)
print("Phase 2: C++ 256QAM LS CE quick scan (20 blocks, 5 errors)", flush=True)
print("="*60, flush=True)

ls_file = os.path.join(results_dir, 'cpp_ls_mcs27.csv')
ls_data = read_csv(ls_file)

for snr in range(30, 41):
    if snr in ls_data and ls_data[snr]['blocks'] >= 15 and ls_data[snr]['bler'] < 0.1:
        print(f"  SNR={snr}dB: BLER={ls_data[snr]['bler']:.4f} (existing)", flush=True)
        continue
    if snr in ls_data and ls_data[snr]['blocks'] >= 20:
        print(f"  SNR={snr}dB: BLER={ls_data[snr]['bler']:.4f} (existing)", flush=True)
        continue
    
    result = run_cpp_point(27, 'ls', snr, 20, 5)
    if result:
        ls_data[snr] = result
        write_csv(ls_file, ls_data)
        print(f"  SNR={snr}dB: BLER={result['bler']:.4f} ({result['errors']}/{result['blocks']})", flush=True)
    else:
        print(f"  SNR={snr}dB: FAILED", flush=True)

# Print summary
print("\n" + "="*60, flush=True)
print("Summary", flush=True)
print("="*60, flush=True)
print(f"\nMCS27 Ideal CE:")
for snr in sorted(ideal_data.keys()):
    if snr >= 24:
        d = ideal_data[snr]
        print(f"  SNR={snr:2d}dB: BLER={d['bler']:.4f} ({d['errors']:3d}/{d['blocks']:3d})")

print(f"\nMCS27 LS CE:")
for snr in sorted(ls_data.keys()):
    if snr >= 24:
        d = ls_data[snr]
        print(f"  SNR={snr:2d}dB: BLER={d['bler']:.4f} ({d['errors']:3d}/{d['blocks']:3d})")
