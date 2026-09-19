import subprocess
import os
import sys

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'
sionna_script = '/workspace/nr-link-simulator/scripts/sionna_tdl_crosscheck.py'
os.makedirs(results_dir, exist_ok=True)

DELAY_SPREAD = 100
MAX_DOPPLER = 10
MAX_BLOCKS = 50
TARGET_ERRORS = 20
BATCH_SIZE = 2
SNR_STEP = 1

# Common C++ args positions:
# mcs, sinr_start, sinr_end, step, perfect_csi, max_blocks, target_errors,
# n_rb, n_rx, n_layers, channel_type, delay_spread_ns, max_doppler,
# dmrs_add_pos, dmrs_duration, mcs_table, outfile

def run_cpp_sim(mcs, ce_type, snr_start, snr_end, name):
    out_file = os.path.join(results_dir, f"{name}.csv")
    if os.path.exists(out_file) and os.path.getsize(out_file) > 50:
        with open(out_file, 'r') as f:
            lines = f.readlines()
        if len(lines) > 5:
            print(f"  [SKIP] {name} exists")
            return True
    
    perfect_csi = 1 if ce_type == 'ideal' else 0
    
    cmd = [
        cpp_bin, str(mcs), str(snr_start), str(snr_end), str(SNR_STEP),
        str(perfect_csi), str(MAX_BLOCKS), str(TARGET_ERRORS),
        '30', '4', '4', '1', str(DELAY_SPREAD), str(MAX_DOPPLER),
        '1', '0', '2', out_file
    ]
    
    print(f"  [RUN] C++ {name}: MCS{mcs} {ce_type} SNR {snr_start}-{snr_end}dB", flush=True)
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True,
                              cwd='/workspace/nr-link-simulator', timeout=600)
        print(f"  [DONE] {name}", flush=True)
        return True
    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT] {name}", flush=True)
        return False
    except subprocess.CalledProcessError as e:
        print(f"  [ERROR] {name}: {e}", flush=True)
        return False

def run_sionna_sim(mcs, ce_type, snr_start, snr_end, name):
    out_file = os.path.join(results_dir, f"{name}.csv")
    if os.path.exists(out_file) and os.path.getsize(out_file) > 50:
        with open(out_file, 'r') as f:
            lines = f.readlines()
        if len(lines) > 5:
            print(f"  [SKIP] {name} exists")
            return True
    
    cmd = [
        'python', sionna_script,
        '--mcs', str(mcs),
        '--snr-start', str(snr_start),
        '--snr-end', str(snr_end),
        '--snr-step', str(SNR_STEP),
        '--ce-type', ce_type,
        '--delay-spread', str(DELAY_SPREAD),
        '--doppler', str(MAX_DOPPLER),
        '--max-blocks', str(MAX_BLOCKS),
        '--target-errors', str(TARGET_ERRORS),
        '--batch-size', str(BATCH_SIZE),
        '--output', out_file
    ]
    
    print(f"  [RUN] Sionna {name}: MCS{mcs} {ce_type} SNR {snr_start}-{snr_end}dB", flush=True)
    try:
        process = subprocess.Popen(cmd, cwd='/workspace/nr-link-simulator',
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, bufsize=1)
        for line in process.stdout:
            if 'SNR=' in line and 'BLER=' in line:
                print(f"    {line.strip()}", flush=True)
        process.wait()
        if process.returncode == 0:
            print(f"  [DONE] {name}", flush=True)
            return True
        else:
            print(f"  [ERROR] {name} returned {process.returncode}", flush=True)
            return False
    except Exception as e:
        print(f"  [ERROR] {name}: {e}", flush=True)
        return False

# All simulation configs: (mcs, ce_type, snr_start, snr_end, cpp_name, sionna_name)
sim_configs = [
    # QPSK MCS4
    (4, 'ideal', -2, 12, 'cpp_ideal_mcs4', 'sionna_ideal_mcs4'),
    (4, 'ls', 4, 22, 'cpp_ls_mcs4', 'sionna_ls_mcs4'),
    # 16QAM MCS8
    (8, 'ideal', 8, 22, 'cpp_ideal_mcs8', 'sionna_ideal_mcs8'),
    (8, 'ls', 14, 30, 'cpp_ls_mcs8', 'sionna_ls_mcs8'),
    # 64QAM MCS15
    (15, 'ideal', 16, 32, 'cpp_ideal_mcs15', 'sionna_ideal_mcs15'),
    (15, 'ls', 22, 40, 'cpp_ls_mcs15', 'sionna_ls_mcs15'),
    # 256QAM MCS27
    (27, 'ideal', 24, 40, 'cpp_ideal_mcs27', 'sionna_ideal_mcs27'),
    (27, 'ls', 30, 48, 'cpp_ls_mcs27', 'sionna_ls_mcs27'),
]

print("="*60, flush=True)
print("Starting all C++ simulations", flush=True)
print("="*60, flush=True)

for mcs, ce, s_start, s_end, cpp_name, sionna_name in sim_configs:
    run_cpp_sim(mcs, ce, s_start, s_end, cpp_name)

print("\n" + "="*60, flush=True)
print("Starting all Sionna simulations", flush=True)
print("="*60, flush=True)

for mcs, ce, s_start, s_end, cpp_name, sionna_name in sim_configs:
    run_sionna_sim(mcs, ce, s_start, s_end, sionna_name)

print("\n" + "="*60, flush=True)
print("All simulations complete!", flush=True)
print("="*60, flush=True)

import glob
print("\nGenerated files:", flush=True)
for f in sorted(glob.glob(os.path.join(results_dir, '*.csv'))):
    size = os.path.getsize(f)
    print(f"  {os.path.basename(f)}: {size} bytes", flush=True)
