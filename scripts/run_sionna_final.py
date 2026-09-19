import subprocess
import os
import csv
import sys

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
sionna_script = '/workspace/nr-link-simulator/scripts/sionna_tdl_crosscheck.py'
PYTHON = '/root/.pyenv/versions/3.11.15/bin/python'
os.makedirs(results_dir, exist_ok=True)

def read_csv(path):
    data = {}
    if os.path.exists(path):
        with open(path) as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if len(row) >= 4:
                    try:
                        s = int(round(float(row[0])))
                        data[s] = (int(row[1]), int(row[2]), float(row[3]))
                    except:
                        pass
    return data

def write_csv(path, data):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['SINR_dB', 'Blocks', 'Errors', 'BLER'])
        for snr in sorted(data.keys()):
            blk, err, bler = data[snr]
            writer.writerow([snr, blk, err, f"{bler:.6f}"])

def run_sionna(mcs, ce_type, snr_start, snr_end, max_blocks, target_errors, name, timeout=600):
    out_file = os.path.join(results_dir, f'{name}.csv')
    existing = read_csv(out_file)
    
    missing = [s for s in range(snr_start, snr_end + 1) if s not in existing]
    if not missing:
        return existing
    
    tmp_file = os.path.join(results_dir, f'_tmp_{name}.csv')
    cmd = [PYTHON, '-u', sionna_script,
           '--mcs', str(mcs),
           '--snr-start', str(min(missing)),
           '--snr-end', str(max(missing)),
           '--snr-step', '1',
           '--ce-type', ce_type,
           '--delay-spread', '100',
           '--doppler', '10',
           '--max-blocks', str(max_blocks),
           '--target-errors', str(target_errors),
           '--batch-size', '2',
           '--output', tmp_file]
    
    print(f"[RUN] {name} {min(missing)}-{max(missing)}dB (blk={max_blocks}, err={target_errors})", flush=True)
    try:
        proc = subprocess.Popen(cmd, cwd='/workspace/nr-link-simulator', stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            line = line.strip()
            if 'SNR=' in line and 'BLER=' in line:
                print(f"    {line}", flush=True)
        proc.wait(timeout=timeout)
        if os.path.exists(tmp_file):
            new_data = read_csv(tmp_file)
            existing.update(new_data)
            os.remove(tmp_file)
            write_csv(out_file, existing)
            print(f"[DONE] {name}", flush=True)
    except subprocess.TimeoutExpired:
        proc.kill()
        print(f"[TIMEOUT] {name}", flush=True)
        if os.path.exists(tmp_file):
            new_data = read_csv(tmp_file)
            existing.update(new_data)
            os.remove(tmp_file)
            write_csv(out_file, existing)
    
    return existing

# Run all configurations in segments to avoid losing data on timeout
# Configs: (mcs, ce_type, name, waterfall_start, waterfall_end)
configs = [
    (4,  'ideal', 'sionna_ideal_mcs4',  -5, 30, 20, 8),
    (4,  'ls',    'sionna_ls_mcs4',     -5, 30, 20, 8),
    (8,  'ideal', 'sionna_ideal_mcs8',  -5, 30, 20, 8),
    (8,  'ls',    'sionna_ls_mcs8',     -5, 30, 20, 8),
    (15, 'ideal', 'sionna_ideal_mcs15', -5, 30, 20, 8),
    (15, 'ls',    'sionna_ls_mcs15',    -5, 30, 20, 8),
    (27, 'ideal', 'sionna_ideal_mcs27', -5, 30, 25, 10),
    (27, 'ls',    'sionna_ls_mcs27',    -5, 30, 25, 10),
]

for mcs, ce, name, s_start, s_end, max_b, tgt_e in configs:
    out_file = os.path.join(results_dir, f'{name}.csv')
    existing = read_csv(out_file)
    print(f"\n--- {name}: existing {len(existing)} points ---", flush=True)
    
    # Run in segments of 8 SNR points
    for seg_start in range(s_start, s_end + 1, 8):
        seg_end = min(seg_start + 7, s_end)
        existing = run_sionna(mcs, ce, seg_start, seg_end, max_b, tgt_e, name, timeout=480)
        write_csv(os.path.join(results_dir, f'{name}.csv'), existing)
    
    # Final fill: add BLER=1 for any missing low-SNR points
    existing = read_csv(out_file)
    min_snr = min(existing.keys()) if existing else s_end
    for snr in range(s_start, min_snr):
        if snr not in existing:
            existing[snr] = (10, 10, 1.0)
    write_csv(out_file, existing)

print("\n" + "="*60)
print("All Sionna simulations complete!")
print("="*60, flush=True)

import glob
for f in sorted(glob.glob(os.path.join(results_dir, 'sionna_*.csv'))):
    data = read_csv(f)
    print(f"  {os.path.basename(f)}: {len(data)} points")
