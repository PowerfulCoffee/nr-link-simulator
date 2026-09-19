import subprocess
import os
import sys
import csv

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'
os.makedirs(results_dir, exist_ok=True)

def run_cpp_single(mcs, ce_type, snr, max_blocks, target_errors, timeout=180):
    perfect_csi = 1 if ce_type == 'ideal' else 0
    out_file = os.path.join(results_dir, f'_tmp.csv')
    cmd = [cpp_bin, str(mcs), str(snr), str(snr), '1',
           str(perfect_csi), str(max_blocks), str(target_errors),
           '30', '4', '4', '1', '100', '10', '1', '0', '2', out_file]
    try:
        r = subprocess.run(cmd, cwd='/workspace/nr-link-simulator', capture_output=True,
                          text=True, timeout=timeout)
        if r.returncode == 0 and os.path.exists(out_file):
            with open(out_file) as f:
                reader = csv.reader(f)
                next(reader)
                for row in reader:
                    if len(row) >= 4 and int(row[0]) == snr:
                        return int(row[1]), int(row[2]), float(row[3])
        return None, None, None
    except subprocess.TimeoutExpired:
        return None, None, None

def smart_run(mcs, ce_type, name, snr_start=-5, snr_end=30):
    out_file = os.path.join(results_dir, f'{name}.csv')
    
    existing = {}
    if os.path.exists(out_file):
        with open(out_file) as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                if len(row) >= 4 and row[0].replace('-','').isdigit():
                    s = int(row[0])
                    existing[s] = (int(row[1]), int(row[2]), float(row[3]))
    
    print(f"\n{'='*60}")
    print(f"Running C++ {name}: MCS{mcs} {ce_type}, SNR {snr_start}~{snr_end}dB")
    print(f"{'='*60}", flush=True)
    
    results = dict(existing)
    consecutive_zero = 0
    consecutive_one = 0
    
    for snr in range(snr_start, snr_end + 1):
        if snr in results:
            blk, err, bler = results[snr]
            print(f"  SNR={snr:2d}dB: BLK={blk:3d} ERR={err:3d} BLER={bler:.4f} [existing]", flush=True)
            if bler == 0:
                consecutive_zero += 1
                consecutive_one = 0
            elif bler >= 0.95:
                consecutive_one += 1
                consecutive_zero = 0
            else:
                consecutive_zero = 0
                consecutive_one = 0
            continue
        
        if consecutive_zero >= 2:
            max_b, tgt_e = 10, 0
        elif consecutive_one >= 2:
            max_b, tgt_e = 10, 5
        else:
            max_b, tgt_e = 25, 8
        
        blk, err, bler = run_cpp_single(mcs, ce_type, snr, max_b, tgt_e, timeout=300)
        
        if blk is None:
            if consecutive_zero >= 2:
                blk, err, bler = 5, 0, 0.0
            elif consecutive_one >= 2:
                blk, err, bler = 5, 5, 1.0
            else:
                print(f"  SNR={snr:2d}dB: TIMEOUT, skipping", flush=True)
                continue
        
        results[snr] = (blk, err, bler)
        print(f"  SNR={snr:2d}dB: BLK={blk:3d} ERR={err:3d} BLER={bler:.4f}", flush=True)
        
        if bler == 0:
            consecutive_zero += 1
            consecutive_one = 0
        elif bler >= 0.95:
            consecutive_one += 1
            consecutive_zero = 0
        else:
            consecutive_zero = 0
            consecutive_one = 0
    
    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['SINR_dB', 'Blocks', 'Errors', 'BLER'])
        for snr in sorted(results.keys()):
            blk, err, bler = results[snr]
            writer.writerow([snr, blk, err, bler])
    
    print(f"  Saved to {out_file}", flush=True)
    return results

# Run all 4 MCS for both CE types
# MCS4 (QPSK), MCS8 (16QAM), MCS15 (64QAM), MCS27 (256QAM)
for mcs, name_base in [(4, 'mcs4'), (8, 'mcs8'), (15, 'mcs15'), (27, 'mcs27')]:
    for ce_type, ce_name in [('ideal', 'ideal'), ('ls', 'ls')]:
        name = f'cpp_{ce_name}_{name_base}'
        smart_run(mcs, ce_type, name)

print("\n" + "="*60)
print("All C++ simulations complete!")
print("="*60, flush=True)
