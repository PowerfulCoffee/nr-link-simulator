import subprocess
import os
import csv

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'
os.makedirs(results_dir, exist_ok=True)

def run_cpp_range(mcs, ce_type, snr_start, snr_end, max_blocks, target_errors, name, timeout=900):
    out_file = os.path.join(results_dir, f'{name}.csv')
    perfect_csi = 1 if ce_type == 'ideal' else 0
    
    existing = {}
    if os.path.exists(out_file):
        with open(out_file) as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if len(row) >= 4:
                    try:
                        existing[int(row[0])] = (int(row[1]), int(row[2]), float(row[3]))
                    except:
                        pass
    
    need_run = any(s not in existing for s in range(snr_start, snr_end + 1))
    if not need_run:
        print(f"[SKIP] {name} {snr_start}-{snr_end}dB", flush=True)
        return existing
    
    tmp_file = os.path.join(results_dir, f'_tmp.csv')
    cmd = [cpp_bin, str(mcs), str(snr_start), str(snr_end), '1',
           str(perfect_csi), str(max_blocks), str(target_errors),
           '30', '4', '4', '1', '100', '10', '1', '0', '2', tmp_file]
    
    print(f"[RUN] {name} MCS{mcs} {ce_type} {snr_start}-{snr_end}dB (blk={max_blocks}, err={target_errors})", flush=True)
    try:
        proc = subprocess.Popen(cmd, cwd='/workspace/nr-link-simulator', stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            line = line.strip()
            if 'SINR' in line and 'BLER' in line and ('dB' in line or any(c.isdigit() for c in line[:10])):
                print(f"    {line}", flush=True)
        proc.wait(timeout=timeout)
        if os.path.exists(tmp_file):
            with open(tmp_file) as f:
                reader = csv.reader(f)
                next(reader, None)
                for row in reader:
                    if len(row) >= 4:
                        try:
                            existing[int(row[0])] = (int(row[1]), int(row[2]), float(row[3]))
                        except:
                            pass
            os.remove(tmp_file)
            print(f"[DONE] {name}", flush=True)
    except subprocess.TimeoutExpired:
        proc.kill()
        print(f"[TIMEOUT] {name}", flush=True)
        if os.path.exists(tmp_file):
            with open(tmp_file) as f:
                reader = csv.reader(f)
                next(reader, None)
                for row in reader:
                    if len(row) >= 4:
                        try:
                            existing[int(row[0])] = (int(row[1]), int(row[2]), float(row[3]))
                        except:
                            pass
            os.remove(tmp_file)
    
    return existing

def save_csv(name, data, fill_start=-5, fill_end=30):
    out_file = os.path.join(results_dir, f'{name}.csv')
    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['SINR_dB', 'Blocks', 'Errors', 'BLER'])
        for snr in range(fill_start, fill_end + 1):
            if snr in data:
                blk, err, bler = data[snr]
                writer.writerow([snr, blk, err, f"{bler:.6f}"])

# ====== QPSK MCS4 ======
# Ideal: -2~12dB exists (50blk/20err in waterfall)
data4i = run_cpp_range(4, 'ideal', -5, -3, 10, 5, 'cpp_ideal_mcs4')
data4i = run_cpp_range(4, 'ideal', 13, 30, 8, 0, 'cpp_ideal_mcs4')
save_csv('cpp_ideal_mcs4', data4i)

# LS: 6~18dB exists
data4l = run_cpp_range(4, 'ls', -5, 5, 10, 5, 'cpp_ls_mcs4')
data4l = run_cpp_range(4, 'ls', 19, 30, 8, 0, 'cpp_ls_mcs4')
save_csv('cpp_ls_mcs4', data4l)

# ====== 16QAM MCS8 ======
# Ideal: 8~20dB exists
data8i = run_cpp_range(8, 'ideal', -5, 7, 10, 5, 'cpp_ideal_mcs8')
data8i = run_cpp_range(8, 'ideal', 21, 30, 8, 0, 'cpp_ideal_mcs8')
save_csv('cpp_ideal_mcs8', data8i)

# LS: 14~26dB exists
data8l = run_cpp_range(8, 'ls', -5, 13, 10, 5, 'cpp_ls_mcs8')
data8l = run_cpp_range(8, 'ls', 27, 30, 8, 0, 'cpp_ls_mcs8')
save_csv('cpp_ls_mcs8', data8l)

# ====== 64QAM MCS15 ======
# Ideal: 18~28dB exists (20blk)
data15i = run_cpp_range(15, 'ideal', -5, 17, 10, 5, 'cpp_ideal_mcs15')
data15i = run_cpp_range(15, 'ideal', 29, 30, 8, 0, 'cpp_ideal_mcs15')
save_csv('cpp_ideal_mcs15', data15i)

# LS: 24~30dB exists (20blk)
data15l = run_cpp_range(15, 'ls', -5, 23, 10, 5, 'cpp_ls_mcs15')
save_csv('cpp_ls_mcs15', data15l)

# ====== 256QAM MCS27 ======
# Test showed BLER=0.5 at 30dB Ideal, so need 24-30dB range
data27i = run_cpp_range(27, 'ideal', -5, 23, 10, 5, 'cpp_ideal_mcs27')
data27i = run_cpp_range(27, 'ideal', 24, 30, 30, 8, 'cpp_ideal_mcs27')
save_csv('cpp_ideal_mcs27', data27i)

# LS: will be even higher - start from 28dB to 30dB
data27l = run_cpp_range(27, 'ls', -5, 27, 10, 5, 'cpp_ls_mcs27')
data27l = run_cpp_range(27, 'ls', 28, 30, 30, 8, 'cpp_ls_mcs27')
save_csv('cpp_ls_mcs27', data27l)

print("\n" + "="*60)
print("C++ simulations complete!")
print("="*60, flush=True)

import glob
for f in sorted(glob.glob(os.path.join(results_dir, 'cpp_*.csv'))):
    with open(f) as fh:
        lines = fh.readlines()
    pts = [l for l in lines[1:] if l.strip()]
    print(f"  {os.path.basename(f)}: {len(pts)} points")
