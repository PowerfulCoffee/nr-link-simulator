import subprocess
import os

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'
os.makedirs(results_dir, exist_ok=True)

MAX_BLOCKS = 20
TARGET_ERRORS = 5
SNR_STEP = 1
TIMEOUT = 1200

def run_cpp(mcs, ce_type, snr_start, snr_end, name):
    out_file = os.path.join(results_dir, f"{name}.csv")
    if os.path.exists(out_file):
        with open(out_file) as f:
            lines = f.readlines()
        if len(lines) > 1:
            print(f"[SKIP] {name}")
            return True
    
    perfect_csi = 1 if ce_type == 'ideal' else 0
    cmd = [cpp_bin, str(mcs), str(snr_start), str(snr_end), str(SNR_STEP),
           str(perfect_csi), str(MAX_BLOCKS), str(TARGET_ERRORS),
           '30', '4', '4', '1', '100', '10', '1', '0', '2', out_file]
    print(f"[RUN] {name} MCS{mcs} {ce_type} {snr_start}-{snr_end}dB", flush=True)
    try:
        r = subprocess.run(cmd, cwd='/workspace/nr-link-simulator', capture_output=True,
                          text=True, timeout=TIMEOUT)
        if r.returncode == 0:
            print(f"[DONE] {name}")
            return True
        else:
            print(f"[FAIL] {name}: {r.stderr[:200]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"[TIMEOUT] {name}")
        return False

# 64QAM MCS15 - tighter ranges based on 16QAM results
run_cpp(15, 'ideal', 18, 28, 'cpp_ideal_mcs15')
run_cpp(15, 'ls',    24, 34, 'cpp_ls_mcs15')
# 256QAM MCS27
run_cpp(27, 'ideal', 26, 36, 'cpp_ideal_mcs27')
run_cpp(27, 'ls',    32, 42, 'cpp_ls_mcs27')

print("\nC++ remaining done!")
import glob
for f in sorted(glob.glob(os.path.join(results_dir, 'cpp_*.csv'))):
    with open(f) as fh:
        lines = fh.readlines()
    print(f"  {os.path.basename(f)}: {len(lines)-1} data points")
