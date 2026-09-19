import subprocess
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RESULTS_DIR = "/workspace/nr-link-simulator/results/ce_compare"
SCRIPT = "/workspace/nr-link-simulator/scripts/sionna_tdl_crosscheck.py"

CONFIGS = [
    # MCS4 LS (QPSK) - re-run with fixes
    {"mcs": 4, "ce": "ls", "snr_start": 4, "snr_end": 14, "output": "sionna_ls_mcs4.csv",
     "max_blocks": 50, "batch_size": 5},
    # MCS8 LS (16QAM)
    {"mcs": 8, "ce": "ls", "snr_start": 10, "snr_end": 22, "output": "sionna_ls_mcs8.csv",
     "max_blocks": 50, "batch_size": 5},
    # MCS15 Ideal (64QAM)
    {"mcs": 15, "ce": "ideal", "snr_start": 14, "snr_end": 24, "output": "sionna_ideal_mcs15.csv",
     "max_blocks": 50, "batch_size": 2},
    # MCS15 LS (64QAM)
    {"mcs": 15, "ce": "ls", "snr_start": 18, "snr_end": 30, "output": "sionna_ls_mcs15.csv",
     "max_blocks": 50, "batch_size": 2},
]

os.makedirs(RESULTS_DIR, exist_ok=True)

for cfg in CONFIGS:
    out_path = os.path.join(RESULTS_DIR, cfg["output"])
    if os.path.exists(out_path) and cfg["mcs"] != 4:
        print(f"SKIP {cfg['output']} (exists)")
        continue
    
    print(f"\n{'='*60}")
    print(f"Running MCS={cfg['mcs']} CE={cfg['ce']} SNR={cfg['snr_start']}-{cfg['snr_end']}dB")
    print(f"{'='*60}")
    
    cmd = [
        sys.executable, SCRIPT,
        "--mcs", str(cfg["mcs"]),
        "--ce-type", cfg["ce"],
        "--snr-start", str(cfg["snr_start"]),
        "--snr-end", str(cfg["snr_end"]),
        "--snr-step", "1",
        "--max-blocks", str(cfg["max_blocks"]),
        "--batch-size", str(cfg["batch_size"]),
        "--output", out_path,
    ]
    
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in proc.stdout:
        print(line, end='', flush=True)
    proc.wait()
    if proc.returncode != 0:
        print(f"ERROR: exited with code {proc.returncode}")
    else:
        print(f"DONE: {out_path}")

print("\nAll simulations complete!")
