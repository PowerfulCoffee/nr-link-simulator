import subprocess
import sys

sys.path.insert(0, '/workspace/nr-link-simulator/scripts')

configs = [
    {"mcs": 8,  "ce": "ideal", "snr_start": 6,  "snr_end": 16, "label": "16QAM Ideal"},
    {"mcs": 8,  "ce": "ls",    "snr_start": 10, "snr_end": 22, "label": "16QAM LS"},
    {"mcs": 15, "ce": "ideal", "snr_start": 14, "snr_end": 24, "label": "64QAM Ideal"},
    {"mcs": 15, "ce": "ls",    "snr_start": 18, "snr_end": 30, "label": "64QAM LS"},
]

for cfg in configs:
    mcs = cfg["mcs"]
    ce = cfg["ce"]
    s0 = cfg["snr_start"]
    s1 = cfg["snr_end"]
    label = cfg["label"]
    outfile = f"/workspace/nr-link-simulator/results/ce_compare/sionna_{ce}_mcs{mcs}.csv"
    print(f"\n=== Sionna {label} MCS{mcs} ({ce}) SNR=[{s0},{s1}] ===", flush=True)
    
    cmd = [
        sys.executable,
        "/workspace/nr-link-simulator/scripts/sionna_tdl_crosscheck.py",
        "--mcs", str(mcs),
        "--snr-start", str(s0),
        "--snr-end", str(s1),
        "--snr-step", "1",
        "--max-blocks", "20",
        "--target-errors", "10",
        "--batch-size", "4",
        "--ce-type", ce,
        "--delay-spread", "15",
        "--doppler", "5",
        "--output", outfile
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd="/workspace/nr-link-simulator")
    for line in result.stdout.split('\n'):
        if 'SNR=' in line or 'Results saved' in line or 'MCS=' in line:
            print(line, flush=True)
    if result.returncode != 0:
        print(f"STDERR: {result.stderr[-500:]}", flush=True)

print("\n=== ALL SIONNA DONE ===", flush=True)
