import subprocess
import os

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
script_path = '/workspace/nr-link-simulator/scripts/sionna_tdl_crosscheck.py'
os.makedirs(results_dir, exist_ok=True)

DELAY_SPREAD = 100
DOPPLER = 10
MAX_BLOCKS = 100
TARGET_ERRORS = 30
BATCH_SIZE = 2

configs = [
    # QPSK MCS4
    {'mcs': 4, 'ce': 'ideal', 'snr_start': -2, 'snr_end': 10, 'name': 'sionna_ideal_mcs4'},
    {'mcs': 4, 'ce': 'ls', 'snr_start': 4, 'snr_end': 20, 'name': 'sionna_ls_mcs4'},
    # 16QAM MCS8
    {'mcs': 8, 'ce': 'ideal', 'snr_start': 6, 'snr_end': 20, 'name': 'sionna_ideal_mcs8'},
    {'mcs': 8, 'ce': 'ls', 'snr_start': 12, 'snr_end': 26, 'name': 'sionna_ls_mcs8'},
    # 64QAM MCS15
    {'mcs': 15, 'ce': 'ideal', 'snr_start': 14, 'snr_end': 28, 'name': 'sionna_ideal_mcs15'},
    {'mcs': 15, 'ce': 'ls', 'snr_start': 20, 'snr_end': 35, 'name': 'sionna_ls_mcs15'},
    # 256QAM MCS27
    {'mcs': 27, 'ce': 'ideal', 'snr_start': 22, 'snr_end': 35, 'name': 'sionna_ideal_mcs27'},
    {'mcs': 27, 'ce': 'ls', 'snr_start': 28, 'snr_end': 42, 'name': 'sionna_ls_mcs27'},
]

for cfg in configs:
    out_file = os.path.join(results_dir, f"{cfg['name']}.csv")
    
    if os.path.exists(out_file) and os.path.getsize(out_file) > 100:
        print(f"\nSkipping {cfg['name']} (output exists)")
        continue
    
    print(f"\n{'='*60}")
    print(f"Running {cfg['name']}: MCS={cfg['mcs']}, CE={cfg['ce']}, SNR {cfg['snr_start']}-{cfg['snr_end']} dB")
    print(f"Delay spread={DELAY_SPREAD}ns, Doppler={DOPPLER}Hz")
    print(f"{'='*60}")
    
    cmd = [
        'python', script_path,
        '--mcs', str(cfg['mcs']),
        '--snr-start', str(cfg['snr_start']),
        '--snr-end', str(cfg['snr_end']),
        '--snr-step', '1',
        '--ce-type', cfg['ce'],
        '--delay-spread', str(DELAY_SPREAD),
        '--doppler', str(DOPPLER),
        '--max-blocks', str(MAX_BLOCKS),
        '--target-errors', str(TARGET_ERRORS),
        '--batch-size', str(BATCH_SIZE),
        '--output', out_file
    ]
    
    try:
        subprocess.run(cmd, check=True, cwd='/workspace/nr-link-simulator')
        print(f"Results saved to {out_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {cfg['name']}: {e}")

print("\n" + "="*60)
print("All Sionna simulations completed!")
print("="*60)

import glob
for f in sorted(glob.glob(os.path.join(results_dir, 'sionna_*.csv'))):
    print(f"  {os.path.basename(f)}: {os.path.getsize(f)} bytes")
