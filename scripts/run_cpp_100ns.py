import subprocess
import os

results_dir = '/workspace/nr-link-simulator/results/ce_compare_100ns'
cpp_bin = '/workspace/nr-link-simulator/build/examples/pdsch_bler_tdl_mimo'
os.makedirs(results_dir, exist_ok=True)

# Common args: mcs, sinr_start, sinr_end, step, perfect_csi, max_blocks, target_errors,
#             n_rb, n_rx, n_layers, channel_type, delay_spread_ns, max_doppler,
#             dmrs_add_pos, dmrs_duration, mcs_table, outfile
N_RB = 30
N_RX = 4
N_LAYERS = 4
CHANNEL_TYPE = 1
DELAY_SPREAD = 100
MAX_DOPPLER = 10
DMRS_ADD_POS = 1
DMRS_DURATION = 0
MCS_TABLE = 2

cpp_configs = [
    # QPSK MCS4
    {'mcs': 4, 'ce': 'ideal', 'snr_start': -2, 'snr_end': 10, 'name': 'cpp_ideal_mcs4'},
    {'mcs': 4, 'ce': 'ls', 'snr_start': 4, 'snr_end': 22, 'name': 'cpp_ls_mcs4'},
    # 16QAM MCS8
    {'mcs': 8, 'ce': 'ideal', 'snr_start': 6, 'snr_end': 22, 'name': 'cpp_ideal_mcs8'},
    {'mcs': 8, 'ce': 'ls', 'snr_start': 12, 'snr_end': 28, 'name': 'cpp_ls_mcs8'},
    # 64QAM MCS15
    {'mcs': 15, 'ce': 'ideal', 'snr_start': 14, 'snr_end': 30, 'name': 'cpp_ideal_mcs15'},
    {'mcs': 15, 'ce': 'ls', 'snr_start': 20, 'snr_end': 38, 'name': 'cpp_ls_mcs15'},
    # 256QAM MCS27
    {'mcs': 27, 'ce': 'ideal', 'snr_start': 22, 'snr_end': 38, 'name': 'cpp_ideal_mcs27'},
    {'mcs': 27, 'ce': 'ls', 'snr_start': 28, 'snr_end': 45, 'name': 'cpp_ls_mcs27'},
]

for cfg in cpp_configs:
    out_file = os.path.join(results_dir, f"{cfg['name']}.csv")
    
    if os.path.exists(out_file) and os.path.getsize(out_file) > 100:
        print(f"\nSkipping {cfg['name']} (output exists)")
        continue
    
    perfect_csi = 1 if cfg['ce'] == 'ideal' else 0
    max_blocks = 100
    target_errors = 30
    snr_step = 1
    
    print(f"\n{'='*60}")
    print(f"Running C++ {cfg['name']}: MCS={cfg['mcs']}, CE={cfg['ce']}, SNR {cfg['snr_start']}-{cfg['snr_end']} dB")
    print(f"Delay spread={DELAY_SPREAD}ns, Doppler={MAX_DOPPLER}Hz")
    print(f"{'='*60}")
    
    cmd = [
        cpp_bin,
        str(cfg['mcs']),
        str(cfg['snr_start']),
        str(cfg['snr_end']),
        str(snr_step),
        str(perfect_csi),
        str(max_blocks),
        str(target_errors),
        str(N_RB),
        str(N_RX),
        str(N_LAYERS),
        str(CHANNEL_TYPE),
        str(DELAY_SPREAD),
        str(MAX_DOPPLER),
        str(DMRS_ADD_POS),
        str(DMRS_DURATION),
        str(MCS_TABLE),
        out_file
    ]
    
    try:
        subprocess.run(cmd, check=True, cwd='/workspace/nr-link-simulator')
        print(f"Results saved to {out_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {cfg['name']}: {e}")

print("\n" + "="*60)
print("All C++ simulations completed!")
print("="*60)

import glob
for f in sorted(glob.glob(os.path.join(results_dir, 'cpp_*.csv'))):
    print(f"  {os.path.basename(f)}: {os.path.getsize(f)} bytes")
