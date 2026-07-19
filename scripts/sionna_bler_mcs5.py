#!/usr/bin/env python3
import sys
sys.path.insert(0, '/workspace/sionna/src')
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import time
import torch
import numpy as np
from sionna.phy.nr import PUSCHConfig, PUSCHTransmitter, PUSCHReceiver
from sionna.phy.channel import AWGN

torch.manual_seed(42)

pusch_config = PUSCHConfig()
pusch_config.n_size_bwp = 25
pusch_config.num_antenna_ports = 1
pusch_config.num_layers = 1
pusch_config.transform_precoding = False
pusch_config.dmrs.config_type = 1
pusch_config.dmrs.length = 1
pusch_config.dmrs.additional_position = 0
pusch_config.dmrs.dmrs_port_set = [0]
pusch_config.tb.mcs_index = 5
pusch_config.tb.mcs_table = 1

print(f"MCS 5, Qm={pusch_config.tb.num_bits_per_symbol.item()}, R={pusch_config.tb.target_coderate.item():.4f}", flush=True)
print(f"TBS={pusch_config.tb_size}, N_coded={pusch_config.num_coded_bits}", flush=True)

pusch_tx = PUSCHTransmitter(pusch_config, output_domain='freq')
pusch_rx = PUSCHReceiver(pusch_tx, channel_estimator=None, mimo_equalizer='lmmse', detector='llr', return_tb_crc_status=True)
channel = AWGN()

def eval_bler(snr_db, max_blocks=1000, target_errs=50, bsz=50):
    no_val = 1.0 / (10 ** (snr_db / 10.0))
    n_err = 0
    n_tot = 0
    t0 = time.time()
    while n_tot < max_blocks and n_err < target_errs:
        b = min(bsz, max_blocks - n_tot)
        x = pusch_tx(b)
        if isinstance(x, tuple):
            x = x[0]
        y = channel(x, torch.tensor(no_val, dtype=torch.float32))
        res = pusch_rx(y, torch.tensor(no_val, dtype=torch.float32))
        crc = res[1]
        n_err += (crc[0] == False).long().sum().item()
        n_tot += b
    dt = time.time() - t0
    return n_err / n_tot, n_tot, n_err, dt

print(f"\n=== Sionna LS CE (freq domain, Es/N0) ===", flush=True)
print(f"{'SNR(dB)':>8} {'BLER':>8} {'Blocks':>8} {'Errors':>8} {'Time(s)':>8}", flush=True)
print("-" * 44, flush=True)
for snr in [4, 3, 2, 1, 0, -1, -2]:
    bler, nblk, nerr, dt = eval_bler(snr, max_blocks=200, target_errs=30, bsz=20)
    print(f"{snr:>8.1f} {bler:>8.4f} {nblk:>8d} {nerr:>8d} {dt:>8.1f}", flush=True)
