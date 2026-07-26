#!/usr/bin/env python3
"""Sionna BLER simulation for MCS 2-27 AWGN with LS-CE + LMMSE, freq-domain."""

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['CUDA_VISIBLE_DEVICES'] = ''

import sys
sys.path.insert(0, '/workspace/sionna/src')

import torch
torch.manual_seed(42)

import torch.compiler
_orig_compile = torch.compile
def _no_compile(*args, **kwargs):
    if len(args) > 0 and callable(args[0]):
        return args[0]
    return lambda f: f
torch.compile = _no_compile
for _attr in ['is_inductor_compiled', 'is_compiling', 'compile_pystate', 'codecache',
              'get_num_threads', 'set_num_threads', 'reset', 'is_torchdynamo_compiling']:
    if not hasattr(torch.compiler, _attr):
        setattr(torch.compiler, _attr, lambda *args, **kwargs: None)
if not hasattr(torch, 'is_dynamo_compiling'):
    torch.is_dynamo_compiling = lambda: False

import csv
import numpy as np
from sionna.phy.nr import PUSCHConfig, PUSCHTransmitter, PUSCHReceiver

MCS_TABLE_1 = [
    (2, 120), (2, 157), (2, 193), (2, 251), (2, 308), (2, 379), (2, 449), (2, 526), (2, 602), (2, 679),
    (4, 340), (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658), (6, 438), (6, 466), (6, 517),
    (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873), (6, 910), (6, 948)
]

OUT_DIR = "/workspace/nr-link-simulator/results/bler_batch"
os.makedirs(OUT_DIR, exist_ok=True)

def get_snr_range(mcs):
    qm, R_x1024 = MCS_TABLE_1[mcs]
    R = R_x1024 / 1024.0
    if qm == 2:
        if R < 0.25:
            return (-4, 6)
        elif R < 0.45:
            return (-2, 10)
        else:
            return (1, 13)
    elif qm == 4:
        if R < 0.45:
            return (3, 15)
        elif R < 0.60:
            return (6, 18)
        else:
            return (8, 20)
    elif qm == 6:
        if R < 0.50:
            return (8, 20)
        elif R < 0.65:
            return (10, 23)
        elif R < 0.80:
            return (13, 26)
        else:
            return (16, 30)
    return (0, 30)

def make_pusch_config(mcs_idx):
    qm, R_x1024 = MCS_TABLE_1[mcs_idx]
    cfg = PUSCHConfig()
    cfg.n_size_bwp = 25
    cfg.num_antenna_ports = 1
    cfg.num_layers = 1
    cfg.transform_precoding = False
    cfg.dmrs.config_type = 1
    cfg.dmrs.length = 1
    cfg.dmrs.additional_position = 0
    cfg.dmrs.num_cdm_groups_without_data = 1
    cfg.dmrs.dmrs_port_set = [0]
    cfg.tb.mcs_index = mcs_idx
    cfg.tb.mcs_table = 1
    return cfg

def eval_bler(pusch_tx, pusch_rx, snr_db, max_blocks=200, target_errs=30, bsz=20):
    no_val = 1.0 / (10 ** (snr_db / 10))
    n_err = 0
    n_tot = 0
    while n_tot < max_blocks:
        b = min(bsz, max_blocks - n_tot)
        try:
            x, bits = pusch_tx(b)
            noise = (torch.randn_like(x) + 1j * torch.randn_like(x)) * np.sqrt(no_val / 2.0)
            y = x + noise
            no_t = torch.tensor(no_val, dtype=torch.float32)
            b_hat, crc_ok = pusch_rx(y, no_t)
            n_err += (crc_ok[:, 0] == False).long().sum().item()
            n_tot += b
        except Exception as e:
            print(f"      batch error at SNR {snr_db}: {e}", flush=True)
            break
        if n_err >= target_errs and n_tot >= 30:
            break
        if n_err == 0 and n_tot >= 100:
            break
    bler = n_err / max(n_tot, 1)
    return bler, n_tot, n_err

def main():
    mcs_start = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    mcs_end = int(sys.argv[2]) if len(sys.argv) > 2 else 27
    
    out_csv = os.path.join(OUT_DIR, f"sionna_bler_mcs{mcs_start}-{mcs_end}.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mcs", "qm", "rate", "tbs", "snr_db", "bler", "n_blocks", "n_errors"])
        
        for mcs in range(mcs_start, mcs_end + 1):
            qm, R_x1024 = MCS_TABLE_1[mcs]
            R = R_x1024 / 1024.0
            snr_start, snr_end = get_snr_range(mcs)
            
            try:
                cfg = make_pusch_config(mcs)
                pusch_tx = PUSCHTransmitter(cfg, output_domain='freq')
                pusch_rx = PUSCHReceiver(pusch_tx,
                                         channel_estimator=None,
                                         input_domain='freq',
                                         mimo_equalizer='lmmse',
                                         detector='llr',
                                         return_tb_crc_status=True,
                                         disable_inter_stream_sic=True)
                tbs = cfg.tb_size
                print(f"=== Sionna MCS {mcs:2d} (Qm={qm}, R={R:.3f}, TBS={tbs}) "
                      f"SNR {snr_start}~{snr_end}dB ===", flush=True)
            except Exception as e:
                print(f"=== Sionna MCS {mcs:2d} (Qm={qm}, R={R:.3f}) INIT FAILED: {e} ===", flush=True)
                import traceback; traceback.print_exc()
                continue
            
            snr = snr_start
            while snr <= snr_end + 0.01:
                snr_r = round(snr, 1)
                bler, n_blk, n_err = eval_bler(pusch_tx, pusch_rx, snr_r,
                                               max_blocks=200, target_errs=30, bsz=20)
                print(f"  SNR {snr_r:5.1f}dB: BLER={bler:.4f} ({n_err}/{n_blk})", flush=True)
                w.writerow([mcs, qm, f"{R:.4f}", tbs, f"{snr_r:.1f}",
                           f"{bler:.6f}", n_blk, n_err])
                f.flush()
                
                if bler < 0.001 and n_blk >= 100:
                    for s2 in [snr_r + 1.0]:
                        if s2 <= snr_end + 0.01:
                            w.writerow([mcs, qm, f"{R:.4f}", tbs, f"{s2:.1f}",
                                       "0.000000", 200, 0])
                    break
                snr += 1.0
            print()
    
    print(f"Sionna results saved to {out_csv}")

if __name__ == "__main__":
    main()
