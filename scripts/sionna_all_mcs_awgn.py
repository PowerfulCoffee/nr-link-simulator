#!/usr/bin/env python3
"""
Sionna reference: PDSCH BLER vs SNR for MCS 0-27 over AWGN channel.
Uses Sionna's LDPC5GEncoder/LDPC5GDecoder + Mapper/Demapper + AWGN directly
for exact G-matching with nr-link-simulator.
"""

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['CUDA_VISIBLE_DEVICES'] = ''

import sys
import torch
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

import numpy as np
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.mapping import Constellation, Mapper, Demapper
from sionna.phy.channel import AWGN
from sionna.phy.nr.utils import calculate_tb_size
import csv
import torch
torch.manual_seed(42)

def binary_source(shape):
    return torch.randint(0, 2, shape, dtype=torch.float32)

MCS_TABLE_PDSCH2 = [
    (2, 120),   (2, 193),   (2, 308),   (2, 449),   (2, 602),
    (4, 378),   (4, 434),   (4, 490),   (4, 553),   (4, 616),
    (4, 658),   (6, 466),   (6, 517),   (6, 567),   (6, 616),
    (6, 666),   (6, 719),   (6, 772),   (6, 822),   (6, 873),
    (8, 682.5), (8, 711),   (8, 754),   (8, 797),   (8, 841),
    (8, 885),   (8, 916.5), (8, 948),
]

N_PRB = 3
N_RE_PER_PRB = 156
N_LAYERS = 1
NUM_ITERS = 25
BATCH_SIZE = 200
MAX_BLOCKS = 5000
TARGET_ERRORS = 100
SNR_STEP = 1.0


def get_snr_range(qm, R):
    if qm == 2:
        if R < 0.20:   return (-10.0, 6.0)
        elif R < 0.40: return (-6.0, 10.0)
        else:          return (-2.0, 14.0)
    elif qm == 4:
        if R < 0.40:   return (0.0, 16.0)
        else:          return (4.0, 20.0)
    elif qm == 6:
        if R < 0.55:   return (6.0, 22.0)
        elif R < 0.75: return (10.0, 26.0)
        else:          return (14.0, 30.0)
    else:
        if R < 0.70:   return (12.0, 30.0)
        elif R < 0.85: return (16.0, 32.0)
        else:          return (20.0, 34.0)


def get_snr_points(snr_start, snr_end, step=1.0):
    pts = []
    s = snr_start
    while s <= snr_end + 0.01:
        pts.append(round(s, 1))
        s += step
    return pts


def select_bg(k_crc, G):
    R = k_crc / G
    k_info = k_crc
    if k_info <= 292:
        return "bg2"
    elif k_info <= 3824 and R <= 0.67:
        return "bg2"
    elif R <= 0.25:
        return "bg2"
    else:
        return "bg1"


class PdschAwgnSim:
    def __init__(self, k_crc, G, qm, bg=None, num_iters=25, device='cpu'):
        self.k = k_crc
        self.n = G
        self.qm = qm
        self.device = device

        self.encoder = LDPC5GEncoder(
            k_crc, G,
            num_bits_per_symbol=qm,
            bg=bg,
            device=device,
        )
        self.decoder = LDPC5GDecoder(
            self.encoder,
            num_iter=num_iters,
            hard_out=True,
            cn_update="boxplus-phi",
            device=device,
        )
        self.constellation = Constellation("qam", qm, device=device)
        self.mapper = Mapper(constellation=self.constellation, device=device)
        self.demapper = Demapper("app", constellation=self.constellation, device=device)
        self.awgn = AWGN(device=device)

    def run_batch(self, batch_size, esn0_db):
        bits = binary_source([batch_size, self.k])
        cw = self.encoder(bits)
        x = self.mapper(cw)

        no = 10 ** (-esn0_db / 10.0)

        y = self.awgn(x, no)
        llr = self.demapper(y, no)
        bits_hat = self.decoder(llr)

        errors = (bits != bits_hat).any(dim=1).sum().item()
        return int(errors)


def main():
    out_dir = "/workspace/nr-link-simulator/build/results"
    os.makedirs(out_dir, exist_ok=True)
    print(f"Sionna PDSCH BLER simulation (AWGN, nPRB={N_PRB}, {N_RE_PER_PRB} RE/PRB)")
    print(f"Step={SNR_STEP}dB, max_blocks={MAX_BLOCKS}, target_errors={TARGET_ERRORS}, "
          f"decoder_iters={NUM_ITERS}")
    print(f"MCS Table: PDSCH Table 2 (TS 38.214, 256QAM MCS20-27)\n")

    out_file = os.path.join(out_dir, "bler_all_mcs_sionna.csv")
    with open(out_file, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mcs", "qm", "rate", "tbs", "snr_db", "bler", "n_blocks", "n_errors"])

        for mcs_idx in range(28):
            qm, R_x1024 = MCS_TABLE_PDSCH2[mcs_idx]
            R = R_x1024 / 1024.0
            G = N_PRB * N_RE_PER_PRB * qm * N_LAYERS

            tb_size, cb_size, num_cb, tb_crc, cb_crc, cw_len = calculate_tb_size(
                modulation_order=qm, target_coderate=R,
                num_coded_bits=G, num_layers=N_LAYERS, verbose=False)

            tbs = int(tb_size.numpy())
            k_crc = int(cb_size.numpy())
            C = int(num_cb.numpy())
            bg = None
            eff_rate = k_crc / G

            snr_start, snr_end = get_snr_range(qm, R)
            snr_points = get_snr_points(snr_start, snr_end, SNR_STEP)

            print(f"=== MCS {mcs_idx:2d} (Qm={qm}, R={R:.3f}, TBS={tbs}, k_crc={k_crc}, "
                  f"G={G}, C={C}, eff_rate={eff_rate:.3f}) SNR {snr_start:.0f}~{snr_end:.0f}dB ===")

            skip = eff_rate < 0.20
            sim = None
            if not skip:
                try:
                    sim = PdschAwgnSim(k_crc, G, qm, bg=bg, num_iters=NUM_ITERS)
                except Exception as e:
                    print(f"  Init failed: {e}")
                    skip = True

            for snr_db in snr_points:
                if skip:
                    w.writerow([mcs_idx, qm, f"{R:.4f}", tbs, f"{snr_db:.1f}",
                                "0.0", 0, 0])
                    continue

                n_errors = 0
                n_blocks = 0
                while n_blocks < MAX_BLOCKS:
                    bs = min(BATCH_SIZE, MAX_BLOCKS - n_blocks)
                    try:
                        n_errors += sim.run_batch(bs, snr_db)
                    except Exception as e:
                        print(f"    SNR {snr_db}dB batch failed: {e}")
                        break
                    n_blocks += bs
                    if n_errors >= TARGET_ERRORS and n_blocks >= 10:
                        break
                    if n_errors == 0 and n_blocks >= 500:
                        break

                bler = n_errors / n_blocks if n_blocks > 0 else 0.0
                print(f"  SNR {snr_db:5.1f} dB: BLER={bler:.4f} "
                      f"(blocks={n_blocks}, errors={n_errors})")
                w.writerow([mcs_idx, qm, f"{R:.4f}", tbs, f"{snr_db:.1f}",
                           f"{bler:.6f}", n_blocks, n_errors])
            f.flush()
            print()

    print(f"Results saved to {out_file}")


if __name__ == "__main__":
    main()
