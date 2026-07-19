import os
import sys
import csv
import math
import numpy as np
import torch

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import sionna
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Mapper, Demapper, Constellation
from sionna.phy.utils import BinarySource

device = torch.device('cpu')

n_prb = 3
n_re_per_prb = 156
n_layers = 1
num_ldpc_iters = 25
max_blocks = 5000
target_errors = 100
batch_size = 100

MCS_TABLE = []
for i in range(28):
    if i <= 1:
        qm, r = 2, [120, 157][i]
    elif i <= 4:
        qm, r = 2, [193, 251, 308][i-2]
    elif i <= 6:
        qm, r = 2, [379, 449, 526][i-5]
    elif i <= 9:
        qm, r = 4, [308, 340, 378, 434, 490][i-10] if i<=10 else None
    elif i <= 10:
        qm, r = 4, [378, 434, 490, 553, 616, 658][i-5]

MCS_TABLE_T2 = [
    (2, 120),   # MCS0
    (2, 157),   # MCS1
    (2, 193),   # MCS2
    (2, 251),   # MCS3
    (2, 308),   # MCS4
    (2, 379),   # MCS5? wait - no, MCS5 starts 16QAM
]

def mcs_params_table2(mcs):
    table = [
        (2, 120), (2, 157), (2, 193), (2, 251), (2, 308), (2, 379), (2, 449), (2, 526),
        (2, 602), (2, 658),
        (4, 340), (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658),
        (6, 466), (6, 490), (6, 517), (6, 567), (6, 616), (6, 666),
        (8, 682), (8, 711), (8, 754), (8, 797), (8, 841), (8, 885), (8, 910),
        (8, 948)
    ]
    qm, r1024 = table[mcs]
    R = r1024 / 1024.0
    return qm, R, r1024

def get_snr_range(mcs, qm, R):
    if qm == 2:
        if R < 0.20:
            return -10.0, 6.0
        elif R < 0.45:
            return -6.0, 10.0
        else:
            return -4.0, 12.0
    elif qm == 4:
        if R < 0.50:
            return 0.0, 16.0
        else:
            return 2.0, 18.0
    elif qm == 6:
        if R < 0.60:
            return 4.0, 22.0
        elif R < 0.76:
            return 8.0, 24.0
        else:
            return 10.0, 26.0
    else:
        if R < 0.75:
            return 12.0, 28.0
        elif R < 0.87:
            return 16.0, 32.0
        else:
            return 18.0, 34.0

def run_mcs_bler(mcs):
    qm, R, r1024 = mcs_params_table2(mcs)
    n_coded_bits = n_prb * n_re_per_prb * qm * n_layers

    tb_size = calculate_tb_size(
        n_prb,
        n_re_per_prb,
        qm,
        n_layers,
        r1024,
        tb_crc_length=24,
        cb_crc_length=24,
        verbose=False
    )
    k = tb_size['tb_size']
    n_cbs = tb_size['n_cbs']
    cw_lengths = tb_size['cw_lengths']
    tbs = k

    print(f"  MCS{mcs}: Qm={qm}, R={R:.3f}, TBS={tbs}, E={n_coded_bits}, n_cbs={n_cbs}", flush=True)

    try:
        tb_enc = TBEncoder(
            n_prb,
            n_re_per_prb,
            qm,
            n_layers,
            r1024,
            num_bits_per_cw=None,
            tb_crc_length=24,
            cb_crc_length=24,
            verbose=False
        )
        tb_dec = TBDecoder(tb_enc, num_iter=num_ldpc_iters, hard_out=True, verbose=False)
    except Exception as e:
        print(f"  Skipping MCS{mcs}: {e}", flush=True)
        return []

    constellation = Constellation("qam", qm)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)
    source = BinarySource()

    snr_start, snr_end = get_snr_range(mcs, qm, R)
    snr_range = np.arange(snr_start, snr_end + 0.5, 1.0)
    results = []

    for esn0_db in snr_range:
        esn0_lin = 10 ** (esn0_db / 10.0)
        no = 1.0 / esn0_lin

        n_errors = 0
        n_blocks = 0

        while n_blocks < max_blocks and n_errors < target_errors:
            bs = min(batch_size, max_blocks - n_blocks)
            b = source([bs, tbs])

            c = tb_enc(b)
            x = mapper(c)

            noise_std = math.sqrt(no / 2.0)
            noise_r = torch.randn_like(x) * noise_std
            noise_i = torch.randn_like(x) * noise_std
            y = x + torch.complex(noise_r, noise_i)

            llr = demapper([y, no])
            b_hat, tb_crc_ok = tb_dec(llr)

            b_hat_hard = (b_hat > 0).to(torch.float32)
            bit_errors = torch.sum(torch.abs(b_hat_hard - b), dim=1)
            block_errors = torch.sum((bit_errors > 0).to(torch.float32)).item()
            n_errors += int(block_errors)
            n_blocks += bs

        bler = n_errors / n_blocks if n_blocks > 0 else 1.0
        results.append((esn0_db, n_blocks, n_errors, bler))
        print(f"    SNR {esn0_db:5.1f} dB: BLER={bler:.4f} ({n_blocks} blocks, {n_errors} errors)", flush=True)

        if bler < 1e-4 and n_blocks >= 500:
            break

    return results

def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
    os.makedirs(out_dir, exist_ok=True)

    out_csv = os.path.join(out_dir, 'sionna_all_mcs_table2_awgn.csv')
    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['mcs','qm','rate','tbs','snr_db','bler','n_blocks','n_errors'])

    print(f"Sionna AWGN BLER - All MCS Table 2 (256QAM), nPRB={n_prb}, SISO", flush=True)
    print(f"LDPC iters={num_ldpc_iters}, max_blocks={max_blocks}, target_errs={target_errors}", flush=True)
    print()

    for mcs in range(28):
        qm, R, _ = mcs_params_table2(mcs)
        print(f"=== MCS {mcs:2d} (Qm={qm}, R={R:.3f}) ===", flush=True)
        results = run_mcs_bler(mcs)
        with open(out_csv, 'a', newline='') as f:
            writer = csv.writer(f)
            for snr, blocks, errs, bler in results:
                tbs = 0
                writer.writerow([mcs, qm, f"{R:.4f}", tbs, f"{snr:.1f}", f"{bler:.6f}", blocks, errs])
        print(flush=True)

    print(f"\nResults saved to {out_csv}", flush=True)

if __name__ == '__main__':
    main()
