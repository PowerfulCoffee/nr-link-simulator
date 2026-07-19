import os
import sys
import csv
import math
import numpy as np
import torch

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import sionna
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Mapper, Demapper, Constellation
from sionna.phy.utils import BinarySource

n_prb = 3
n_re_per_prb = 156
n_layers = 1
num_ldpc_iters = 25
max_blocks = 5000
target_errors = 100
batch_size = 200

MCS_TABLE = [
    (2, 120),   # MCS0
    (2, 193),   # MCS1
    (2, 308),   # MCS2
    (2, 449),   # MCS3
    (2, 602),   # MCS4
    (4, 378),   # MCS5
    (4, 434),   # MCS6
    (4, 490),   # MCS7
    (4, 553),   # MCS8
    (4, 616),   # MCS9
    (4, 658),   # MCS10
    (6, 466),   # MCS11
    (6, 517),   # MCS12
    (6, 567),   # MCS13
    (6, 616),   # MCS14
    (6, 666),   # MCS15
    (6, 719),   # MCS16
    (6, 772),   # MCS17
    (6, 822),   # MCS18
    (6, 873),   # MCS19
    (8, 682.5), # MCS20
    (8, 711),   # MCS21
    (8, 754),   # MCS22
    (8, 797),   # MCS23
    (8, 841),   # MCS24
    (8, 885),   # MCS25
    (8, 916.5), # MCS26
    (8, 948),   # MCS27
]

CPP_TBS = {
    0:104, 1:176, 2:288, 3:408, 4:552,
    5:704, 6:808, 7:888, 8:1032, 9:1128, 10:1224,
    11:1288, 12:1416, 13:1608, 14:1736, 15:1864,
    16:2024, 17:2152, 18:2280, 19:2408,
    20:2472, 21:2600, 22:2792, 23:2976, 24:3104, 25:3240, 26:3368, 27:3496
}

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
    qm, r1024 = MCS_TABLE[mcs]
    R = r1024 / 1024.0
    tbs_cpp = CPP_TBS[mcs]

    try:
        tbs_info = calculate_tb_size(
            n_prb, n_re_per_prb, qm, n_layers, r1024,
            tb_crc_length=24, cb_crc_length=24, verbose=False
        )
        tbs = int(tbs_info['tb_size'].numpy()) if torch.is_tensor(tbs_info['tb_size']) else int(tbs_info['tb_size'])
        n_cbs = int(tbs_info['n_cbs'])
    except Exception as e:
        tbs = tbs_cpp
        n_cbs = 0

    print(f"  MCS{mcs}: Qm={qm}, R={R:.3f}, TBS_cpp={tbs_cpp}, TBS_sionna={tbs}", flush=True)
    use_tbs = tbs_cpp

    try:
        tb_enc = TBEncoder(
            n_prb, n_re_per_prb, qm, n_layers, r1024,
            num_bits_per_cw=None, tb_crc_length=24, cb_crc_length=24, verbose=False
        )
        tb_dec = TBDecoder(tb_enc, num_iter=num_ldpc_iters, hard_out=True, verbose=False)
    except Exception as e:
        err_str = str(e)
        if "rate" in err_str.lower() or "1/5" in err_str or "0.2" in err_str:
            print(f"  Skipping MCS{mcs}: Sionna does not support R<1/5 (low rate requires repetition coding)", flush=True)
        else:
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
            current_bs = min(batch_size, max_blocks - n_blocks)
            b = source([current_bs, use_tbs])

            c = tb_enc(b)
            x = mapper(c)

            noise_std = math.sqrt(no / 2.0)
            noise_r = torch.randn_like(x) * noise_std
            noise_i = torch.randn_like(x) * noise_std
            y = x + torch.complex(noise_r, noise_i)

            llr = demapper([y, no])
            b_hat, tb_crc_ok = tb_dec(llr)

            b_hat_hard = (b_hat > 0).to(torch.float32)
            b_target = b.to(torch.float32)
            if b_hat_hard.shape != b_target.shape:
                min_len = min(b_hat_hard.shape[1], b_target.shape[1])
                b_hat_hard = b_hat_hard[:, :min_len]
                b_target = b_target[:, :min_len]
            bit_errors = torch.sum(torch.abs(b_hat_hard - b_target), dim=1)
            block_errors = torch.sum((bit_errors > 0).to(torch.float32)).item()
            n_errors += int(block_errors)
            n_blocks += current_bs

        bler = n_errors / n_blocks if n_blocks > 0 else 1.0
        results.append((esn0_db, n_blocks, n_errors, bler))
        print(f"    SNR {esn0_db:5.1f} dB: BLER={bler:.4f} ({n_blocks} blocks, {n_errors} errors)", flush=True)

        if bler < 1e-4 and n_blocks >= 500:
            break

    return results

def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
    os.makedirs(out_dir, exist_ok=True)

    out_csv = os.path.join(out_dir, 'sionna_all_mcs_matched_awgn.csv')
    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['mcs','qm','rate','tbs','snr_db','bler','n_blocks','n_errors'])

    print(f"Sionna AWGN BLER - All MCS, nPRB={n_prb}, n_re_per_prb={n_re_per_prb}, SISO", flush=True)
    print(f"LDPC iters={num_ldpc_iters}, max_blocks={max_blocks}, target_errs={target_errors}", flush=True)
    print()

    for mcs in range(28):
        qm, r1024 = MCS_TABLE[mcs]
        R = r1024 / 1024.0
        tbs_cpp = CPP_TBS[mcs]
        print(f"=== MCS {mcs:2d} (Qm={qm}, R={R:.3f}, TBS={tbs_cpp}) ===", flush=True)
        results = run_mcs_bler(mcs)
        with open(out_csv, 'a', newline='') as f:
            writer = csv.writer(f)
            for snr, blocks, errs, bler in results:
                writer.writerow([mcs, qm, f"{R:.4f}", tbs_cpp, f"{snr:.1f}", f"{bler:.6f}", blocks, errs])
        f.flush = os.fsync(f.fileno()) if hasattr(f, 'fileno') else None
        print(flush=True)

    print(f"\nResults saved to {out_csv}", flush=True)

if __name__ == '__main__':
    main()
