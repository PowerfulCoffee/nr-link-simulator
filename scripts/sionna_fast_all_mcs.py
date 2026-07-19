import os
import csv
import math
import numpy as np
import torch

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from sionna.phy.nr import TBEncoder, TBDecoder
from sionna.phy.mapping import Mapper, Demapper, Constellation

n_prb = 3
n_re_per_prb = 156
n_layers = 1
num_ldpc_iters = 25
max_blocks = 2000
target_errors = 50
batch_size = 400

MCS_TABLE = [
    (2, 120), (2, 193), (2, 308), (2, 449), (2, 602),
    (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658),
    (6, 466), (6, 517), (6, 567), (6, 616), (6, 666),
    (6, 719), (6, 772), (6, 822), (6, 873),
    (8, 682.5), (8, 711), (8, 754), (8, 797), (8, 841), (8, 885), (8, 916.5), (8, 948),
]

CPP_TBS = {
    0:104, 1:176, 2:288, 3:408, 4:552,
    5:704, 6:808, 7:888, 8:1032, 9:1128, 10:1224,
    11:1288, 12:1416, 13:1608, 14:1736, 15:1864,
    16:2024, 17:2152, 18:2280, 19:2408,
    20:2472, 21:2600, 22:2792, 23:2976, 24:3104, 25:3240, 26:3368, 27:3496
}

def get_snr_range(qm, R):
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
    tbs = CPP_TBS[mcs]
    n_coded_bits = n_prb * n_re_per_prb * qm * n_layers
    eff_rate = tbs / n_coded_bits

    print(f"  Qm={qm}, R={R:.3f}, TBS={tbs}, n_coded_bits={n_coded_bits}, eff_rate={eff_rate:.3f}", flush=True)

    if eff_rate < 1.0/5.0:
        print(f"  Skipping (R={eff_rate:.3f} < 0.2: Sionna doesn't support repetition)", flush=True)
        return []

    tb_enc = TBEncoder(
        target_tb_size=tbs, num_coded_bits=n_coded_bits, target_coderate=R,
        num_bits_per_symbol=qm, num_layers=n_layers, channel_type='PUSCH',
        use_scrambler=True, verbose=False
    )
    tb_dec = TBDecoder(tb_enc, num_iter=num_ldpc_iters, hard_out=True, verbose=False)
    constellation = Constellation("qam", qm)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)

    snr_start, snr_end = get_snr_range(qm, R)
    snr_range = np.arange(snr_start, snr_end + 0.5, 2.0)
    results = []

    for esn0_db in snr_range:
        no_val = 1.0 / (10 ** (esn0_db / 10.0))
        no_t = torch.tensor(no_val, dtype=torch.float32)
        n_errors = 0
        n_blocks = 0

        while n_blocks < max_blocks and n_errors < target_errors:
            current_bs = min(batch_size, max_blocks - n_blocks)
            b = torch.randint(0, 2, (current_bs, tbs), dtype=torch.float32)
            c = tb_enc(b)
            x = mapper(c)
            noise_std = math.sqrt(no_val / 2.0)
            noise = torch.complex(
                torch.randn_like(x, dtype=torch.float32) * noise_std,
                torch.randn_like(x, dtype=torch.float32) * noise_std
            ).to(x.dtype)
            y = x + noise
            llr = demapper(y, no_t)
            _, crc_ok = tb_dec(llr)
            block_ok = crc_ok.to(torch.float32)
            block_errors = current_bs - int(torch.sum(block_ok).item())
            n_errors += block_errors
            n_blocks += current_bs

        bler = n_errors / n_blocks if n_blocks > 0 else 1.0
        results.append((esn0_db, n_blocks, n_errors, bler))
        print(f"    SNR {esn0_db:5.1f} dB: BLER={bler:.4f} ({n_blocks} blocks, {n_errors} errors)", flush=True)

        if bler < 1e-4 and n_blocks >= 400:
            break

    return results

def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, 'sionna_all_mcs_awgn_fast.csv')
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['mcs','qm','rate','tbs','snr_db','bler','n_blocks','n_errors'])

    print(f"Sionna AWGN BLER (fast) - {max_blocks} max blocks, {target_errors} target errors", flush=True)
    for mcs in range(28):
        qm, r1024 = MCS_TABLE[mcs]
        R = r1024 / 1024.0
        tbs = CPP_TBS[mcs]
        print(f"=== MCS {mcs:2d} (Qm={qm}, R={R:.3f}, TBS={tbs}) ===", flush=True)
        results = run_mcs_bler(mcs)
        with open(out_csv, 'a', newline='') as f:
            w = csv.writer(f)
            for snr, blocks, errs, bler in results:
                w.writerow([mcs, qm, f"{R:.4f}", tbs, f"{snr:.1f}", f"{bler:.6f}", blocks, errs])
            f.flush()
        print(flush=True)

    print(f"\nAll done! Results saved to {out_csv}", flush=True)

if __name__ == '__main__':
    main()
