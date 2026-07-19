import os, csv, math, numpy as np, torch
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from sionna.phy.nr import TBEncoder, TBDecoder
from sionna.phy.mapping import Mapper, Demapper, Constellation

n_prb = 3; n_re_per_prb = 156; n_layers = 1
num_ldpc_iters = 25; max_blocks = 3000; target_errors = 80; batch_size = 400

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

refine_mcs = [2, 4, 7, 9, 10, 11, 13, 15, 18, 20, 22, 24, 26, 27]
snr_offsets = list(range(-3, 4))

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, 'sionna_refined_key_mcs.csv')
with open(out_csv, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['mcs','qm','rate','tbs','snr_db','bler','n_blocks','n_errors'])

for mcs in refine_mcs:
    qm, r1024 = MCS_TABLE[mcs]
    R = r1024 / 1024.0
    tbs = CPP_TBS[mcs]
    n_coded_bits = n_prb * n_re_per_prb * qm * n_layers
    print(f"=== Refining MCS{mcs} (Qm={qm}, R={R:.3f}) ===", flush=True)

    tb_enc = TBEncoder(target_tb_size=tbs, num_coded_bits=n_coded_bits, target_coderate=R,
        num_bits_per_symbol=qm, num_layers=n_layers, channel_type='PUSCH', use_scrambler=True, verbose=False)
    tb_dec = TBDecoder(tb_enc, num_iter=num_ldpc_iters, hard_out=True, verbose=False)
    constellation = Constellation("qam", qm)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)

    cpp_bler10_est = {
        2:0.71, 4:3.85, 7:6.87, 9:8.80, 10:9.26, 11:10.92,
        13:12.88, 15:14.44, 18:16.92, 20:18.91, 22:20.89, 24:22.88, 26:24.84, 27:25.90
    }
    center = cpp_bler10_est.get(mcs, 10.0)
    snr_range = [center + off for off in snr_offsets]

    for esn0_db in snr_range:
        no_val = 1.0 / (10 ** (esn0_db / 10.0))
        no_t = torch.tensor(no_val, dtype=torch.float32)
        n_errors = 0; n_blocks = 0
        while n_blocks < max_blocks and n_errors < target_errors:
            bs = min(batch_size, max_blocks - n_blocks)
            b = torch.randint(0, 2, (bs, tbs), dtype=torch.float32)
            c = tb_enc(b); x = mapper(c)
            ns = math.sqrt(no_val/2.0)
            noise = torch.complex(torch.randn_like(x,dtype=torch.float32)*ns, torch.randn_like(x,dtype=torch.float32)*ns).to(x.dtype)
            y = x + noise
            llr = demapper(y, no_t)
            _, crc_ok = tb_dec(llr)
            n_errors += bs - int(torch.sum(crc_ok.to(torch.float32)).item())
            n_blocks += bs
        bler = n_errors / n_blocks
        print(f"  SNR {esn0_db:.1f} dB: BLER={bler:.4f} ({n_blocks} blocks, {n_errors} errors)", flush=True)
        with open(out_csv, 'a', newline='') as f:
            w = csv.writer(f)
            w.writerow([mcs, qm, f"{R:.4f}", tbs, f"{esn0_db:.1f}", f"{bler:.6f}", n_blocks, n_errors])
        if bler < 1e-4 and n_blocks >= 400: break

print(f"\nDone! Refined results saved to {out_csv}", flush=True)
