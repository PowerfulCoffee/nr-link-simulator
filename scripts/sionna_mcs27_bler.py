import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.mapping import Constellation, Mapper, Demapper
import time

k_info = 2472
k_crc = k_info + 24
n_coded = 2808
Qm = 6

crc_enc = CRCEncoder("CRC24A")
crc_dec = CRCDecoder(crc_enc)
ldpc_enc = LDPC5GEncoder(k_crc, n_coded, num_bits_per_symbol=Qm)
ldpc_dec = LDPC5GDecoder(ldpc_enc, num_iter=20, hard_out=True, return_info=False)
constellation = Constellation("qam", Qm)
mapper = Mapper(constellation=constellation)
demapper = Demapper("app", constellation=constellation)

snrs = [16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0]
max_blocks = 5000
min_errors = 200
results = {}

print("Sionna MCS27 fine BLER sweep (k={}, n={}, Qm=6, AWGN, Es/N0)".format(k_info, n_coded), flush=True)
for snr_db in snrs:
    no = 10 ** (-snr_db/10)
    n_errors = 0
    n_total = 0
    t0 = time.time()
    while n_total < max_blocks and n_errors < min_errors:
        bs = min(200, max_blocks - n_total)
        bits = torch.randint(0, 2, (bs, k_info), dtype=torch.float32)
        bits_crc = crc_enc(bits)
        cw = ldpc_enc(bits_crc)
        x = mapper(cw)
        noise = torch.randn_like(x) * np.sqrt(no/2) + 1j*torch.randn_like(x) * np.sqrt(no/2)
        y = x + noise
        llr = demapper(y, no)
        cw_hat = ldpc_dec(llr)
        _, crc_valid = crc_dec(cw_hat)
        n_errors += (crc_valid < 0.5).sum().item()
        n_total += bs
    bler = n_errors / n_total
    dt = time.time() - t0
    results[snr_db] = (n_errors, n_total, bler, dt)
    print(f"SNR={snr_db:5.1f} dB: BLER={bler:.4f} ({n_errors}/{n_total}, {dt:.1f}s)", flush=True)
    if bler < 1e-3 and snr_db > 17:
        break

os.makedirs("results", exist_ok=True)
with open("results/sionna_mcs27_awgn.txt", "w") as f:
    f.write("snr_db,n_errors,n_total,bler\n")
    for snr_db, (ne, nt, bler, _) in sorted(results.items()):
        f.write(f"{snr_db},{ne},{nt},{bler:.6f}\n")
print("Done.")
