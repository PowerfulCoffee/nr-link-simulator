import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.mapping import Constellation, Mapper, Demapper
from sionna.phy.utils import hard_decisions

k_info = 2472
R = 910/1024
k_crc = k_info + 24
n_coded = 2808
print(f"k_info={k_info}, k_crc={k_crc}, n_coded(G)={n_coded}, R_eff={k_crc/n_coded:.4f}")

crc_enc = CRCEncoder("CRC24A")
crc_dec = CRCDecoder(crc_enc)
ldpc_enc = LDPC5GEncoder(k_crc, n_coded, num_bits_per_symbol=6)
ldpc_dec = LDPC5GDecoder(ldpc_enc, num_iter=20, hard_out=False, return_info=False)
constellation = Constellation("qam", 6)
mapper = Mapper(constellation=constellation)
demapper = Demapper("app", constellation=constellation)

print(f"Constellation power: {torch.mean(torch.abs(constellation.points)**2):.4f}")

bits = torch.randint(0, 2, (1, k_info), dtype=torch.float32)
bits_crc = crc_enc(bits)
cw = ldpc_enc(bits_crc)
x = mapper(cw)
print(f"Codeword bits: {cw.shape[1]}, symbols: {x.shape[1]}, x_power={torch.mean(torch.abs(x)**2):.4f}")

for snr_db in [20, 25, 30, 40]:
    noise_var = 10 ** (-snr_db/10)
    n_ok = 0
    N = 10
    for _ in range(N):
        bits = torch.randint(0, 2, (1, k_info), dtype=torch.float32)
        bits_crc = crc_enc(bits)
        cw = ldpc_enc(bits_crc)
        x = mapper(cw)
        noise = torch.randn_like(x) * np.sqrt(noise_var/2) + 1j*torch.randn_like(x) * np.sqrt(noise_var/2)
        y = x + noise
        llr = demapper(y, noise_var)
        cw_hat = ldpc_dec(llr)
        bits_hat, crc_valid = crc_dec(cw_hat)
        if crc_valid.item() > 0.5:
            n_ok += 1
    print(f"SNR={snr_db}dB: {n_ok}/{N} passed")
