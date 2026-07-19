import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.mapping import Constellation, Mapper, Demapper
from sionna.phy.utils import hard_decisions

k_info = 2472
k_crc = k_info + 24
n_coded = 2808

crc_enc = CRCEncoder("CRC24A")
crc_dec = CRCDecoder(crc_enc)
ldpc_enc = LDPC5GEncoder(k_crc, n_coded, num_bits_per_symbol=6)
ldpc_dec = LDPC5GDecoder(ldpc_enc, num_iter=50, hard_out=True, return_info=False)
constellation = Constellation("qam", 6)
mapper = Mapper(constellation=constellation)
demapper = Demapper("app", constellation=constellation)

bits = torch.randint(0, 2, (1, k_info), dtype=torch.float32)
bits_crc = crc_enc(bits)
cw = ldpc_enc(bits_crc)
x = mapper(cw)

y = x
no = torch.tensor(1e-4, dtype=torch.float32)
llr = demapper(y, no)
cw_hat = ldpc_dec(llr)
print(f"cw_hat shape: {cw_hat.shape}")
errs = (cw_hat.long() != bits_crc.long()).sum().item()
print(f"Hard errors in k_crc={k_crc} bits: {errs}")
if errs > 0:
    err_idx = (cw_hat.long() != bits_crc.long())[0].nonzero()
    print(f"Error positions: {err_idx[:20].flatten().tolist()}")
bits_hat, crc_valid = crc_dec(cw_hat)
print(f"CRC valid: {crc_valid.item()}, info bit errors: {(bits_hat.long() != bits.long()).sum().item()}")
