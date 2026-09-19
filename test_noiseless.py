import torch
import numpy as np
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper

torch.manual_seed(42)
np.random.seed(42)

qm = 2
target_coderate = 449/1024
n_prb = 20
n_sc = n_prb*12
n_ofdm_sym = 14
n_layers = 1
pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in [2,7]]
n_pdsch_re = n_sc * len(pdsch_symbols)
G = n_pdsch_re * qm * n_layers

tb_config = calculate_tb_size(modulation_order=qm, target_coderate=target_coderate,
                              num_coded_bits=G, num_layers=n_layers, verbose=False)
tb_size = int(tb_config[0])
cw_length = int(torch.sum(tb_config[5]))
print(f"tb_size={tb_size}, G={G}, cw_length={cw_length}")

constellation = Constellation("qam", num_bits_per_symbol=qm, normalize=True)
mapper = Mapper(constellation=constellation)
demapper = Demapper("app", constellation=constellation)
tb_enc = TBEncoder(tb_size, G, target_coderate, qm, num_layers=n_layers, channel_type="PDSCH", verbose=False)
tb_dec = TBDecoder(tb_enc, num_bp_iter=20)

bits = torch.randint(0,2,(1,tb_size))
cw = tb_enc(bits)
print(f"cw shape: {cw.shape}, cw[:20]: {cw[0,:20]}")
syms = mapper(cw)
print(f"syms shape: {syms.shape}, power: {torch.mean(torch.abs(syms)**2):.3f}")
x = syms.reshape(-1)
llr = demapper(x.unsqueeze(0), no=1e-10)
print(f"llr shape: {llr.shape}, llr[:20]: {llr[0,:20]}")
dec, crc_status = tb_dec(llr.reshape(1,-1))
print(f"crc_status: {crc_status}")
dec_bits = (dec[0] < 0).long()
print(f"bits match: {torch.equal(bits[0], dec_bits[0])}, bit errors: {torch.sum(bits[0]!=dec_bits[0]).item()}")
