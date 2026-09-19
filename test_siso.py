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
dmrs_symbols = [2,7]
pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in dmrs_symbols]
n_pdsch_re = n_sc * len(pdsch_symbols)
G = n_pdsch_re * qm * n_layers

tb_config = calculate_tb_size(
    modulation_order=qm, target_coderate=target_coderate,
    num_coded_bits=G, num_layers=n_layers, verbose=True)
tb_size = int(tb_config[0])
cw_length = int(torch.sum(tb_config[5]))
print(f"tb_size={tb_size}, G={G}, cw_length={cw_length}")

constellation = Constellation("qam", num_bits_per_symbol=qm, normalize=True)
mapper = Mapper(constellation=constellation)
demapper = Demapper("app", constellation=constellation)
tb_enc = TBEncoder(tb_size, G, target_coderate, qm, num_layers=n_layers, channel_type="PDSCH", verbose=False)
tb_dec = TBDecoder(tb_enc, num_bp_iter=20)

for snr_db in [0, 5, 10, 20]:
    noise_var = 1.0/(10**(snr_db/10))
    n_errors = 0
    n_blocks = 10
    for _ in range(n_blocks):
        bits = torch.randint(0,2,(1,tb_size))
        cw = tb_enc(bits)
        syms = mapper(cw)
        x = syms.reshape(n_pdsch_re)
        y = x + (torch.randn_like(x) + 1j*torch.randn_like(x))*np.sqrt(noise_var/2)
        llr = demapper(y.unsqueeze(0), no=noise_var)
        dec, _ = tb_dec(llr.reshape(1,-1))
        dec_bits = (dec[0] < 0).long()
        if not torch.equal(bits[0], dec_bits[0]):
            n_errors += 1
    print(f"SNR={snr_db}dB: BLER={n_errors/n_blocks:.3f}")
