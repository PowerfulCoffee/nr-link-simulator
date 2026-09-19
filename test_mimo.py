import torch
import numpy as np
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper

torch.manual_seed(42)
np.random.seed(42)

qm = 2
target_coderate = 363/1024
n_prb = 30
n_sc = n_prb*12
n_ofdm_sym = 14
n_layers = 4
n_rx = 4
n_tx = 4
dmrs_symbols = [2,7]
pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in dmrs_symbols]
n_dmrs_per_prb = 24
dmrs_sc_per_layer = [list(range(l, n_sc, 4)) for l in range(n_layers)]

tb_config = calculate_tb_size(modulation_order=qm, target_coderate=target_coderate,
                              num_prbs=n_prb, num_ofdm_symbols=n_ofdm_sym,
                              num_dmrs_per_prb=n_dmrs_per_prb, num_layers=n_layers, verbose=False)
tb_size = int(tb_config[0])
cw_length = int(torch.sum(tb_config[5]))
n_pdsch_re_per_layer = n_sc * len(pdsch_symbols)
G = n_pdsch_re_per_layer * qm * n_layers
print(f"tb_size={tb_size}, G={G}, cw_length={cw_length}, n_pdsch_re/layer={n_pdsch_re_per_layer}")

constellation = Constellation("qam", num_bits_per_symbol=qm, normalize=True)
mapper = Mapper(constellation=constellation)
demapper = Demapper("app", constellation=constellation)
tb_enc = TBEncoder(tb_size, G, target_coderate, qm, num_layers=n_layers, channel_type="PDSCH", verbose=False)
tb_dec = TBDecoder(tb_enc, num_bp_iter=20)

bits = torch.randint(0,2,(1,tb_size))
cw = tb_enc(bits)
syms = mapper(cw)
print(f"syms shape after mapper: {syms.shape}")
syms = syms.reshape(1, n_layers, n_sc, len(pdsch_symbols))
syms = syms.permute(0, 1, 3, 2)
print(f"syms shape after reshape: {syms.shape}")

tx_grid = torch.zeros(1, n_layers, n_ofdm_sym, n_sc, dtype=torch.complex64)
for i, s in enumerate(pdsch_symbols):
    tx_grid[:, :, s, :] = syms[:, :, i, :]

h_time = (torch.randn(1, n_ofdm_sym, n_sc, n_rx, n_tx) + 
         1j*torch.randn(1, n_ofdm_sym, n_sc, n_rx, n_tx))/np.sqrt(2.0)

rx_signal = torch.zeros(1, n_rx, n_ofdm_sym, n_sc, dtype=torch.complex64)
for sym in range(n_ofdm_sym):
    for sc in range(n_sc):
        x = tx_grid[0, :, sym, sc]
        H = h_time[0, sym, sc, :, :]
        rx_signal[0, :, sym, sc] = H @ x

rx_data = torch.zeros(1, n_rx, len(pdsch_symbols), n_sc, dtype=torch.complex64)
h_data = torch.zeros(1, n_rx, len(pdsch_symbols), n_sc, n_layers, dtype=torch.complex64)
for i, s in enumerate(pdsch_symbols):
    rx_data[:, :, i, :] = rx_signal[:, :, s, :]
    h_data[:, :, i, :, :] = h_time[:, s, :, :, :].permute(0,2,1,3)

print(f"rx_data shape: {rx_data.shape}, h_data shape: {h_data.shape}")

H_t = h_data.permute(0,2,3,1,4).reshape(len(pdsch_symbols)*n_sc, n_rx, n_layers)
y_t = rx_data.permute(0,2,3,1).reshape(len(pdsch_symbols)*n_sc, n_rx)
x_hat = torch.zeros_like(y_t[:, :n_layers])
for k in range(H_t.shape[0]):
    Hk = H_t[k]
    yk = y_t[k]
    x_hat[k] = torch.linalg.lstsq(Hk, yk.unsqueeze(1)).solution.squeeze(1)

x_hat = x_hat.reshape(1, len(pdsch_symbols), n_sc, n_layers).permute(0,3,1,2)
rx_syms = x_hat.reshape(1, n_layers*len(pdsch_symbols)*n_sc)[:, :cw_length]
print(f"rx_syms shape: {rx_syms.shape}")

x_ref = syms.permute(0,1,2,3).reshape(1, n_layers*len(pdsch_symbols)*n_sc)[:, :cw_length]
print(f"x_ref shape: {x_ref.shape}")
print(f"x power: {torch.mean(torch.abs(x_ref)**2):.4f}")
print(f"x_hat power: {torch.mean(torch.abs(rx_syms)**2):.4f}")
err = rx_syms - x_ref
print(f"NMSE (noiseless): {torch.mean(torch.abs(err)**2)/torch.mean(torch.abs(x_ref)**2):.6f}")
print(f"First 5 x: {x_ref[0,:5]}")
print(f"First 5 x_hat: {rx_syms[0,:5]}")

llr = demapper(rx_syms, no=1e-10)
dec_bits, crc_ok = tb_dec(llr.reshape(1,-1))
print(f"crc_ok: {crc_ok}")
print(f"bit match: {torch.equal(bits[0], dec_bits[0].long())}")
