import torch
import numpy as np
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper

torch.manual_seed(42); np.random.seed(42)
qm=2; R=363/1024; n_prb=30; n_sc=360; n_ofdm_sym=14
n_layers=4; n_rx=4; n_tx=4
dmrs_symbols=[2,7]
pdsch_symbols=[s for s in range(n_ofdm_sym) if s not in dmrs_symbols]

tb = calculate_tb_size(modulation_order=qm, target_coderate=R, num_prbs=n_prb,
    num_ofdm_symbols=n_ofdm_sym, num_dmrs_per_prb=24, num_layers=n_layers, verbose=False)
tb_size=int(tb[0]); cw_length=int(torch.sum(tb[5]))
G = n_sc*len(pdsch_symbols)*qm*n_layers
print(f"tb={tb_size}, G={G}, cw={cw_length}")

const = Constellation("qam", num_bits_per_symbol=qm, normalize=True)
map_ = Mapper(constellation=const)
demap = Demapper("app", constellation=const)
enc = TBEncoder(tb_size, G, R, qm, num_layers=n_layers, channel_type="PDSCH", verbose=False)
dec = TBDecoder(enc, num_bp_iter=20)

bits = torch.randint(0,2,(1,tb_size))
cw = enc(bits)
syms = map_(cw)

syms_grid = syms.reshape(1, n_layers, n_sc, len(pdsch_symbols)).permute(0,1,3,2)

tx_grid = torch.zeros(1, n_layers, n_ofdm_sym, n_sc, dtype=torch.complex64)
for i, s in enumerate(pdsch_symbols):
    tx_grid[:, :, s, :] = syms_grid[:, :, i, :]

h_time = (torch.randn(1,n_ofdm_sym,n_sc,n_rx,n_tx) + 1j*torch.randn(1,n_ofdm_sym,n_sc,n_rx,n_tx))/np.sqrt(2)
rx_signal = torch.zeros(1,n_rx,n_ofdm_sym,n_sc,dtype=torch.complex64)
for sym in range(n_ofdm_sym):
    for sc in range(n_sc):
        rx_signal[0,:,sym,sc] = h_time[0,sym,sc,:,:] @ tx_grid[0,:,sym,sc]

rx_data = torch.zeros(1,n_rx,len(pdsch_symbols),n_sc,dtype=torch.complex64)
h_data = torch.zeros(1,n_rx,len(pdsch_symbols),n_sc,n_layers,dtype=torch.complex64)
for i,s in enumerate(pdsch_symbols):
    rx_data[:,:,i,:] = rx_signal[:,:,s,:]
    h_data[:,:,i,:,:] = h_time[:,s,:,:,:].permute(0,2,1,3)

H_t = h_data.permute(0,2,3,1,4).reshape(len(pdsch_symbols)*n_sc, n_rx, n_layers)
y_t = rx_data.permute(0,2,3,1).reshape(len(pdsch_symbols)*n_sc, n_rx)
x_hat = torch.zeros_like(y_t[:,:n_layers])
for k in range(H_t.shape[0]):
    x_hat[k] = torch.linalg.lstsq(H_t[k], y_t[k].unsqueeze(1)).solution.squeeze(1)

x_hat_grid = x_hat.reshape(1,len(pdsch_symbols),n_sc,n_layers).permute(0,3,1,2)
rx_syms = x_hat_grid.reshape(1, n_layers*len(pdsch_symbols)*n_sc)[:,:cw_length]

print(f"syms shape: {syms.shape}, rx_syms shape: {rx_syms.shape}")
print(f"max |syms - rx_syms|: {torch.max(torch.abs(syms - rx_syms)):.10f}")
print(f"syms[:5]: {syms[0,:5]}")
print(f"rx_syms[:5]: {rx_syms[0,:5]}")
print(f"syms is rx_syms: {torch.allclose(syms, rx_syms)}")

llr_direct = demap(syms, no=1e-10)
d_direct, crc_direct = dec(llr_direct)
print(f"direct: crc={crc_direct}, match={torch.equal(bits[0], d_direct[0].long())}")

llr_rx = demap(rx_syms, no=1e-10)
d_rx, crc_rx = dec(llr_rx)
print(f"rx: crc={crc_rx}, match={torch.equal(bits[0], d_rx[0].long())}")
