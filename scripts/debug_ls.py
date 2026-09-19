import torch
import numpy as np
import sys
sys.path.insert(0, '/workspace/nr-link-simulator/scripts')
from sionna_tdl_crosscheck import (
    build_dmrs_grid, generate_tdl_channel, mmse_equalize
)

n_prb = 30
n_sc = n_prb * 12
n_ofdm_sym = 14
n_layers = 4
n_rx = 4
n_tx = 4
dmrs_symbols = [2, 7]
pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in dmrs_symbols]

dmrs_ref, dmrs_sc_per_layer = build_dmrs_grid(n_sc, n_ofdm_sym, n_layers, dmrs_symbols)
dmrs_ref = dmrs_ref.to(torch.complex64)

bs = 2
snr = 100
noise_var = 1.0 / (10 ** (snr / 10.0))

syms_dummy = torch.randn(bs, n_layers, len(pdsch_symbols), n_sc) + 1j * torch.randn(bs, n_layers, len(pdsch_symbols), n_sc)
syms_dummy = syms_dummy / syms_dummy.abs().clamp(min=1e-8)

tx_grid = torch.zeros(bs, n_tx, n_ofdm_sym, n_sc, dtype=torch.complex64)
for i, s in enumerate(pdsch_symbols):
    tx_grid[:, :n_layers, s, :] = syms_dummy[:, :, i, :]
for l in range(n_layers):
    for sym in dmrs_symbols:
        for sc in dmrs_sc_per_layer[l]:
            tx_grid[:, l, sym, sc] = dmrs_ref[l, sym, sc]

h_t = generate_tdl_channel(bs, n_sc, n_ofdm_sym, n_rx, n_tx, 15, 5, scs=15000, seed=42)
print(f"h_t shape: {h_t.shape}")
print(f"h_t power per element: {h_t.abs().pow(2).mean().item():.4f}")

rx_signal = torch.einsum('brscl,blsc->brsc', h_t[:, :, :, :, :n_layers], tx_grid)
noise = (torch.randn_like(rx_signal) + 1j * torch.randn_like(rx_signal)) * np.sqrt(noise_var / 2.0)
rx_signal = rx_signal + noise

exec(open('/workspace/nr-link-simulator/scripts/sionna_tdl_crosscheck.py').read().split('def main():')[0])

h_ideal = h_t[:, :, :, :, :n_layers]
h_est_ls = ls_channel_estimate(rx_signal, dmrs_ref, dmrs_symbols, dmrs_sc_per_layer, n_sc, n_ofdm_sym, n_rx, n_layers)

print(f"\nh_est_ls shape: {h_est_ls.shape}")
print(f"h_ideal shape: {h_ideal.shape}")

h_ls_pdsch = h_est_ls[:, :, pdsch_symbols, :, :]
h_id_pdsch = h_ideal[:, :, pdsch_symbols, :, :]

mse = (h_ls_pdsch - h_id_pdsch).abs().pow(2).mean().item()
sig_pow = h_id_pdsch.abs().pow(2).mean().item()
nmse_db = 10 * np.log10(mse / sig_pow)
print(f"\nNMSE of LS CE vs Ideal (SNR={snr}dB): {nmse_db:.2f} dB")

for sym_idx, sym in enumerate(pdsch_symbols):
    err = (h_ls_pdsch[:, :, sym_idx, :, :] - h_id_pdsch[:, :, sym_idx, :, :]).abs().pow(2).mean().item()
    pow_s = h_id_pdsch[:, :, sym_idx, :, :].abs().pow(2).mean().item()
    nmse = 10 * np.log10(err / pow_s) if err > 0 else -100
    print(f"  Symbol {sym}: NMSE = {nmse:.2f} dB")

for l in range(n_layers):
    err = (h_ls_pdsch[:, :, :, :, l] - h_id_pdsch[:, :, :, :, l]).abs().pow(2).mean().item()
    pow_s = h_id_pdsch[:, :, :, :, l].abs().pow(2).mean().item()
    nmse = 10 * np.log10(err / pow_s) if err > 0 else -100
    print(f"  Layer {l}: NMSE = {nmse:.2f} dB")

for sym in dmrs_symbols:
    err = (h_est_ls[:, :, sym, :, :] - h_ideal[:, :, sym, :, :]).abs().pow(2).mean().item()
    pow_s = h_ideal[:, :, sym, :, :].abs().pow(2).mean().item()
    nmse = 10 * np.log10(err / pow_s) if err > 0 else -100
    print(f"  DMRS Symbol {sym}: NMSE = {nmse:.2f} dB")
