#!/usr/bin/env python3
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper

def main():
    n_prb = 30
    n_tx = 4
    n_rx = 4
    n_layers = 4
    n_sc = n_prb * 12
    n_ofdm_sym = 14
    scs = 30e3
    dmrs_symbols = [2, 7]
    pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in dmrs_symbols]
    num_dmrs_per_prb = len(dmrs_symbols) * 12

    mcs = 4
    Qm, R_x1024 = 2, 602
    target_coderate = R_x1024 / 1024.0

    tb_size, cb_size, num_cb, tb_crc_len, cb_crc_len, cw_length = calculate_tb_size(
        modulation_order=Qm,
        target_coderate=target_coderate,
        num_prbs=n_prb,
        num_ofdm_symbols=n_ofdm_sym,
        num_dmrs_per_prb=num_dmrs_per_prb,
        num_layers=n_layers,
        num_ov=0,
        tb_scaling=1.0,
        return_cw_length=True,
        verbose=False
    )
    k_info = int(tb_size.item())
    G = int(cw_length.sum().item())
    n_symbols_per_layer = G // (Qm * n_layers)

    print(f"4x4 Flat Rayleigh (1 tap) MCS{mcs}: k_info={k_info}, G={G}")

    device = "cpu"
    tb_enc = TBEncoder(target_tb_size=k_info, num_coded_bits=G, target_coderate=target_coderate,
                       num_bits_per_symbol=Qm, num_layers=n_layers, channel_type='PDSCH',
                       use_scrambler=True, verbose=False)
    tb_dec = TBDecoder(tb_enc, num_iter=20, hard_out=True)
    constellation = Constellation("qam", Qm)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)

    def solve_via_cholesky(A, B):
        L = torch.linalg.cholesky(A)
        Y = torch.linalg.solve_triangular(L, B, upper=False)
        X = torch.linalg.solve_triangular(L.conj().mT, Y, upper=True)
        return X

    def lmmse_equalize(H, y, no, n_layers):
        sqrt_no = torch.sqrt(torch.clamp(no, min=1e-12))
        H_tilde = H / sqrt_no.unsqueeze(-1).unsqueeze(-1)
        y_tilde = y / sqrt_no.unsqueeze(-1)
        Hh = H_tilde.conj().mT
        A = Hh @ H_tilde + torch.eye(n_layers, dtype=H.dtype, device=H.device).unsqueeze(0)
        Hhy = Hh @ y_tilde.unsqueeze(-1)
        x_hat = solve_via_cholesky(A, Hhy).squeeze(-1)
        A_inv = torch.linalg.inv(A)
        diag_A_inv = torch.diagonal(A_inv, dim1=-2, dim2=-1).real
        d = 1.0 - diag_A_inv
        d = torch.where(torch.abs(d) < 1e-12, torch.full_like(d, 1e-12), d)
        x_hat = x_hat / d
        eff_var = diag_A_inv / d
        return x_hat, eff_var

    sinr_points = np.arange(-4.0, 12.0, 2.0)
    batch_size = 5
    max_blocks = 50
    target_errors = 20

    for sinr_db in sinr_points:
        n_blocks = 0
        n_errors = 0
        while n_blocks < max_blocks and n_errors < target_errors:
            bs = min(batch_size, max_blocks - n_blocks)
            bits = torch.randint(0, 2, (bs, k_info), dtype=torch.float32, device=device)
            c = tb_enc(bits)
            symbols = mapper(c)
            x = symbols.reshape(bs, n_symbols_per_layer, n_layers)

            tx_grid = torch.zeros(bs, n_tx, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)
            re_idx = 0
            for sym in pdsch_symbols:
                for sc in range(n_sc):
                    for layer in range(n_layers):
                        tx_grid[:, layer, sym, sc] = x[:, re_idx, layer]
                    re_idx += 1

            esn0_lin = 10.0 ** (sinr_db / 10.0)
            noise_var = 1.0 / esn0_lin
            sigma_dim = np.sqrt(noise_var / 2.0)

            # 1-tap Rayleigh: H is same for all sc and sym
            H_flat = torch.randn(bs, n_rx, n_layers, dtype=torch.complex64) / np.sqrt(2.0)
            h_freq = H_flat.unsqueeze(1).unsqueeze(2).expand(-1, n_sc, n_ofdm_sym, -1, -1)

            rx_grid = torch.zeros_like(tx_grid)
            for rx in range(n_rx):
                for tx in range(n_tx):
                    h = h_freq[:, :, :, rx, tx].permute(0, 2, 1)
                    rx_grid[:, rx, :, :] += tx_grid[:, tx, :, :] * h
            rx_grid = rx_grid + torch.randn_like(rx_grid) * sigma_dim + 1j * torch.randn_like(rx_grid) * sigma_dim

            no_est = torch.tensor(noise_var, dtype=torch.float32, device=device).expand(bs)

            eq_list = []
            nv_list = []
            for sym in pdsch_symbols:
                for sc in range(n_sc):
                    H_sc = h_freq[:, sc, sym, :, :n_layers]
                    y_sc = rx_grid[:, :, sym, sc]
                    x_eq, eff_nv = lmmse_equalize(H_sc, y_sc, no_est, n_layers)
                    for l in range(n_layers):
                        eq_list.append(x_eq[:, l])
                        nv_list.append(eff_nv[:, l])
            eq_syms = torch.stack(eq_list, dim=1)
            nv = torch.stack(nv_list, dim=1)

            llr_cat = demapper(eq_syms, nv)[:, :G]
            b_hat, crc_ok = tb_dec(llr_cat)
            n_errors += int((crc_ok < 0.5).sum().item())
            n_blocks += bs

        bler = n_errors / n_blocks
        print(f"SNR {sinr_db:5.1f} dB: BLER = {bler:.4f} ({n_blocks} blocks, {n_errors} errors)")

if __name__ == "__main__":
    main()
