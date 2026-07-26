#!/usr/bin/env python3
import os
import argparse
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mcs', type=int, required=True)
    parser.add_argument('--snr-start', type=float, required=True)
    parser.add_argument('--snr-end', type=float, required=True)
    parser.add_argument('--snr-step', type=float, required=True)
    parser.add_argument('--max-blocks', type=int, default=100)
    parser.add_argument('--target-errors', type=int, default=30)
    parser.add_argument('--n-prb', type=int, default=30)
    parser.add_argument('--n-rx', type=int, default=4)
    parser.add_argument('--n-layers', type=int, default=4)
    parser.add_argument('--channel', type=str, default="tdl", choices=["tdl", "awgn"])
    parser.add_argument('--delay-spread', type=float, default=100)
    parser.add_argument('--max-doppler', type=float, default=5.0)
    parser.add_argument('--dmrs-add-pos', type=int, default=1)
    parser.add_argument('--dmrs-duration', type=int, default=0)
    parser.add_argument('--ce-type', type=str, default="ideal", choices=["ideal", "ls"])
    parser.add_argument('--ldpc-iters', type=int, default=20)
    parser.add_argument('--batch-size', type=int, default=5)
    parser.add_argument('--mcs-table', type=int, default=2)
    parser.add_argument('--output', type=str, default=None)
    args = parser.parse_args()

    n_prb = args.n_prb
    n_tx = args.n_layers
    n_rx = args.n_rx
    n_layers = args.n_layers
    n_sc = n_prb * 12
    n_ofdm_sym = 14
    scs = 30e3

    if args.dmrs_add_pos == 1 and args.dmrs_duration == 0:
        dmrs_symbols = [2, 7]
    elif args.dmrs_add_pos == 1 and args.dmrs_duration == 1:
        dmrs_symbols = [2, 11]
    else:
        dmrs_symbols = [2]
    pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in dmrs_symbols]
    num_dmrs_per_prb = len(dmrs_symbols) * 12

    if args.mcs_table == 2:
        mcs_table2 = [
            (2,120),(2,193),(2,308),(2,449),(2,602),
            (4,378),(4,434),(4,490),(4,553),(4,616),(4,658),
            (6,466),(6,517),(6,567),(6,616),(6,666),(6,719),(6,772),(6,822),(6,873),
            (8,682.5),(8,711),(8,754),(8,797),(8,841),(8,885),(8,916.5),(8,948),
            (-1,0)
        ]
        Qm, R_x1024 = mcs_table2[args.mcs]
    else:
        mcs_table1 = [
            (2,78),(2,193),(2,449),
            (4,378),(4,490),(4,616),
            (6,466),(6,567),(6,666),(6,772),(6,873),
            (8,711),(8,797),(8,885),(8,948),
            (-1,0)
        ]
        Qm, R_x1024 = mcs_table1[args.mcs]
    target_coderate = R_x1024 / 1024.0
    mod_name = "qam"

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

    print(f"Configuration: MCS={args.mcs} (table {args.mcs_table}, {Qm}QAM, R={target_coderate:.4f}), "
          f"{n_tx}x{n_rx} {n_layers} layers, {n_prb} PRB, CE={args.ce_type}")
    print(f"  k_info={k_info}, G={G}, n_sym/layer={n_symbols_per_layer}")

    if args.channel == "tdl":
        print(f"  TDL-A channel, DS={args.delay_spread}ns, fd={args.max_doppler}Hz")
        taps_powers_db = np.array([
            -13.4, 0.0, -2.2, -4.0, -6.0, -8.2, -9.9, -10.5,
            -7.5, -15.9, -6.6, -16.7, -12.4, -15.2, -10.8,
            -11.3, -12.7, -16.2, -18.3, -18.9, -16.6, -19.9, -29.7])
        taps_norm_delays = np.array([
            0.0, 0.3819, 0.4025, 0.5868, 0.4610, 0.5375, 0.6708, 0.5750,
            0.7618, 1.5375, 1.8978, 2.2242, 2.1718, 2.4942, 2.5119,
            3.0582, 4.0810, 4.4579, 4.5695, 4.7966, 5.0066, 5.3043, 9.6586])
        n_taps = len(taps_powers_db)
        taps_powers_lin = 10.0 ** (taps_powers_db / 10.0)
        taps_powers_norm = taps_powers_lin / taps_powers_lin.sum()
        delay_spread_s = args.delay_spread * 1e-9
        taps_delays = taps_norm_delays * delay_spread_s

        half_sc = n_sc // 2
        freq_axis = np.zeros(n_sc)
        for sc in range(n_sc):
            if sc < half_sc:
                freq_axis[sc] = (sc - half_sc) * scs
            else:
                freq_axis[sc] = (sc - half_sc + 1) * scs

        phase_freq = np.exp(-1j * 2.0 * np.pi * freq_axis[:, None] * taps_delays[None, :])
        phase_freq_t = torch.tensor(phase_freq, dtype=torch.complex64)

    device = "cpu"
    tb_enc = TBEncoder(
        target_tb_size=k_info,
        num_coded_bits=G,
        target_coderate=target_coderate,
        num_bits_per_symbol=Qm,
        num_layers=n_layers,
        channel_type='PDSCH',
        use_scrambler=True,
        verbose=False
    )
    tb_dec = TBDecoder(tb_enc, num_iter=args.ldpc_iters, hard_out=True)
    constellation = Constellation(mod_name, Qm)
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

    results = []
    sinr_points = np.arange(args.snr_start, args.snr_end + 0.001, args.snr_step)
    batch_size = args.batch_size

    for sinr_db in sinr_points:
        n_blocks = 0
        n_errors = 0
        while n_blocks < args.max_blocks and n_errors < args.target_errors:
            bs = min(batch_size, args.max_blocks - n_blocks)
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

            if args.channel == "awgn":
                h_freq = torch.eye(n_layers, dtype=torch.complex64, device=device).unsqueeze(0).unsqueeze(0).unsqueeze(0).expand(bs, n_sc, n_ofdm_sym, -1, -1)
            else:
                tap_std = torch.sqrt(torch.tensor(taps_powers_norm, dtype=torch.float32))
                h_freq = torch.zeros(bs, n_sc, n_ofdm_sym, n_rx, n_tx, dtype=torch.complex64, device=device)
                for rx in range(n_rx):
                    for tx in range(n_tx):
                        taps_c = torch.randn(bs, n_taps, dtype=torch.complex64) * tap_std
                        h_sym = (taps_c.unsqueeze(1) * phase_freq_t.unsqueeze(0)).sum(dim=-1)
                        for sym in range(n_ofdm_sym):
                            h_freq[:, :, sym, rx, tx] = h_sym

            rx_grid = torch.zeros_like(tx_grid)
            for rx in range(n_rx):
                for tx in range(n_tx):
                    h = h_freq[:, :, :, rx, tx].permute(0, 2, 1)
                    rx_grid[:, rx, :, :] += tx_grid[:, tx, :, :] * h
            rx_grid = rx_grid + torch.randn_like(rx_grid) * sigma_dim + 1j * torch.randn_like(rx_grid) * sigma_dim

            no_est = torch.tensor(noise_var, dtype=torch.float32, device=device).expand(bs)
            h_est = h_freq

            eq_list = []
            nv_list = []
            for sym in pdsch_symbols:
                for sc in range(n_sc):
                    H_sc = h_est[:, sc, sym, :, :n_layers]
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

        bler = n_errors / n_blocks if n_blocks > 0 else 0.0
        results.append((sinr_db, n_blocks, n_errors, bler))
        print(f"SNR {sinr_db:5.1f} dB: BLER = {bler:.4f} ({n_blocks} blocks, {n_errors} errors)")

    if args.output:
        import csv
        with open(args.output, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['SINR_dB', 'Blocks', 'Errors', 'BLER'])
            for r in results:
                w.writerow(r)
        print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()
