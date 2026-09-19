import torch
import numpy as np
from sionna_tdl_crosscheck import (MCS_TABLE_2, build_dmrs_grid, get_fft_size,
                                   ls_channel_estimate, mmse_equalize)
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper

torch.manual_seed(42)

def run_test(n_layers, n_rx, n_tx, snr, ce_type='ideal', use_mimo=True, n_blocks=50):
    mcs = 4
    qm, R_x1024 = MCS_TABLE_2[mcs]
    target_coderate = R_x1024 / 1024.0
    n_prb = 30
    n_sc = n_prb * 12
    n_ofdm_sym = 14
    scs = 15000
    fft_size = get_fft_size(n_sc)

    dmrs_symbols = [2, 7]
    pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in dmrs_symbols]
    n_dmrs_per_prb = 12 * len(dmrs_symbols)

    tb_config = calculate_tb_size(
        modulation_order=qm, target_coderate=target_coderate,
        num_prbs=n_prb, num_ofdm_symbols=n_ofdm_sym,
        num_dmrs_per_prb=n_dmrs_per_prb, num_layers=n_layers, verbose=False)
    tb_size = int(tb_config[0])
    cw_length = int(torch.sum(tb_config[5]))

    constellation = Constellation("qam", num_bits_per_symbol=qm, normalize=True)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)
    tb_enc = TBEncoder(tb_size, cw_length, target_coderate, qm,
                       num_layers=n_layers, channel_type="PDSCH", verbose=False)
    tb_dec = TBDecoder(tb_enc, num_bp_iter=20)

    dmrs_ref, dmrs_sc_per_layer = build_dmrs_grid(n_sc, n_ofdm_sym, n_layers, dmrs_symbols)
    dmrs_ref = dmrs_ref.to(torch.complex64)
    n_re_per_layer = len(pdsch_symbols) * n_sc

    noise_var = 1.0 / (10 ** (snr / 10.0))
    bs = 5
    n_err = 0
    n_tot = 0
    chan_seed = 42

    for _ in range(n_blocks // bs):
        bits = torch.randint(0, 2, (bs, tb_size))
        cw = tb_enc(bits)
        syms = mapper(cw)
        syms = syms.reshape(bs, n_re_per_layer, n_layers)
        syms = syms.permute(0, 2, 1).reshape(bs, n_layers, len(pdsch_symbols), n_sc)

        tx_grid = torch.zeros(bs, n_tx, n_ofdm_sym, n_sc, dtype=torch.complex64)
        for i, s in enumerate(pdsch_symbols):
            tx_grid[:, :n_layers, s, :] = syms[:, :, i, :]
        for l in range(n_layers):
            for sym in dmrs_symbols:
                for sc in dmrs_sc_per_layer[l]:
                    tx_grid[:, l, sym, sc] = dmrs_ref[l, sym, sc]

        if use_mimo:
            from sionna_tdl_crosscheck import generate_tdl_channel
            h_t = generate_tdl_channel(bs, n_sc, n_ofdm_sym, n_rx, n_tx, 15, 5, scs, seed=chan_seed)
            chan_seed += 1
            h_used = h_t[:, :, :, :, :n_layers]
            rx_signal = torch.einsum('brscl,blsc->brsc', h_used, tx_grid)
        else:
            # Identity channel
            h_used = torch.zeros(bs, n_rx, n_ofdm_sym, n_sc, n_layers, dtype=torch.complex64)
            for r in range(n_rx):
                for t in range(n_layers):
                    if r == t:
                        h_used[:, r, :, :, t] = 1.0
            rx_signal = torch.einsum('brscl,blsc->brsc', h_used, tx_grid)

        noise = (torch.randn_like(rx_signal) + 1j * torch.randn_like(rx_signal)) * np.sqrt(noise_var / 2.0)
        rx_signal = rx_signal + noise

        if ce_type == 'ideal':
            h_est = h_used
        else:
            h_est = ls_channel_estimate(rx_signal, dmrs_ref, dmrs_symbols, dmrs_sc_per_layer,
                                        n_sc, n_ofdm_sym, n_rx, n_layers, scs)

        rx_data = torch.zeros(bs, n_rx, len(pdsch_symbols), n_sc, dtype=torch.complex64)
        h_data = torch.zeros(bs, n_rx, len(pdsch_symbols), n_sc, n_layers, dtype=torch.complex64)
        for i, s in enumerate(pdsch_symbols):
            rx_data[:, :, i, :] = rx_signal[:, :, s, :]
            h_data[:, :, i, :, :] = h_est[:, :, s, :, :]

        x_hat, eff_noise = mmse_equalize(rx_data, h_data, noise_var, n_layers)
        x_hat = x_hat.reshape(bs, n_layers, n_re_per_layer).permute(0, 2, 1).reshape(bs, n_layers * n_re_per_layer)
        rx_noise = eff_noise.reshape(bs, n_layers, n_re_per_layer).permute(0, 2, 1).reshape(bs, n_layers * n_re_per_layer)
        rx_syms = x_hat[:, :cw_length]
        rx_noise = rx_noise[:, :cw_length]

        llr = demapper(rx_syms, no=rx_noise)
        dec_out, crc_ok = tb_dec(llr.reshape(bs, cw_length))
        crc_ok = crc_ok[0].bool()
        n_err += bs - crc_ok.sum().item()
        n_tot += bs

    return n_err / n_tot

print("=== Identity Channel (AWGN) ===")
for snr in [4, 6, 8, 10]:
    bler = run_test(1, 1, 1, snr, 'ideal', use_mimo=False, n_blocks=50)
    print(f"SISO 1x1 Ideal SNR={snr}dB: BLER={bler:.3f}")

print("\n=== 4x4 Identity Channel ===")
for snr in [8, 10, 12]:
    bler = run_test(4, 4, 4, snr, 'ideal', use_mimo=False, n_blocks=20)
    print(f"MIMO 4x4 Ideal SNR={snr}dB: BLER={bler:.3f}")
