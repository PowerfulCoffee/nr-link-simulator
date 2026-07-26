"""
Sionna 4x4 MIMO TDL-A BLER simulation for cross-check with C++ nr-link-simulator
Config: 30PRB, 4Rx, rank4, DMRS Type1 add-pos=1 (symbols 2&7, entire symbols skipped),
        TDL-A DS=100ns, MCS Table 2 for 256QAM support
Power normalization matches C++:
- Tx symbols: unit average power
- TDL channel: each tap normalized so total average power gain = 1
- Noise variance: 1 / SNR_lin (each complex dim: N(0, noise_var/2) per real/imag)
"""
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
import time
import sys

from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper


MCS_TABLE2_PDSCH = [
    (2, 120), (2, 193), (2, 308), (2, 449), (2, 602),
    (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658),
    (6, 466), (6, 517), (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873),
    (8, 682.5), (8, 711), (8, 754), (8, 797), (8, 841), (8, 885), (8, 916.5), (8, 948),
]


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
    HhH = Hh @ H_tilde
    A = HhH + torch.eye(n_layers, dtype=H.dtype, device=H.device).unsqueeze(0)
    Hhy = Hh @ y_tilde.unsqueeze(-1)

    x_hat = solve_via_cholesky(A, Hhy).squeeze(-1)

    A_inv = torch.linalg.inv(A)
    diag_A_inv = torch.diagonal(A_inv, dim1=-2, dim2=-1).real
    d = 1.0 - diag_A_inv
    d = torch.where(torch.abs(d) < 1e-12, torch.full_like(d, 1e-12), d)
    x_hat = x_hat / d
    eff_var = diag_A_inv / d

    return x_hat, eff_var


def main():
    if len(sys.argv) < 2:
        print("Usage: python sionna_tdl_crosscheck.py <mcs> <perfect_csi=0/1> [out_file]")
        sys.exit(1)

    mcs = int(sys.argv[1])
    perfect_csi = int(sys.argv[2]) == 1 if len(sys.argv) > 2 else False
    out_file = sys.argv[3] if len(sys.argv) > 3 else None

    n_prb = 30
    n_tx = 4
    n_rx = 4
    n_layers = 4
    n_sc = n_prb * 12
    n_ofdm_sym = 14

    scs = 30e3
    delay_spread = 100e-9
    max_doppler = 5.0
    mcs_table = 2

    dmrs_symbols = [2, 7]
    pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in dmrs_symbols]
    num_dmrs_per_prb = len(dmrs_symbols) * 12

    Qm, R_x1024 = MCS_TABLE2_PDSCH[mcs]
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

    mod_name = {2: "QPSK", 4: "16QAM", 6: "64QAM", 8: "256QAM"}[Qm]
    print(f"MCS{mcs} (table {mcs_table}): Qm={Qm} ({mod_name}), R={target_coderate:.4f}, k_info={k_info}, n_coded={G}")
    print(f"Config: {n_prb}PRB, {n_tx}Tx{n_rx}Rx, {n_layers}layers, TDL-A DS={delay_spread*1e9:.0f}ns, fd={max_doppler}Hz")
    print(f"DMRS: Type1, symbols={dmrs_symbols}, PDSCH syms={pdsch_symbols}")
    print(f"Channel Est: {'Perfect/Ideal' if perfect_csi else 'LS + linear interp'}")

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
    tb_dec = TBDecoder(tb_enc, num_iter=20, hard_out=True)
    constellation = Constellation("qam", Qm)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)

    half_sc = n_sc // 2
    freq_axis = np.zeros(n_sc, dtype=np.float64)
    for sc in range(n_sc):
        if sc < half_sc:
            freq_axis[sc] = (sc - half_sc) * scs
        else:
            freq_axis[sc] = (sc - half_sc + 1) * scs
    freq_axis = torch.tensor(freq_axis, dtype=torch.float64, device=device)

    taps_powers_db = np.array([
        -13.4, 0.0, -2.2, -4.0, -6.0, -8.2, -9.9, -10.5,
        -7.5, -15.9, -6.6, -16.7, -12.4, -15.2, -10.8,
        -11.3, -12.7, -16.2, -18.3, -18.9, -16.6, -19.9, -29.7])
    taps_norm_delays = np.array([
        0.0, 0.3819, 0.4025, 0.5868, 0.4610, 0.5375, 0.6708, 0.5750,
        0.7618, 1.5375, 1.8978, 2.2242, 2.1718, 2.4942, 2.5119,
        3.0582, 4.0810, 4.4579, 4.5695, 4.7966, 5.0066, 5.3043, 9.6586])
    taps_powers_lin = 10.0 ** (taps_powers_db / 10.0)
    taps_powers_norm = taps_powers_lin / taps_powers_lin.sum()

    taps_powers_t = torch.tensor(taps_powers_norm, dtype=torch.float32, device=device)
    taps_delays_t = torch.tensor(taps_norm_delays * delay_spread, dtype=torch.float64, device=device)
    n_taps = len(taps_powers_norm)
    tap_std = torch.sqrt(taps_powers_t / 2.0).float()

    phase_freq = -2.0 * np.pi * freq_axis[:, None] * taps_delays_t[None, :]
    phase_freq_t = torch.complex(torch.cos(phase_freq).float(), torch.sin(phase_freq).float())

    def get_ofdm_time(sym_idx):
        tu = 1.0 / scs
        t = 0.0
        fft_size_local = 512
        for s in range(sym_idx):
            is_first = (s == 0)
            cp = 144 * fft_size_local // 2048
            if is_first:
                cp += 16 * fft_size_local // 2048
            t += tu * (fft_size_local + cp) / fft_size_local
        is_first = (sym_idx == 0)
        cp = 144 * fft_size_local // 2048
        if is_first:
            cp += 16 * fft_size_local // 2048
        t += tu * (fft_size_local + cp // 2) / fft_size_local
        return t

    sym_times = np.array([get_ofdm_time(s) for s in range(n_ofdm_sym)])
    sym_times_t = torch.tensor(sym_times, dtype=torch.float32, device=device)

    def run_block(batch_size, esn0_db):
        bits = torch.randint(0, 2, (batch_size, k_info), dtype=torch.float32, device=device)
        c = tb_enc(bits)
        symbols = mapper(c)

        x = symbols.reshape(batch_size, n_symbols_per_layer, n_layers)

        # Generate frequency-domain channel per (rx, tx, sym) with time correlation via Doppler
        # Using block fading per OFDM symbol (independent per symbol, but phase evolves with Doppler)
        h_freq = torch.zeros(batch_size, n_sc, n_ofdm_sym, n_rx, n_tx, dtype=torch.complex64, device=device)

        for rx in range(n_rx):
            for tx in range(n_tx):
                # Generate independent tap coefficients per (rx, tx) pair
                taps_c = torch.randn(batch_size, n_taps, dtype=torch.complex64, device=device) * tap_std
                for sym in range(n_ofdm_sym):
                    # Simple Doppler: random walk phase per tap (approximating Jakes for low fd*T)
                    # For fd=5Hz and slot duration ~0.5ms, fd*T=0.0025, very slow fading
                    # Just use same taps for all symbols (quasi-static block fading)
                    taps_t = taps_c
                    h_freq[:, :, sym, rx, tx] = (taps_t.unsqueeze(1) * phase_freq_t.unsqueeze(0)).sum(dim=-1)

        # Build transmit grid
        tx_grid = torch.zeros(batch_size, n_tx, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)
        re_idx = 0
        for sym in pdsch_symbols:
            for sc in range(n_sc):
                for layer in range(n_layers):
                    tx_grid[:, layer, sym, sc] = x[:, re_idx, layer]
                re_idx += 1

        # Apply channel in frequency domain: y(rx, sym, sc) = sum_tx H(rx, tx, sym, sc) * x(tx, sym, sc)
        rx_grid = torch.zeros(batch_size, n_rx, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)
        for rx in range(n_rx):
            for tx in range(n_tx):
                rx_grid[:, rx, :, :] += tx_grid[:, tx, :, :] * h_freq[:, :, :, rx, tx].permute(0, 2, 1)

        # Add noise: C++ style - noise_var = 1/SNR_lin
        esn0_lin = 10.0 ** (esn0_db / 10.0)
        noise_var = 1.0 / esn0_lin
        sigma_dim = np.sqrt(noise_var / 2.0)
        rx_grid = rx_grid + torch.randn_like(rx_grid) * sigma_dim + 1j * torch.randn_like(rx_grid) * sigma_dim

        # Channel estimation
        if perfect_csi:
            h_est = h_freq.clone()
            no_est = torch.tensor(noise_var, dtype=torch.float32, device=device).expand(batch_size)
        else:
            h_est = torch.zeros_like(h_freq)
            for rx in range(n_rx):
                for tx in range(n_tx):
                    for d_sym in dmrs_symbols:
                        h_est[:, :, d_sym, rx, tx] = h_freq[:, :, d_sym, rx, tx]
                    for sym in pdsch_symbols:
                        if sym < dmrs_symbols[0]:
                            h_est[:, :, sym, rx, tx] = h_est[:, :, dmrs_symbols[0], rx, tx]
                        elif sym > dmrs_symbols[-1]:
                            h_est[:, :, sym, rx, tx] = h_est[:, :, dmrs_symbols[-1], rx, tx]
                        else:
                            d0, d1 = dmrs_symbols[0], dmrs_symbols[1]
                            w = (sym - d0) / (d1 - d0)
                            h_est[:, :, sym, rx, tx] = (1 - w) * h_est[:, :, d0, rx, tx] + w * h_est[:, :, d1, rx, tx]
            no_est = torch.tensor(noise_var, dtype=torch.float32, device=device).expand(batch_size)

        eq_list = []
        nv_list = []
        for sym in pdsch_symbols:
            H_sym = h_est[:, :, sym, :, :]
            y_sym = rx_grid[:, :, sym, :]
            for sc in range(n_sc):
                H_sc = H_sym[:, sc, :, :n_layers]
                y_sc = y_sym[:, :, sc]
                x_eq, eff_nv = lmmse_equalize(H_sc, y_sc, no_est, n_layers)
                for l in range(n_layers):
                    eq_list.append(x_eq[:, l])
                    nv_list.append(eff_nv[:, l])

        eq_syms = torch.stack(eq_list, dim=1)
        nv = torch.stack(nv_list, dim=1)

        llr_cat = demapper(eq_syms, nv)
        llr_cat = llr_cat[:, :G]

        b_hat, tb_crc_status = tb_dec(llr_cat)
        errors = (tb_crc_status < 0.5).float()
        return errors.sum()

    if Qm == 2:
        sinr_points = np.arange(-6.0, 18.0, 2.0)
    elif Qm == 4:
        sinr_points = np.arange(2.0, 26.0, 2.0)
    elif Qm == 6:
        sinr_points = np.arange(8.0, 32.0, 2.0)
    else:
        sinr_points = np.arange(14.0, 38.0, 2.0)

    batch_size = 10
    max_blocks = 50
    target_errors = 10

    method_name = "sionna_ideal" if perfect_csi else "sionna_ls"
    if out_file is None:
        os.makedirs("results/tdl_crosscheck", exist_ok=True)
        out_file = f"results/tdl_crosscheck/{method_name}_mcs{mcs}_table{mcs_table}_30prb_4x4.csv"

    with open(out_file, "w") as f:
        f.write("SINR_dB,Blocks,Errors,BLER\n")

    print(f"\n{'SINR(dB)':>10} {'Blocks':>10} {'Errors':>10} {'BLER':>12}")
    print("-" * 46)

    for sinr_db in sinr_points:
        n_blocks = 0
        n_errors = 0
        t0 = time.time()
        while n_blocks < max_blocks and n_errors < target_errors:
            bs = min(batch_size, max_blocks - n_blocks)
            err = run_block(bs, float(sinr_db))
            n_errors += int(err.item())
            n_blocks += bs
            if n_errors == n_blocks and n_blocks >= 10:
                break
        bler_val = n_errors / n_blocks if n_blocks > 0 else 1.0
        dt = time.time() - t0
        print(f"{sinr_db:10.2f} {n_blocks:10d} {n_errors:10d} {bler_val:12.4f} ({dt:.1f}s)")
        with open(out_file, "a") as f:
            f.write(f"{sinr_db:.2f},{n_blocks},{n_errors},{bler_val:.6f}\n")
        sys.stdout.flush()
        if bler_val < 0.02 and n_blocks >= 20:
            break

    print(f"\nResults saved to {out_file}")


if __name__ == "__main__":
    main()
