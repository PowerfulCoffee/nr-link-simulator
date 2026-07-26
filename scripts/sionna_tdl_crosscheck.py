"""
Sionna 4x4 MIMO TDL-A BLER simulation for cross-check with C++ nr-link-simulator
Config: 30PRB, 4Rx, rank4, DMRS Type1 add-pos=1 (symbols 2&7), TDL-A DS=100ns
Supports both Ideal CSI and LS channel estimation with LMMSE equalization (whitening + Cholesky)
MCS Table 2 for 256QAM support
"""
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
import time
import sys

from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.mapping import Constellation, Mapper, Demapper


MCS_TABLE2_PDSCH = [
    (2, 120), (2, 193), (2, 308), (2, 449), (2, 602),
    (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658),
    (6, 466), (6, 517), (6, 567), (6, 616), (6, 666), (6, 719), (6, 772), (6, 822), (6, 873),
    (8, 682.5), (8, 711), (8, 754), (8, 797), (8, 841), (8, 885), (8, 916.5), (8, 948),
]


def get_mcs_params(mcs, n_prb, n_layers, dmrs_symbols, n_ofdm_sym=14, n_sc_per_rb=12):
    Qm, R_x1024 = MCS_TABLE2_PDSCH[mcs]
    target_coderate = R_x1024 / 1024.0

    n_sc = n_prb * n_sc_per_rb
    n_pdsch_sym = n_ofdm_sym - len(dmrs_symbols)
    n_re_per_prb = n_pdsch_sym * n_sc_per_rb

    n_re_tbs = min(156, n_re_per_prb) * n_prb
    tb_size_approx = target_coderate * n_re_tbs * Qm * n_layers
    if tb_size_approx <= 3824:
        n_small = max(3, int(np.floor(np.log2(tb_size_approx)) - 6))
        n_info_q = max(24, int(2**n_small * np.floor(tb_size_approx / 2**n_small)))
    else:
        n_large = int(np.floor(np.log2(tb_size_approx - 24)) - 5)
        n_info_q = max(3840, int(2**n_large * round((tb_size_approx - 24) / 2**n_large)))

    if n_info_q <= 3824:
        num_cb = 1
    elif target_coderate <= 0.25:
        num_cb = int(np.ceil((n_info_q + 24) / 3816.0))
    elif n_info_q > 8424:
        num_cb = int(np.ceil((n_info_q + 24) / 8424.0))
    else:
        num_cb = 1

    crc_len = 24 if n_info_q > 3824 else 16
    k_info = n_info_q
    k_crc = k_info + crc_len

    n_dmrs_re_per_rb_per_sym = 6
    n_pdsch_re_total = 0
    for sym in range(n_ofdm_sym):
        if sym not in dmrs_symbols:
            n_pdsch_re_total += n_sc
    n_coded_total = n_pdsch_re_total * Qm * n_layers
    G = n_coded_total
    G = ((G // (Qm * n_layers)) * (Qm * n_layers))
    n_symbols_per_layer = G // (Qm * n_layers)

    return Qm, target_coderate, k_info, crc_len, G, n_symbols_per_layer, n_re_per_prb


def solve_via_cholesky(A, B):
    L = torch.linalg.cholesky(A)
    Y = torch.linalg.solve_triangular(L, B, upper=False)
    X = torch.linalg.solve_triangular(L.conj().transpose(-1, -2), Y, upper=True)
    return X


def invert_via_cholesky(A):
    L = torch.linalg.cholesky(A)
    L_inv = torch.linalg.solve_triangular(L, torch.eye(A.shape[-1], device=A.device, dtype=A.dtype), upper=False)
    return L_inv.conj().transpose(-1, -2) @ L_inv


def lmmse_equalize(H, y, no, n_layers):
    """
    LMMSE equalizer matching C++ MmseEqualizer:
    - Whitening: H_tilde = H/sqrt(no), y_tilde = y/sqrt(no)
    - Solve: (H^H H + I) x_hat = H^H y_tilde
    - Bias compensation: x_hat/diag(I - A^{-1})
    - Returns: equalized symbols (batch, n_layers), effective noise var per layer (batch, n_layers)
    """
    sqrt_no = torch.sqrt(torch.clamp(no, min=1e-12))
    H_tilde = H / sqrt_no.unsqueeze(-1).unsqueeze(-1)
    y_tilde = y / sqrt_no.unsqueeze(-1)

    batch = H.shape[0]
    n_rx = H.shape[1]
    n_l = n_layers

    Hh = H_tilde.conj().transpose(-1, -2)
    HhH = Hh @ H_tilde
    A = HhH + torch.eye(n_l, device=H.device, dtype=H.dtype).unsqueeze(0)

    try:
        L, info = torch.linalg.cholesky_ex(A)
        if info.sum() > 0:
            A_reg = A.clone()
            for i in range(n_l):
                A_reg[:, i, i] += 1e-8
            L = torch.linalg.cholesky(A_reg)
            A = A_reg
        Hhy = Hh @ y_tilde.unsqueeze(-1)
        Y = torch.linalg.solve_triangular(L, Hhy, upper=False)
        x_hat = torch.linalg.solve_triangular(L.conj().transpose(-1, -2), Y, upper=True).squeeze(-1)

        L_inv = torch.linalg.solve_triangular(L, torch.eye(n_l, device=H.device, dtype=H.dtype).unsqueeze(0).expand(batch, -1, -1), upper=False)
        A_inv = L_inv.conj().transpose(-1, -2) @ L_inv

        eq = torch.zeros_like(x_hat)
        eff_noise = torch.zeros(batch, n_l, device=H.device, dtype=torch.float32)
        for l in range(n_l):
            err_var = A_inv[:, l, l].real.clamp(min=1e-12)
            d = (1.0 - err_var).clamp(min=1e-12, max=1.0 - 1e-12)
            eq[:, l] = x_hat[:, l] / d
            eff_noise[:, l] = (err_var / d).float()
    except Exception:
        eq = torch.zeros(batch, n_l, device=H.device, dtype=H.dtype)
        eff_noise = torch.ones(batch, n_l, device=H.device, dtype=torch.float32)

    return eq, eff_noise


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
    fft_size = 512
    delay_spread = 100e-9
    max_doppler = 5.0
    mcs_table = 2

    dmrs_symbols = [2, 7]

    Qm, target_coderate, k_info, crc_len, G, n_symbols_per_layer, n_re_per_prb = \
        get_mcs_params(mcs, n_prb, n_layers, dmrs_symbols)

    k_crc = k_info + crc_len
    n_coded = G

    mod_name = {2: "QPSK", 4: "16QAM", 6: "64QAM", 8: "256QAM"}[Qm]
    print(f"MCS{mcs} (table {mcs_table}): Qm={Qm} ({mod_name}), R={target_coderate:.4f}, k_info={k_info}, n_coded={n_coded}")
    print(f"Config: {n_prb}PRB, {n_tx}Tx{n_rx}Rx, {n_layers}layers, TDL-A DS={delay_spread*1e9:.0f}ns, fd={max_doppler}Hz")
    print(f"DMRS: Type1, symbols={dmrs_symbols} (add-pos=1)")
    print(f"Channel Est: {'Perfect/Ideal' if perfect_csi else 'LS + linear interp'}")

    device = "cpu"
    ldpc_enc = LDPC5GEncoder(k_crc, n_coded, num_bits_per_symbol=Qm)
    ldpc_dec = LDPC5GDecoder(ldpc_enc, num_iter=20, hard_out=True, return_info=False)
    crc_enc = CRCEncoder("CRC24A")
    crc_dec = CRCDecoder(crc_enc)
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
        for s in range(sym_idx):
            is_first = (s == 0)
            cp = 144 * fft_size // 2048
            if is_first:
                cp += 16 * fft_size // 2048
            t += tu * (fft_size + cp) / fft_size
        is_first = (sym_idx == 0)
        cp = 144 * fft_size // 2048
        if is_first:
            cp += 16 * fft_size // 2048
        t += tu * (fft_size + cp // 2) / fft_size
        return t

    sym_times = np.array([get_ofdm_time(s) for s in range(n_ofdm_sym)])
    sym_times_t = torch.tensor(sym_times, dtype=torch.float32, device=device)

    dmrs_sc_type1 = []
    for sc_start in [0, 1]:
        sc_list = list(range(sc_start, n_sc, 2))
        dmrs_sc_type1.append(sc_list)

    def run_block(batch_size, esn0_db):
        no = 10.0 ** (-esn0_db / 10.0)
        no_t = torch.tensor(no, dtype=torch.float32, device=device)

        b = torch.randint(0, 2, (batch_size, k_info), dtype=torch.float32, device=device)
        b_crc = crc_enc(b)
        c = ldpc_enc(b_crc)
        x_flat = mapper(c)

        x_per_re = x_flat.reshape(batch_size, n_symbols_per_layer, n_layers)

        n_total = batch_size * n_ofdm_sym * n_rx * n_tx
        taps_re = torch.randn(n_total, n_taps, device=device) * tap_std
        taps_im = torch.randn(n_total, n_taps, device=device) * tap_std
        taps = torch.complex(taps_re, taps_im)
        taps = taps.reshape(batch_size, n_ofdm_sym, n_rx, n_tx, n_taps)

        if max_doppler > 0:
            aoa = torch.rand(batch_size, n_ofdm_sym, n_rx, n_tx, n_taps, device=device) * 2 * np.pi
            doppler_phase = 2.0 * np.pi * max_doppler * sym_times_t[None, :, None, None, None] * torch.cos(aoa)
            doppler_rot = torch.complex(torch.cos(doppler_phase), torch.sin(doppler_phase))
            taps = taps * doppler_rot

        h_freq = torch.einsum('brttp,fp->brtfs', taps.double(), phase_freq_t).cfloat()
        h_perfect = h_freq.permute(0, 2, 4, 3, 1)

        rx_grid = torch.zeros(batch_size, n_rx, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)
        tx_grid = torch.zeros(batch_size, n_tx, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)

        re_idx = 0
        for sym in range(n_ofdm_sym):
            if sym in dmrs_symbols:
                for port in range(n_tx):
                    for sc in dmrs_sc_type1[port % 2]:
                        seed = (sym * n_sc + sc + port * 1000) % 2**31
                        rng = torch.Generator(device=device)
                        rng.manual_seed(int(seed))
                        dmrs_val = torch.randn(1, generator=rng, device=device) + 1j*torch.randn(1, generator=rng, device=device)
                        dmrs_val = dmrs_val[0] * np.sqrt(0.5)
                        tx_grid[:, port, sym, sc] = dmrs_val
            else:
                for sc in range(n_sc):
                    for layer in range(n_layers):
                        if re_idx < n_symbols_per_layer:
                            tx_grid[:, layer, sym, sc] = x_per_re[:, re_idx, layer]
                    re_idx += 1

        for sym in range(n_ofdm_sym):
            for tx_ant in range(n_tx):
                for rx_ant in range(n_rx):
                    rx_grid[:, rx_ant, sym, :] += h_perfect[:, tx_ant, :, rx_ant, sym] * tx_grid[:, tx_ant, sym, :]

        noise_std = np.sqrt(no / 2.0)
        noise = (torch.randn_like(rx_grid, dtype=torch.float32) +
                 1j*torch.randn_like(rx_grid, dtype=torch.float32)) * noise_std
        rx_grid = rx_grid + noise

        if perfect_csi:
            h_est = h_perfect
            no_est = no_t
        else:
            h_est = torch.zeros_like(h_perfect)
            for dmrs_sym in dmrs_symbols:
                for port in range(n_tx):
                    for sc in dmrs_sc_type1[port % 2]:
                        ref_val = tx_grid[:, port, dmrs_sym, sc:sc+1]
                        if torch.abs(ref_val[0,0]) > 1e-10:
                            for rx_ant in range(n_rx):
                                h_est[:, port, sc:sc+1, rx_ant, dmrs_sym] = rx_grid[:, rx_ant:rx_ant+1, dmrs_sym, sc] / ref_val

            for dmrs_sym in dmrs_symbols:
                for port in range(n_tx):
                    for rx_ant in range(n_rx):
                        dmrs_scs = dmrs_sc_type1[port % 2]
                        for i in range(len(dmrs_scs)):
                            sc_curr = dmrs_scs[i]
                            if i == 0:
                                sc_next = dmrs_scs[i+1]
                                h_next = h_est[:, port, sc_next, rx_ant, dmrs_sym]
                                slope = h_next - h_est[:, port, sc_curr, rx_ant, dmrs_sym]
                                slope /= (sc_next - sc_curr)
                                for sc in range(0, sc_curr):
                                    h_est[:, port, sc, rx_ant, dmrs_sym] = h_est[:, port, sc_curr, rx_ant, dmrs_sym] - slope * (sc_curr - sc)
                            elif i == len(dmrs_scs) - 1:
                                sc_prev = dmrs_scs[i-1]
                                h_prev = h_est[:, port, sc_prev, rx_ant, dmrs_sym]
                                slope = h_est[:, port, sc_curr, rx_ant, dmrs_sym] - h_prev
                                slope /= (sc_curr - sc_prev)
                                for sc in range(sc_curr+1, n_sc):
                                    h_est[:, port, sc, rx_ant, dmrs_sym] = h_est[:, port, sc_curr, rx_ant, dmrs_sym] + slope * (sc - sc_curr)
                            else:
                                sc_prev = dmrs_scs[i-1]
                                sc_next = dmrs_scs[i+1]
                                h_prev = h_est[:, port, sc_prev, rx_ant, dmrs_sym]
                                h_next = h_est[:, port, sc_next, rx_ant, dmrs_sym]
                                for sc in range(sc_prev+1, sc_curr):
                                    alpha = (sc - sc_prev) / (sc_curr - sc_prev)
                                    h_est[:, port, sc, rx_ant, dmrs_sym] = h_prev * (1-alpha) + h_est[:, port, sc_curr, rx_ant, dmrs_sym] * alpha
                                for sc in range(sc_curr+1, sc_next):
                                    alpha = (sc - sc_curr) / (sc_next - sc_curr)
                                    h_est[:, port, sc, rx_ant, dmrs_sym] = h_est[:, port, sc_curr, rx_ant, dmrs_sym] * (1-alpha) + h_next * alpha

            if len(dmrs_symbols) >= 2 and max_doppler > 0:
                s0 = dmrs_symbols[0]
                s1 = dmrs_symbols[1]
                cross_corr = torch.complex(torch.zeros(batch_size, n_tx, n_rx, device=device),
                                           torch.zeros(batch_size, n_tx, n_rx, device=device))
                for port in range(n_tx):
                    for rx_ant in range(n_rx):
                        cross_corr[:, port, rx_ant] = (h_est[:, port, :, rx_ant, s1] * torch.conj(h_est[:, port, :, rx_ant, s0])).sum(dim=1)
                dmrs_time_diff = sym_times[s1] - sym_times[s0]
                est_doppler = torch.angle(cross_corr.mean(dim=(1,2))) / (2 * np.pi * dmrs_time_diff)

                for port in range(n_tx):
                    for rx_ant in range(n_rx):
                        h_a = h_est[:, port, :, rx_ant, s0]
                        h_b = h_est[:, port, :, rx_ant, s1]
                        for sym in range(n_ofdm_sym):
                            t_sym = sym_times[sym]
                            if sym <= s0:
                                dt = t_sym - sym_times[s0]
                                phase = -2 * np.pi * est_doppler[:, None] * dt
                                comp = torch.complex(torch.cos(phase), torch.sin(phase))
                                h_est[:, port, :, rx_ant, sym] = h_a * comp
                            elif sym >= s1:
                                dt = t_sym - sym_times[s1]
                                phase = -2 * np.pi * est_doppler[:, None] * dt
                                comp = torch.complex(torch.cos(phase), torch.sin(phase))
                                h_est[:, port, :, rx_ant, sym] = h_b * comp
                            else:
                                ta = sym_times[s0]
                                tb = sym_times[s1]
                                alpha = (t_sym - ta) / (tb - ta)
                                dta = t_sym - ta
                                dtb = t_sym - tb
                                phase_a = -2 * np.pi * est_doppler[:, None] * dta
                                phase_b = -2 * np.pi * est_doppler[:, None] * dtb
                                comp_a = torch.complex(torch.cos(phase_a), torch.sin(phase_a))
                                comp_b = torch.complex(torch.cos(phase_b), torch.sin(phase_b))
                                h_est[:, port, :, rx_ant, sym] = h_a * comp_a * (1-alpha) + h_b * comp_b * alpha

            no_est = no_t * 2.0

        llr_all = []
        noise_vars_all = []
        pdsch_sym_idx = 0
        pdsch_sc_idx = 0

        for sym in range(n_ofdm_sym):
            if sym in dmrs_symbols:
                continue
            H_sym = h_est[:, :, :, :, sym].permute(0, 2, 3, 1)
            y_sym = rx_grid[:, :, sym, :]

            for sc in range(n_sc):
                H_sc = H_sym[:, sc, :, :]
                y_sc = y_sym[:, :, sc]

                x_eq, eff_nv = lmmse_equalize(H_sc, y_sc, no_est.expand(batch_size), n_layers)

                for l in range(n_layers):
                    llr_all.append(x_eq[:, l])
                    noise_vars_all.append(eff_nv[:, l])

                pdsch_sc_idx += 1
                if pdsch_sc_idx >= n_sc:
                    pdsch_sc_idx = 0
                    pdsch_sym_idx += 1

        n_needed = n_coded // Qm
        if len(llr_all) > n_needed:
            llr_all = llr_all[:n_needed]
            noise_vars_all = noise_vars_all[:n_needed]

        eq_syms = torch.stack(llr_all, dim=1)
        nv = torch.stack(noise_vars_all, dim=1)

        llr_cat = demapper(eq_syms, nv)
        llr_cat = llr_cat[:, :n_coded]

        c_hat = ldpc_dec(llr_cat)
        _, crc_valid = crc_dec(c_hat)
        errors = (crc_valid < 0.5).float()
        return errors.sum()

    if Qm == 2:
        sinr_points = np.arange(-6.0, 18.0, 1.0)
    elif Qm == 4:
        sinr_points = np.arange(2.0, 26.0, 1.0)
    elif Qm == 6:
        sinr_points = np.arange(8.0, 32.0, 1.0)
    else:
        sinr_points = np.arange(14.0, 38.0, 1.0)

    batch_size = 20
    max_blocks = 200
    target_errors = 50

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
        bler_val = n_errors / n_blocks if n_blocks > 0 else 1.0
        dt = time.time() - t0
        print(f"{sinr_db:10.2f} {n_blocks:10d} {n_errors:10d} {bler_val:12.4f} ({dt:.1f}s)")
        with open(out_file, "a") as f:
            f.write(f"{sinr_db:.2f},{n_blocks},{n_errors},{bler_val:.6f}\n")
        sys.stdout.flush()
        if bler_val < 0.01 and n_errors >= 10:
            break

    print(f"\nResults saved to {out_file}")


if __name__ == "__main__":
    main()
