import torch
import numpy as np
from sionna.phy.nr import TBEncoder, TBDecoder, calculate_tb_size
from sionna.phy.mapping import Constellation, Mapper, Demapper
import argparse

MCS_TABLE_2 = [
    (2, 120), (2, 193), (2, 308), (2, 449), (2, 602),
    (4, 378), (4, 434), (4, 490), (4, 553), (4, 616), (4, 658),
    (6, 466), (6, 517), (6, 567), (6, 616), (6, 666), (6, 719),
    (6, 772), (6, 822), (6, 873),
    (8, 682), (8, 711), (8, 754), (8, 797), (8, 841), (8, 885), (8, 916), (8, 948),
]

TDL_A_DELAYS_NS = [0.0, 0.3819, 0.4025, 0.5868, 0.4610, 0.5375, 0.6708, 0.5750,
                   0.7618, 1.5375, 1.8978, 2.2242, 2.1718, 2.4942, 2.5119,
                   3.0582, 4.0810, 4.4579, 4.5695, 4.7966, 5.0066, 5.3043, 9.6586]
TDL_A_POWERS_DB = [-13.4, 0.0, -2.2, -4.0, -6.0, -8.2, -9.9, -10.5,
                   -7.5, -15.9, -6.6, -16.7, -12.4, -15.2, -10.8,
                   -11.3, -12.7, -16.2, -18.3, -18.9, -16.6, -19.9, -29.7]

def generate_pn_sequence(cinit, length):
    Nc = 1600
    N = Nc + length
    x1 = np.zeros(N + 31, dtype=np.uint8)
    x2 = np.zeros(N + 31, dtype=np.uint8)
    x1[0] = 1
    for i in range(31):
        x2[i] = (cinit >> i) & 1
    for n in range(N):
        x1[n + 31] = (x1[n + 3] + x1[n]) % 2
        x2[n + 31] = (x2[n + 3] + x2[n + 2] + x2[n + 1] + x2[n]) % 2
    seq = np.zeros(length, dtype=np.int8)
    for n in range(length):
        seq[n] = (x1[n + Nc] + x2[n + Nc]) % 2
    return seq

def get_dmrs_symbol(c_seq, m):
    c0 = c_seq[2 * m]
    c1 = c_seq[2 * m + 1]
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    re = (1 - 2 * c0) * inv_sqrt2
    im = (1 - 2 * c1) * inv_sqrt2
    return complex(re, im)

def build_dmrs_grid(n_sc, n_sym, n_layers, dmrs_symbols, slot_idx=0, n_symb_slot=14, N_id=0):
    grid = np.zeros((n_layers, n_sym, n_sc), dtype=np.complex64)
    dmrs_sc_per_layer = [list(range(l, n_sc, n_layers)) for l in range(n_layers)]
    for sym in dmrs_symbols:
        cinit = (1 << 17) * (n_symb_slot * slot_idx + sym + 1) * (2 * N_id + 1) + 2 * N_id
        c_seq = generate_pn_sequence(cinit, 2 * n_sc)
        for l in range(n_layers):
            for sc in dmrs_sc_per_layer[l]:
                grid[l, sym, sc] = get_dmrs_symbol(c_seq, sc)
    return torch.from_numpy(grid), dmrs_sc_per_layer

def get_fft_size(n_sc):
    fft_size = 1
    while fft_size < n_sc:
        fft_size <<= 1
    return fft_size

def get_cp_len(fft_size, is_first):
    cp = 144 * fft_size // 2048
    return cp + (16 * fft_size // 2048 if is_first else 0)

def get_ofdm_symbol_time(scs, fft_size, sym_idx):
    tu = 1.0 / scs
    t = 0.0
    for s in range(sym_idx):
        cp = get_cp_len(fft_size, s == 0)
        t += tu * (fft_size + cp) / fft_size
    cp = get_cp_len(fft_size, sym_idx == 0)
    t += tu * (fft_size + cp / 2) / fft_size
    return t

def generate_tdl_channel(bs, n_sc, n_sym, n_rx, n_tx, delay_spread_ns, max_doppler_hz, scs=15000, seed=42):
    n_taps = len(TDL_A_DELAYS_NS)
    fft_size = get_fft_size(n_sc)
    half_sc = n_sc // 2
    
    delays = np.array(TDL_A_DELAYS_NS, dtype=np.float64) * delay_spread_ns * 1e-9
    powers_lin = 10.0 ** (np.array(TDL_A_POWERS_DB, dtype=np.float64) / 10.0)
    pwr_sum = powers_lin.sum()
    
    rng = np.random.default_rng(seed)
    
    tap_scale = np.sqrt(powers_lin / (2.0 * pwr_sum))
    
    h_np = np.zeros((bs, n_rx, n_tx, n_sym, n_sc), dtype=np.complex64)
    
    sc_indices = np.arange(n_sc)
    freqs = np.where(sc_indices < half_sc,
                     (sc_indices - half_sc) * scs,
                     (sc_indices - half_sc + 1) * scs).astype(np.float64)
    
    sym_times = np.zeros(n_sym, dtype=np.float64)
    t_start = 0.0
    for sym in range(n_sym):
        is_first = (sym == 0)
        cp_len = get_cp_len(fft_size, is_first)
        tu = 1.0 / scs
        sym_dur = tu * (fft_size + cp_len) / fft_size
        sym_times[sym] = t_start + sym_dur / 2.0
        t_start += sym_dur
    
    for b in range(bs):
        tap_coeffs_init = (rng.standard_normal((n_rx, n_tx, n_taps)) + 
                          1j * rng.standard_normal((n_rx, n_tx, n_taps))) * tap_scale[np.newaxis, np.newaxis, :]
        
        if max_doppler_hz > 0:
            aoa = rng.uniform(0, 2 * np.pi, (n_rx, n_tx, n_taps))
            cos_aoa = np.cos(aoa)
            dop_phase_init = rng.uniform(0, 2 * np.pi, (n_rx, n_tx, n_taps))
        else:
            cos_aoa = np.zeros((n_rx, n_tx, n_taps))
            dop_phase_init = np.zeros((n_rx, n_tx, n_taps))
        
        for sym in range(n_sym):
            t_mid = sym_times[sym]
            
            fd = max_doppler_hz * cos_aoa
            dop_phase = 2.0 * np.pi * fd * t_mid + dop_phase_init
            tap_at_sym = tap_coeffs_init * np.exp(1j * dop_phase)
            
            phase_freq = -2.0 * np.pi * np.outer(delays, freqs)
            exp_phase = np.exp(1j * phase_freq)
            
            h_freq = tap_at_sym @ exp_phase
            h_np[b, :, :, sym, :] = h_freq.astype(np.complex64)
    
    h = torch.from_numpy(h_np).permute(0, 1, 3, 4, 2).contiguous()
    return h

def ls_channel_estimate(rx, dmrs_ref, dmrs_syms, dmrs_sc_per_layer, n_sc, n_sym, n_rx, n_layers, scs=15000):
    bs = rx.shape[0]
    device = rx.device
    h = torch.zeros(bs, n_rx, n_sym, n_sc, n_layers, dtype=torch.complex64, device=device)

    for sym in dmrs_syms:
        for rx_ant in range(n_rx):
            for l in range(n_layers):
                scl = dmrs_sc_per_layer[l]
                for i, sc in enumerate(scl):
                    ref = dmrs_ref[l, sym, sc]
                    if abs(ref) > 1e-10:
                        h[:, rx_ant, sym, sc, l] = rx[:, rx_ant, sym, sc] / ref
                    else:
                        h[:, rx_ant, sym, sc, l] = 1.0

        for rx_ant in range(n_rx):
            for l in range(n_layers):
                scl = dmrs_sc_per_layer[l]
                if len(scl) < 2:
                    continue
                first_sc = scl[0]
                last_sc = scl[-1]
                h_ant = h[:, rx_ant, sym, :, l]
                h_ant[:, :first_sc] = h_ant[:, first_sc:first_sc+1]
                for i in range(len(scl) - 1):
                    sa, sb = scl[i], scl[i+1]
                    if sb - sa > 1:
                        ha = h_ant[:, sa:sa+1]
                        hb = h_ant[:, sb:sb+1]
                        for off in range(1, sb - sa):
                            alpha = off / (sb - sa)
                            h_ant[:, sa+off] = (ha[:, 0] * (1 - alpha) + hb[:, 0] * alpha)
                h_ant[:, last_sc:] = h_ant[:, last_sc:last_sc+1]

    fft_size = get_fft_size(n_sc)
    
    if len(dmrs_syms) >= 2:
        sym_a, sym_b = dmrs_syms[0], dmrs_syms[1]
        t_a = get_ofdm_symbol_time(scs, fft_size, sym_a)
        t_b = get_ofdm_symbol_time(scs, fft_size, sym_b)
        cross_corr = torch.zeros(bs, dtype=torch.complex64, device=device)
        for rx_ant in range(n_rx):
            for l in range(n_layers):
                ha = h[:, rx_ant, sym_a, :, l]
                hb = h[:, rx_ant, sym_b, :, l]
                cross_corr += (hb * ha.conj()).sum(dim=-1)
        phase_diff = torch.angle(cross_corr)
        delta_t = t_b - t_a
        est_doppler = phase_diff / (2 * np.pi * delta_t)
    else:
        est_doppler = torch.zeros(bs, device=device)

    for rx_ant in range(n_rx):
        for l in range(n_layers):
            h_ant = h[:, rx_ant, :, :, l]
            first_dmrs = dmrs_syms[0]
            for sym in range(first_dmrs):
                h_ant[:, sym, :] = h_ant[:, first_dmrs, :]
            for i in range(len(dmrs_syms) - 1):
                sa, sb = dmrs_syms[i], dmrs_syms[i+1]
                dist = sb - sa
                if dist <= 0: continue
                for sym in range(sa + 1, sb):
                    alpha = (sym - sa) / dist
                    h_ant[:, sym, :] = h_ant[:, sa, :] * (1 - alpha) + h_ant[:, sb, :] * alpha
            last_dmrs = dmrs_syms[-1]
            for sym in range(last_dmrs + 1, n_sym):
                h_ant[:, sym, :] = h_ant[:, last_dmrs, :]
    return h

def mmse_equalize(rx, h_est, noise_var, n_layers):
    bs, n_rx, n_pdsch_sym, n_sc, _ = h_est.shape
    n_re = n_pdsch_sym * n_sc
    rx_flat = rx.permute(0, 2, 3, 1).reshape(bs, n_re, n_rx)
    H = h_est.permute(0, 2, 3, 1, 4).reshape(bs, n_re, n_rx, n_layers)
    sqrt_no = torch.sqrt(torch.tensor(noise_var, dtype=torch.float32))
    H_tilde = H / sqrt_no
    y_tilde = rx_flat / sqrt_no
    Hh = H_tilde.conj().transpose(-2, -1)
    HhH = Hh @ H_tilde
    eye = torch.eye(n_layers, dtype=torch.complex64, device=H.device).view(1, 1, n_layers, n_layers)
    A = HhH + eye
    Hhy = Hh @ y_tilde.unsqueeze(-1)
    
    A_inv = torch.linalg.inv(A)
    x_hat = torch.linalg.solve(A, Hhy).squeeze(-1)
    
    eff_noise = torch.zeros(bs, n_re, n_layers, dtype=torch.float32)
    x_eq = torch.zeros_like(x_hat)
    for l in range(n_layers):
        err_var = A_inv[..., l, l].real
        err_var = torch.clamp(err_var, min=1e-12)
        d = 1.0 - err_var
        d = torch.where(d.abs() < 1e-12, torch.full_like(d, 1e-12), d)
        x_eq[..., l] = x_hat[..., l] / d
        eff_noise[..., l] = err_var / d
    
    x_eq = x_eq.reshape(bs, n_pdsch_sym, n_sc, n_layers).permute(0, 3, 1, 2)
    eff_noise = eff_noise.reshape(bs, n_pdsch_sym, n_sc, n_layers).permute(0, 3, 1, 2)
    return x_eq, eff_noise

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mcs', type=int, default=4)
    parser.add_argument('--snr-start', type=float, default=0)
    parser.add_argument('--snr-end', type=float, default=20)
    parser.add_argument('--snr-step', type=float, default=1)
    parser.add_argument('--max-blocks', type=int, default=50)
    parser.add_argument('--target-errors', type=int, default=20)
    parser.add_argument('--batch-size', type=int, default=2)
    parser.add_argument('--ce-type', type=str, default='ideal', choices=['ideal','ls'])
    parser.add_argument('--delay-spread', type=float, default=15)
    parser.add_argument('--doppler', type=float, default=5)
    parser.add_argument('--output', type=str)
    args = parser.parse_args()

    qm, R_x1024 = MCS_TABLE_2[args.mcs]
    target_coderate = R_x1024 / 1024.0

    n_prb = 30
    n_sc = n_prb * 12
    n_ofdm_sym = 14
    n_layers = 4
    n_rx = 4
    n_tx = 4

    dmrs_symbols = [2, 7]
    pdsch_symbols = [s for s in range(n_ofdm_sym) if s not in dmrs_symbols]
    n_dmrs_per_prb = 12 * len(dmrs_symbols)

    tb_config = calculate_tb_size(
        modulation_order=qm,
        target_coderate=target_coderate,
        num_prbs=n_prb,
        num_ofdm_symbols=n_ofdm_sym,
        num_dmrs_per_prb=n_dmrs_per_prb,
        num_layers=n_layers,
        verbose=False
    )
    tb_size = int(tb_config[0])
    cb_sizes = tb_config[5]
    cw_length = int(torch.sum(cb_sizes))

    constellation = Constellation("qam", num_bits_per_symbol=qm, normalize=True)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)
    tb_enc = TBEncoder(
        tb_size, cw_length, target_coderate, qm,
        num_layers=n_layers, channel_type="PDSCH", verbose=False
    )
    tb_dec = TBDecoder(tb_enc, num_bp_iter=20)

    dmrs_ref_np, dmrs_sc_per_layer = build_dmrs_grid(n_sc, n_ofdm_sym, n_layers, dmrs_symbols)
    dmrs_ref = dmrs_ref_np.to(torch.complex64)

    print(f"MCS={args.mcs} Qm={qm} R={target_coderate:.3f} TB={tb_size} G={cw_length}", flush=True)
    print(f"SNR range: {args.snr_start} to {args.snr_end} dB, step={args.snr_step} dB, CE={args.ce_type}", flush=True)
    print(f"Channel: TDL-A, DS={args.delay_spread}ns, fd={args.doppler}Hz", flush=True)

    results = []
    snr = args.snr_start
    block_seed = 12345
    while snr <= args.snr_end + 1e-9:
        noise_var = 1.0 / (10 ** (snr / 10.0))
        n_errors = 0
        n_blocks = 0
        chan_seed = 42
        while n_blocks < args.max_blocks and n_errors < args.target_errors:
            bs = min(args.batch_size, args.max_blocks - n_blocks)
            bits = torch.randint(0, 2, (bs, tb_size))
            cw = tb_enc(bits)
            syms = mapper(cw)
            n_re_per_layer = len(pdsch_symbols) * n_sc
            syms = syms.reshape(bs, n_re_per_layer, n_layers)
            syms = syms.permute(0, 2, 1).reshape(bs, n_layers, len(pdsch_symbols), n_sc)

            tx_grid = torch.zeros(bs, n_tx, n_ofdm_sym, n_sc, dtype=torch.complex64)
            for i, s in enumerate(pdsch_symbols):
                tx_grid[:, :n_layers, s, :] = syms[:, :, i, :]
            for l in range(n_layers):
                for sym in dmrs_symbols:
                    for sc in dmrs_sc_per_layer[l]:
                        tx_grid[:, l, sym, sc] = dmrs_ref[l, sym, sc]

            h_t = generate_tdl_channel(bs, n_sc, n_ofdm_sym, n_rx, n_tx,
                                        args.delay_spread, args.doppler, scs=15000, seed=chan_seed)
            chan_seed += 1
            
            rx_signal = torch.einsum('brscl,blsc->brsc', h_t[:, :, :, :, :n_layers], tx_grid)
            noise = (torch.randn_like(rx_signal) + 1j * torch.randn_like(rx_signal)) * np.sqrt(noise_var / 2.0)
            rx_signal = rx_signal + noise

            if args.ce_type == 'ideal':
                h_est = h_t[:, :, :, :, :n_layers]
            else:
                h_est = ls_channel_estimate(rx_signal, dmrs_ref, dmrs_symbols,
                                            dmrs_sc_per_layer, n_sc, n_ofdm_sym, n_rx, n_layers)

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
            llr = llr.reshape(bs, cw_length)
            dec_out, crc_ok = tb_dec(llr)
            if isinstance(dec_out, (list, tuple)):
                dec_bits = dec_out[0].long()
            else:
                dec_bits = dec_out.long()
            for b in range(bs):
                n_blocks += 1
                if not crc_ok[b] or not torch.equal(bits[b], dec_bits[b]):
                    n_errors += 1

        bler = n_errors / n_blocks
        print(f"SNR={snr:.1f} dB: BLER={bler:.4f} ({n_errors}/{n_blocks})", flush=True)
        results.append((snr, n_blocks, n_errors, bler))
        snr += args.snr_step

    if args.output:
        with open(args.output, 'w') as f:
            f.write("SINR_dB,Blocks,Errors,BLER\n")
            for snr_v, blocks, errors, bler_v in results:
                f.write(f"{snr_v:.1f},{blocks},{errors},{bler_v:.6f}\n")
        print(f"Results saved to {args.output}")

if __name__ == '__main__':
    main()
