"""
Sionna 4x4 MIMO TDL-A BLER simulation
Supports LS channel estimation (linear interpolation) and Perfect CSI
Time-varying channel with Doppler
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

def get_mcs_params(mcs):
    if mcs == 3:
        return 2, 308, 2, 24, 96, 256, 340
    elif mcs == 10:
        return 4, 680, 4, 24, 480, 1024, 680
    elif mcs == 17:
        return 6, 1304, 6, 24, 840, 1024, 680
    elif mcs == 24:
        return 8, 1992, 8, 24, 1120, 1024, 680
    else:
        raise ValueError(f"Unsupported MCS {mcs}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python sionna_tdl_mimo.py <mcs> <perfect_csi=0/1> [out_file]")
        sys.exit(1)
    
    mcs = int(sys.argv[1])
    perfect_csi = int(sys.argv[2]) == 1 if len(sys.argv) > 2 else False
    out_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    Qm, k_info, n_layers, crc_len, n_re_per_prb, n_prb_test, n_info_test = get_mcs_params(mcs)
    
    n_prb = 25
    n_tx = n_layers
    n_rx = 4
    n_sc = n_prb * 12
    n_ofdm_sym = 14
    
    scs = 30e3
    fft_size = 1024
    delay_spread = 100e-9
    max_doppler = 70.0
    
    dmrs_type1_pos_add1 = [2, 11]
    
    crc_enc = CRCEncoder("CRC24A")
    crc_dec = CRCDecoder(crc_enc)
    
    k_crc = k_info + crc_len
    n_coded_per_rb = n_re_per_prb * Qm * n_layers
    
    n_pdsch_re_total = 0
    for sym in range(n_ofdm_sym):
        if sym not in dmrs_type1_pos_add1:
            n_pdsch_re_total += n_sc
    n_coded_total = n_pdsch_re_total * Qm
    n_coded = ((n_coded_total // (Qm * n_layers)) * Qm * n_layers)
    n_symbols_per_layer = n_coded // (Qm * n_layers)
    n_coded = n_symbols_per_layer * Qm * n_layers
    
    R = k_crc / n_coded
    print(f"MCS{mcs}: Qm={Qm}, k_info={k_info}, n_coded={n_coded}, R={R:.3f}")
    print(f"Config: {n_prb}PRB, {n_tx}Tx{n_rx}Rx, {n_layers}layers, TDL-A DS={delay_spread*1e9:.0f}ns, fd={max_doppler}Hz")
    print(f"Channel Est: {'Perfect/Ideal' if perfect_csi else 'LS + linear interp'}")
    
    device = "cpu"
    ldpc_enc = LDPC5GEncoder(k_crc, n_coded, num_bits_per_symbol=Qm)
    ldpc_dec = LDPC5GDecoder(ldpc_enc, num_iter=20, hard_out=True, return_info=False)
    constellation = Constellation("qam", Qm)
    mapper = Mapper(constellation=constellation)
    demapper = Demapper("app", constellation=constellation)
    
    half_sc = n_sc // 2
    freq_axis = []
    for sc in range(n_sc):
        if sc < half_sc:
            freq = (sc - half_sc) * scs
        else:
            freq = (sc - half_sc + 1) * scs
        freq_axis.append(freq)
    freq_axis = np.array(freq_axis, dtype=np.float64)
    
    taps_powers_db = np.array([
        -13.4, 0.0, -2.2, -4.0, -6.0, -8.2, -9.9, -10.5,
        -7.5, -15.9, -6.6, -16.7, -12.4, -15.2, -10.8,
        -11.3, -12.7, -16.2, -18.3, -18.9, -16.6, -19.9, -29.7])
    taps_norm_delays = np.array([
        0.0, 0.3819, 0.4025, 0.5868, 0.4610, 0.5375, 0.6708, 0.5750,
        0.7618, 1.5375, 1.8978, 2.2242, 2.1718, 2.4942, 2.5119,
        3.0582, 4.0810, 4.4579, 4.5695, 4.7966, 5.0066, 5.3043, 9.6586])
    taps_powers_lin = 10.0 ** (taps_powers_db / 10.0)
    pwr_sum = taps_powers_lin.sum()
    taps_powers_norm = taps_powers_lin / pwr_sum
    
    taps_powers_norm_t = torch.tensor(taps_powers_norm, dtype=torch.float32, device=device)
    taps_norm_delays_t = torch.tensor(taps_norm_delays, dtype=torch.float64, device=device)
    freq_axis_t = torch.tensor(freq_axis, dtype=torch.float64, device=device)
    
    delays = taps_norm_delays_t * delay_spread
    phase_freq = -2.0 * np.pi * freq_axis_t[:, None] * delays[None, :]
    phase_freq_t = torch.complex(torch.cos(phase_freq).float(), torch.sin(phase_freq).float())
    
    n_taps = len(taps_powers_norm)
    tap_std = torch.sqrt(taps_powers_norm_t / 2.0).float()
    
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
    dmrs_time_diffs = sym_times[dmrs_type1_pos_add1[1]] - sym_times[dmrs_type1_pos_add1[0]]
    
    dmrs_sc_type1 = []
    for sc_start in [0, 1]:
        scs_list = []
        sc = sc_start
        while sc < n_sc:
            scs_list.append(sc)
            sc += 2
        dmrs_sc_type1.append(scs_list)
    
    def run_block(batch_size, esn0_db):
        no = 10.0 ** (-esn0_db / 10.0)
        
        b = torch.randint(0, 2, (batch_size, k_info), dtype=torch.float32, device=device)
        b_crc = crc_enc(b)
        c = ldpc_enc(b_crc)
        x = mapper(c)
        x = x.reshape(batch_size, n_layers, n_symbols_per_layer)
        
        n_total = batch_size * n_ofdm_sym * n_rx * n_tx
        taps_re = torch.randn(n_total, n_taps, device=device) * tap_std
        taps_im = torch.randn(n_total, n_taps, device=device) * tap_std
        taps = torch.complex(taps_re, taps_im)
        taps = taps.reshape(batch_size, n_ofdm_sym, n_rx, n_tx, n_taps)
        
        aoa = torch.rand(batch_size, n_ofdm_sym, n_rx, n_tx, n_taps, device=device) * 2 * np.pi
        doppler_phase = 2.0 * np.pi * max_doppler * torch.tensor(sym_times, device=device).float()[None, :, None, None, None] * torch.cos(aoa)
        doppler_rot = torch.complex(torch.cos(doppler_phase), torch.sin(doppler_phase))
        taps = taps * doppler_rot
        
        h_freq = torch.einsum('brttp,fp->brtfs', taps, phase_freq_t)
        h_perfect = h_freq.permute(0, 2, 3, 4, 1)
        
        rx_grid = torch.zeros(batch_size, n_rx, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)
        tx_grid = torch.zeros(batch_size, n_tx, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)
        
        pdsch_sym_idx = 0
        pdsch_sc_idx = 0
        for sym in range(n_ofdm_sym):
            if sym in dmrs_type1_pos_add1:
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
                        if pdsch_sym_idx < n_symbols_per_layer:
                            tx_grid[:, layer, sym, sc] = x[:, layer, pdsch_sym_idx]
                            pdsch_sc_idx += 1
                            if pdsch_sc_idx >= n_sc:
                                pdsch_sc_idx = 0
                                pdsch_sym_idx += 1
        
        for sym in range(n_ofdm_sym):
            for rx in range(n_rx):
                for tx in range(n_tx):
                    rx_grid[:, rx, sym, :] += h_perfect[:, tx, :, rx, sym] * tx_grid[:, tx, sym, :]
        
        noise_std = np.sqrt(no / 2.0)
        noise = (torch.randn_like(rx_grid, dtype=torch.float32) + 
                 1j*torch.randn_like(rx_grid, dtype=torch.float32)) * noise_std
        rx_grid = rx_grid + noise
        
        if perfect_csi:
            h_est = h_perfect
            no_est = no
        else:
            h_est = torch.zeros_like(h_perfect)
            for dmrs_sym in dmrs_type1_pos_add1:
                for port in range(n_tx):
                    for sc in dmrs_sc_type1[port % 2]:
                        ref_val = tx_grid[:, port, dmrs_sym, sc:sc+1]
                        if torch.abs(ref_val[0,0]) > 1e-10:
                            for rx in range(n_rx):
                                h_est[:, port, sc:sc+1, rx, dmrs_sym] = rx_grid[:, rx:rx+1, dmrs_sym, sc] / ref_val
            
            for dmrs_sym in dmrs_type1_pos_add1:
                for port in range(n_tx):
                    for rx in range(n_rx):
                        dmrs_scs = dmrs_sc_type1[port % 2]
                        for i in range(len(dmrs_scs)):
                            sc_curr = dmrs_scs[i]
                            if i == 0:
                                sc_next = dmrs_scs[i+1]
                                h_next = h_est[:, port, sc_next, rx, dmrs_sym]
                                slope = h_next - h_est[:, port, sc_curr, rx, dmrs_sym]
                                slope /= (sc_next - sc_curr)
                                for sc in range(0, sc_curr):
                                    h_est[:, port, sc, rx, dmrs_sym] = h_est[:, port, sc_curr, rx, dmrs_sym] - slope * (sc_curr - sc)
                            elif i == len(dmrs_scs) - 1:
                                sc_prev = dmrs_scs[i-1]
                                h_prev = h_est[:, port, sc_prev, rx, dmrs_sym]
                                slope = h_est[:, port, sc_curr, rx, dmrs_sym] - h_prev
                                slope /= (sc_curr - sc_prev)
                                for sc in range(sc_curr+1, n_sc):
                                    h_est[:, port, sc, rx, dmrs_sym] = h_est[:, port, sc_curr, rx, dmrs_sym] + slope * (sc - sc_curr)
                            else:
                                sc_prev = dmrs_scs[i-1]
                                sc_next = dmrs_scs[i+1]
                                h_prev = h_est[:, port, sc_prev, rx, dmrs_sym]
                                h_next = h_est[:, port, sc_next, rx, dmrs_sym]
                                for sc in range(sc_prev+1, sc_curr):
                                    alpha = (sc - sc_prev) / (sc_curr - sc_prev)
                                    h_est[:, port, sc, rx, dmrs_sym] = h_prev * (1-alpha) + h_est[:, port, sc_curr, rx, dmrs_sym] * alpha
                                for sc in range(sc_curr+1, sc_next):
                                    alpha = (sc - sc_curr) / (sc_next - sc_curr)
                                    h_est[:, port, sc, rx, dmrs_sym] = h_est[:, port, sc_curr, rx, dmrs_sym] * (1-alpha) + h_next * alpha
            
            cross_corr = torch.complex(torch.zeros(batch_size, n_tx, n_rx, device=device), torch.zeros(batch_size, n_tx, n_rx, device=device))
            ref_pow = torch.zeros(batch_size, n_tx, n_rx, device=device)
            s0 = dmrs_type1_pos_add1[0]
            s1 = dmrs_type1_pos_add1[1]
            for port in range(n_tx):
                for rx in range(n_rx):
                    cross_corr[:, port, rx] = (h_est[:, port, :, rx, s1] * torch.conj(h_est[:, port, :, rx, s0])).sum(dim=1)
                    ref_pow[:, port, rx] = (torch.abs(h_est[:, port, :, rx, s0])**2).sum(dim=1)
            
            est_doppler = torch.angle(cross_corr.mean(dim=(1,2))) / (2 * np.pi * dmrs_time_diffs)
            
            for port in range(n_tx):
                for rx in range(n_rx):
                    h_a = h_est[:, port, :, rx, s0]
                    h_b = h_est[:, port, :, rx, s1]
                    for sym in range(n_ofdm_sym):
                        t_sym = sym_times[sym]
                        if sym <= s0:
                            dt = t_sym - sym_times[s0]
                            phase = -2 * np.pi * est_doppler[:, None] * dt
                            comp = torch.complex(torch.cos(phase), torch.sin(phase))
                            h_est[:, port, :, rx, sym] = h_a * comp
                        elif sym >= s1:
                            dt = t_sym - sym_times[s1]
                            phase = -2 * np.pi * est_doppler[:, None] * dt
                            comp = torch.complex(torch.cos(phase), torch.sin(phase))
                            h_est[:, port, :, rx, sym] = h_b * comp
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
                            h_est[:, port, :, rx, sym] = h_a * comp_a * (1-alpha) + h_b * comp_b * alpha
            
            no_est = no * 2.0
        
        llr_all = []
        pdsch_sym_idx = 0
        pdsch_sc_idx = 0
        for sym in range(n_ofdm_sym):
            if sym in dmrs_type1_pos_add1:
                continue
            for sc in range(n_sc):
                y_sc = rx_grid[:, :, sym, sc]
                H_sc = h_est[:, :, sc, :, sym].permute(0, 2, 1)
                
                H_H = torch.conj(H_sc.permute(0, 2, 1))
                Hty = torch.bmm(H_H, y_sc.unsqueeze(2)).squeeze(2)
                HtH = torch.bmm(H_H, H_sc)
                
                for layer in range(n_layers):
                    sig = Hty[:, layer]
                    var = HtH[:, layer, layer]
                    if torch.abs(var[0]) > 1e-10:
                        y_eq = sig / var
                        no_eq = no_est / (torch.abs(var)**2 + 1e-12) * (torch.abs(HtH[:, layer, :])**2).sum(dim=1)
                        llr = demapper(y_eq, no_eq.real.float())
                        llr_all.append(llr)
                    else:
                        llr_all.append(torch.zeros(batch_size, Qm, device=device))
                
                pdsch_sc_idx += 1
                if pdsch_sc_idx >= n_sc:
                    pdsch_sc_idx = 0
                    pdsch_sym_idx += 1
        
        if len(llr_all) * Qm > n_coded:
            llr_all = llr_all[:n_coded//Qm]
        llr_cat = torch.stack(llr_all, dim=2).reshape(batch_size, -1)[:, :n_coded]
        c_hat = ldpc_dec(llr_cat)
        _, crc_valid = crc_dec(c_hat)
        errors = (crc_valid < 0.5).float()
        return errors.sum()
    
    batch_size = 20
    if mcs == 3:
        sinr_points = np.arange(-5.0, 15.0, 1.0)
    elif mcs == 10:
        sinr_points = np.arange(5.0, 25.0, 1.0)
    elif mcs == 17:
        sinr_points = np.arange(10.0, 30.0, 1.0)
    else:
        sinr_points = np.arange(15.0, 35.0, 1.0)
    
    max_blocks = 200
    target_errors = 50
    
    method_name = "sionna_perfect" if perfect_csi else "sionna_ls"
    if out_file is None:
        os.makedirs("results/tdl_mimo", exist_ok=True)
        out_file = f"results/tdl_mimo/{method_name}_mcs{mcs}.csv"
    
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
        if bler_val < 0.01 and n_errors >= 10:
            break
    
    print(f"\nResults saved to {out_file}")

if __name__ == "__main__":
    main()
