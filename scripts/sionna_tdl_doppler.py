"""
Sionna SISO TDL-A BLER simulation with Doppler
Supports LS channel estimation (linear interp w/o Doppler est) and Perfect CSI
Matches C++ simulator parameters exactly
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

MCS_PARAMS = {
    3:  {"Qm": 2, "k_info": 1800, "n_coded": 7200},
    10: {"Qm": 4, "k_info": 4736, "n_coded": 14400},
    17: {"Qm": 6, "k_info": 9224, "n_coded": 21600},
    24: {"Qm": 6, "k_info": 16392, "n_coded": 21600},
}

def main():
    if len(sys.argv) < 2:
        print("Usage: python sionna_tdl_doppler.py <mcs> <perfect_csi=0/1> [out_file]")
        sys.exit(1)
    
    mcs = int(sys.argv[1])
    perfect_csi = int(sys.argv[2]) == 1 if len(sys.argv) > 2 else False
    out_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    if mcs not in MCS_PARAMS:
        raise ValueError(f"Unsupported MCS {mcs}")
    
    p = MCS_PARAMS[mcs]
    Qm = p["Qm"]
    k_info = p["k_info"]
    n_coded = p["n_coded"]
    crc_len = 24
    k_crc = k_info + crc_len
    n_symbols = n_coded // Qm
    
    n_prb = 25
    n_sc = n_prb * 12
    n_ofdm_sym = 14
    
    scs = 30e3
    fft_size = 1024
    delay_spread = 100e-9
    max_doppler = 70.0
    
    dmrs_pos = [2, 11]
    
    R = k_crc / n_coded
    print(f"MCS{mcs}: Qm={Qm}, k_info={k_info}, n_coded={n_coded}, R={R:.3f}")
    print(f"Config: {n_prb}PRB SISO, TDL-A DS={delay_spread*1e9:.0f}ns, fd={max_doppler}Hz, DMRS add-pos=1")
    print(f"Channel Est: {'Perfect/Ideal' if perfect_csi else 'LS + linear interp (no Doppler est)'}")
    
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
    n_taps = len(taps_powers_norm)
    
    taps_powers_norm_t = torch.tensor(taps_powers_norm, dtype=torch.float32, device=device)
    taps_norm_delays_t = torch.tensor(taps_norm_delays, dtype=torch.float64, device=device)
    freq_axis_t = torch.tensor(freq_axis, dtype=torch.float64, device=device)
    
    delays = taps_norm_delays_t * delay_spread
    phase_freq = -2.0 * np.pi * freq_axis_t[:, None] * delays[None, :]
    phase_freq_t = torch.complex(torch.cos(phase_freq).float(), torch.sin(phase_freq).float())
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
    
    dmrs_scs = list(range(0, n_sc, 2))
    
    def run_block(batch_size, esn0_db):
        no = 10.0 ** (-esn0_db / 10.0)
        
        b = torch.randint(0, 2, (batch_size, k_info), dtype=torch.float32, device=device)
        b_crc = crc_enc(b)
        c = ldpc_enc(b_crc)
        x = mapper(c)
        x = x[:, :n_symbols]
        
        taps_re = torch.randn(batch_size, n_ofdm_sym, n_taps, device=device) * tap_std
        taps_im = torch.randn(batch_size, n_ofdm_sym, n_taps, device=device) * tap_std
        taps = torch.complex(taps_re, taps_im)
        
        aoa = torch.rand(batch_size, n_ofdm_sym, n_taps, device=device) * 2 * np.pi
        doppler_phase = 2.0 * np.pi * max_doppler * torch.tensor(sym_times, device=device).float()[None, :, None] * torch.cos(aoa)
        doppler_rot = torch.complex(torch.cos(doppler_phase), torch.sin(doppler_phase))
        taps = taps * doppler_rot
        
        h_freq = torch.einsum('btp,fp->btf', taps, phase_freq_t)
        
        rx_grid = torch.zeros(batch_size, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)
        tx_grid = torch.zeros(batch_size, n_ofdm_sym, n_sc, dtype=torch.complex64, device=device)
        
        pdsch_sym_idx = 0
        pdsch_sc_idx = 0
        
        for sym in range(n_ofdm_sym):
            if sym in dmrs_pos:
                for sc in dmrs_scs:
                    tx_grid[:, sym, sc] = 1.0
            else:
                for sc in range(n_sc):
                    if pdsch_sym_idx < n_symbols:
                        tx_grid[:, sym, sc] = x[:, pdsch_sym_idx]
                        pdsch_sc_idx += 1
                        if pdsch_sc_idx >= n_sc:
                            pdsch_sc_idx = 0
                            pdsch_sym_idx += 1
        
        for sym in range(n_ofdm_sym):
            rx_grid[:, sym, :] = h_freq[:, sym, :] * tx_grid[:, sym, :]
        
        noise_std = np.sqrt(no / 2.0)
        noise = (torch.randn_like(rx_grid, dtype=torch.float32) + 
                 1j*torch.randn_like(rx_grid, dtype=torch.float32)) * noise_std
        rx_grid = rx_grid + noise
        
        if perfect_csi:
            h_est = h_freq
        else:
            h_est = torch.zeros_like(h_freq)
            for dmrs_sym in dmrs_pos:
                for sc in dmrs_scs:
                    h_est[:, dmrs_sym, sc] = rx_grid[:, dmrs_sym, sc] / (tx_grid[:, dmrs_sym, sc] + 1e-12)
            
            for dmrs_sym in dmrs_pos:
                h_est[:, dmrs_sym, 0] = h_est[:, dmrs_sym, dmrs_scs[0]]
                for i in range(len(dmrs_scs)):
                    sc_curr = dmrs_scs[i]
                    if i + 1 < len(dmrs_scs):
                        sc_next = dmrs_scs[i+1]
                        h_curr = h_est[:, dmrs_sym, sc_curr]
                        h_next = h_est[:, dmrs_sym, sc_next]
                        for sc in range(sc_curr + 1, sc_next):
                            alpha = (sc - sc_curr) / (sc_next - sc_curr)
                            h_est[:, dmrs_sym, sc] = h_curr * (1 - alpha) + h_next * alpha
                    if i > 0:
                        sc_prev = dmrs_scs[i-1]
                        h_prev = h_est[:, dmrs_sym, sc_prev]
                        h_curr = h_est[:, dmrs_sym, sc_curr]
                        for sc in range(sc_prev + 1, sc_curr):
                            alpha = (sc - sc_prev) / (sc_curr - sc_prev)
                            h_est[:, dmrs_sym, sc] = h_prev * (1 - alpha) + h_curr * alpha
                last_dmrs = dmrs_scs[-1]
                for sc in range(last_dmrs + 1, n_sc):
                    h_est[:, dmrs_sym, sc] = h_est[:, dmrs_sym, last_dmrs]
            
            s0, s1 = dmrs_pos[0], dmrs_pos[1]
            for sym in range(n_ofdm_sym):
                t_sym = sym_times[sym]
                t0 = sym_times[s0]
                t1 = sym_times[s1]
                if sym <= s0:
                    alpha = 0.0
                elif sym >= s1:
                    alpha = 1.0
                else:
                    alpha = (t_sym - t0) / (t1 - t0)
                for sc in range(n_sc):
                    h_est[:, sym, sc] = h_est[:, s0, sc] * (1 - alpha) + h_est[:, s1, sc] * alpha
        
        llr_list = []
        pdsch_sym_idx = 0
        pdsch_sc_idx = 0
        for sym in range(n_ofdm_sym):
            if sym in dmrs_pos:
                continue
            for sc in range(n_sc):
                if pdsch_sym_idx >= n_symbols:
                    break
                y_eq = rx_grid[:, sym, sc] / (h_est[:, sym, sc] + 1e-12)
                no_eq = no / (torch.abs(h_est[:, sym, sc])**2 + 1e-12)
                llr = demapper(y_eq, no_eq.float())
                llr_list.append(llr)
                pdsch_sc_idx += 1
                if pdsch_sc_idx >= n_sc:
                    pdsch_sc_idx = 0
                    pdsch_sym_idx += 1
        
        llr_cat = torch.stack(llr_list, dim=1).reshape(batch_size, -1)[:, :n_coded]
        c_hat = ldpc_dec(llr_cat)
        _, crc_valid = crc_dec(c_hat)
        errors = (crc_valid < 0.5).float()
        return errors.sum()
    
    batch_size = 50
    if mcs == 3:
        sinr_points = np.arange(-6.0, 8.0, 1.0)
    elif mcs == 10:
        sinr_points = np.arange(2.0, 18.0, 1.0)
    elif mcs == 17:
        sinr_points = np.arange(8.0, 22.0, 1.0)
    else:
        sinr_points = np.arange(14.0, 28.0, 1.0)
    
    max_blocks = 300
    target_errors = 50
    
    method_name = "sionna_perfect" if perfect_csi else "sionna_ls_nodoppler"
    if out_file is None:
        os.makedirs("results/tdl_doppler", exist_ok=True)
        out_file = f"results/tdl_doppler/{method_name}_mcs{mcs}.csv"
    
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
    crc_enc = CRCEncoder("CRC24A")
    crc_dec = CRCDecoder(crc_enc)
    main()
