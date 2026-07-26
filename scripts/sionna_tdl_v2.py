"""
Sionna 3PRB SISO TDL-A BLER with Doppler - simplified flat-RE approach
Based exactly on the verified sionna_mcs27_tdla.py pattern
For dual-DMRS (dmrs_add_pos=1), we flatten only non-DMRS REs sequentially.
LS receiver uses simple linear interpolation between two DMRS symbols.
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

def main():
    if len(sys.argv) < 2:
        print("Usage: python sionna_tdl_v2.py <mcs> <perfect_csi=0/1> [out_file]")
        sys.exit(1)
    
    mcs = int(sys.argv[1])
    perfect_csi = int(sys.argv[2]) == 1 if len(sys.argv) > 2 else False
    out_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    mcs_params = {
        3:  (2, 208,  864),
        10: (4, 576,  1728),
        17: (6, 1128, 2592),
        27: (6, 2280, 2592),
    }
    
    Qm, k_info, n_coded = mcs_params[mcs]
    crc_len = 24
    k_crc = k_info + crc_len
    n_symbols = n_coded // Qm
    
    n_prb = 3
    n_sc = n_prb * 12
    n_ofdm_sym = 14
    scs = 30e3
    delay_spread = 100e-9
    max_doppler = 70.0
    
    dmrs_pos = [2, 11]
    n_dmrs = len(dmrs_pos)
    n_pdsch_sym = n_ofdm_sym - n_dmrs
    n_pdsch_res = n_pdsch_sym * n_sc
    
    print(f"MCS{mcs}: Qm={Qm}, k_info={k_info}, n_coded={n_coded}, R={k_crc/n_coded:.3f}")
    print(f"n_pdsch_res={n_pdsch_res}, n_symbols={n_symbols}")
    
    device = "cpu"
    crc_enc = CRCEncoder("CRC24A")
    crc_dec = CRCDecoder(crc_enc)
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
    phase = -2.0 * np.pi * freq_axis_t[:, None] * delays[None, :]
    phase_t = torch.complex(torch.cos(phase).float(), torch.sin(phase).float())
    tap_std = torch.sqrt(taps_powers_norm_t / 2.0).float()
    
    pdsch_mask = np.ones(n_ofdm_sym * n_sc, dtype=bool)
    for dmrs_sym in dmrs_pos:
        for sc in range(n_sc):
            pdsch_mask[dmrs_sym * n_sc + sc] = False
    
    def get_ofdm_time(sym_idx):
        tu = 1.0 / scs
        t = 0.0
        fft_size = 64
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
        t += tu * (fft_size + cp//2) / fft_size
        return t
    
    sym_times = np.array([get_ofdm_time(s) for s in range(n_ofdm_sym)])
    
    dmrs_re_indices = {}
    for dmrs_sym in dmrs_pos:
        indices = []
        for sc in range(0, n_sc, 2):
            indices.append(dmrs_sym * n_sc + sc)
        dmrs_re_indices[dmrs_sym] = indices
    
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
        
        h_per_sc = torch.einsum('btp,fp->btf', taps, phase_t)
        h_flat = h_per_sc.reshape(batch_size, n_ofdm_sym * n_sc)
        
        tx_grid_flat = torch.zeros(batch_size, n_ofdm_sym * n_sc, dtype=torch.complex64, device=device)
        
        for dmrs_sym in dmrs_pos:
            for idx in dmrs_re_indices[dmrs_sym]:
                tx_grid_flat[:, idx] = 1.0
        
        pdsch_indices = np.where(pdsch_mask)[0]
        tx_grid_flat[:, pdsch_indices] = x
        
        y = h_flat * tx_grid_flat
        noise_std = np.sqrt(no / 2.0)
        noise = (torch.randn_like(y, dtype=torch.float32) + 1j*torch.randn_like(y, dtype=torch.float32)) * noise_std
        y = y + noise
        
        if perfect_csi:
            h = h_flat[:, pdsch_indices]
            y_eq = y[:, pdsch_indices] / h
            no_per_sc = no / (torch.abs(h)**2 + 1e-12)
        else:
            h_est_flat = torch.zeros_like(h_flat)
            for dmrs_sym in dmrs_pos:
                for idx in dmrs_re_indices[dmrs_sym]:
                    h_est_flat[:, idx] = y[:, idx] / (tx_grid_flat[:, idx] + 1e-12)
            
            for dmrs_sym in dmrs_pos:
                dmrs_scs = range(0, n_sc, 2)
                for i in range(len(dmrs_scs)):
                    sc_curr = dmrs_scs[i]
                    idx_curr = dmrs_sym * n_sc + sc_curr
                    if i == 0:
                        for sc in range(0, sc_curr):
                            h_est_flat[:, dmrs_sym * n_sc + sc] = h_est_flat[:, idx_curr]
                    if i + 1 < len(dmrs_scs):
                        sc_next = dmrs_scs[i+1]
                        idx_next = dmrs_sym * n_sc + sc_next
                        for sc in range(sc_curr+1, sc_next):
                            alpha = (sc - sc_curr) / (sc_next - sc_curr)
                            h_est_flat[:, dmrs_sym*n_sc + sc] = h_est_flat[:, idx_curr]*(1-alpha) + h_est_flat[:, idx_next]*alpha
                last_sc = dmrs_scs[-1]
                for sc in range(last_sc+1, n_sc):
                    h_est_flat[:, dmrs_sym*n_sc + sc] = h_est_flat[:, dmrs_sym*n_sc + last_sc]
            
            s0, s1 = dmrs_pos[0], dmrs_pos[1]
            t0 = sym_times[s0]
            t1 = sym_times[s1]
            h_est_pdsch = torch.zeros(batch_size, n_pdsch_res, dtype=torch.complex64, device=device)
            no_per_sc = torch.ones(batch_size, n_pdsch_res, dtype=torch.float32, device=device) * no
            
            p_idx = 0
            for sym in range(n_ofdm_sym):
                if sym in dmrs_pos:
                    continue
                t_sym = sym_times[sym]
                if sym <= s0:
                    alpha = 0.0
                elif sym >= s1:
                    alpha = 1.0
                else:
                    alpha = (t_sym - t0) / (t1 - t0)
                for sc in range(n_sc):
                    idx = sym * n_sc + sc
                    h0 = h_est_flat[:, s0 * n_sc + sc]
                    h1 = h_est_flat[:, s1 * n_sc + sc]
                    h_interp = h0 * (1-alpha) + h1 * alpha
                    h_est_pdsch[:, p_idx] = h_interp
                    no_per_sc[:, p_idx] = no / (torch.abs(h_interp)**2 + 1e-12)
                    p_idx += 1
            
            y_eq = y[:, pdsch_indices] / (h_est_pdsch + 1e-12)
        
        llr = demapper(y_eq, no_per_sc)
        llr = llr.reshape(batch_size, n_coded)
        c_hat = ldpc_dec(llr)
        _, crc_valid = crc_dec(c_hat)
        errors = (crc_valid < 0.5).float()
        return errors.sum()
    
    batch_size = 200
    if mcs == 3:
        sinr_points = np.arange(-6.0, 14.0, 1.0)
    elif mcs == 10:
        sinr_points = np.arange(2.0, 22.0, 1.0)
    elif mcs == 17:
        sinr_points = np.arange(8.0, 28.0, 1.0)
    else:
        sinr_points = np.arange(18.0, 38.0, 1.0)
    
    max_blocks = 1000
    target_errors = 100
    
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
    main()
