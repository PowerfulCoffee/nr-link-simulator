"""
Sionna MCS27 BLER over TDL-A (perfect CSI) - Es/N0
Matches C++ configuration: 3PRB, 15kHz SCS, MCS27 (64QAM, R=910/1024), TBS=2472
TDL-A 30ns DS, perfect CSI, LDPC 20 iterations
Uses Es/N0 directly (no = 10^(-EsN0_dB/10)) matching C++ simulator
PyTorch version for Sionna 2.x
"""
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import torch
import numpy as np
import time

from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.mapping import Constellation, Mapper, Demapper

def main():
    print("=== Sionna MCS27 BLER over TDL-A (perfect CSI, Es/N0) ===")

    k_info = 2472
    crc_len = 24
    k_crc = k_info + crc_len
    R = 910 / 1024
    n_coded = 2808
    Qm = 6
    n_symbols = n_coded // Qm
    assert n_coded == n_symbols * Qm

    n_prb = 3
    n_sc = n_prb * 12
    n_ofdm_sym = 14
    dmrs_sym = 2
    n_pdsch_res = (n_ofdm_sym - 1) * n_sc
    print(f"n_symbols={n_symbols}, n_pdsch_res={n_pdsch_res}")
    assert n_pdsch_res == n_symbols

    delay_spread = 30e-9
    fft_size = 64
    print(f"TDL-A DS={delay_spread*1e9:.0f}ns, Es/N0 sweep")

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
            freq = (sc - half_sc) * 15e3
        else:
            freq = (sc - half_sc + 1) * 15e3
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
    phase = -2.0 * np.pi * freq_axis_t[:, None] * delays[None, :]
    phase_t = torch.complex(torch.cos(phase).float(), torch.sin(phase).float())

    n_taps = len(taps_powers_norm)
    tap_std = torch.sqrt(taps_powers_norm_t / 2.0).float()

    def run_block(batch_size, esn0_db):
        no = 10.0 ** (-esn0_db / 10.0)

        b = torch.randint(0, 2, (batch_size, k_info), dtype=torch.float32, device=device)
        b_crc = crc_enc(b)
        c = ldpc_enc(b_crc)
        x = mapper(c)
        x = x.reshape(batch_size, n_symbols)

        n_total = batch_size * (n_ofdm_sym - 1)
        taps_re = torch.randn(n_total, n_taps, device=device) * tap_std
        taps_im = torch.randn(n_total, n_taps, device=device) * tap_std
        taps = torch.complex(taps_re, taps_im)
        taps = taps.reshape(batch_size, n_ofdm_sym - 1, n_taps)

        h_per_sc = torch.einsum('btp,fp->btf', taps, phase_t)
        h = h_per_sc.reshape(batch_size, -1)[:, :n_symbols]

        y = h * x
        noise_std = np.sqrt(no / 2.0)
        noise = (torch.randn_like(y, dtype=torch.float32) +
                 1j*torch.randn_like(y, dtype=torch.float32)) * noise_std
        y = y + noise

        h_amp = torch.abs(h).float()
        y_eq = y / h
        no_per_sc = no / (h_amp ** 2 + 1e-12)

        llr = demapper(y_eq, no_per_sc)
        llr = llr.reshape(batch_size, n_coded)
        c_hat = ldpc_dec(llr)
        _, crc_valid = crc_dec(c_hat)
        errors = (crc_valid < 0.5).float()
        return errors.sum()

    batch_size = 100
    sinr_points = np.arange(20.0, 42.0, 2.0)
    max_blocks = 1000
    target_errors = 100

    out_file = "sionna_mcs27_tdla_perfect_csi.csv"
    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", out_file)
    with open(out_path, "w") as f:
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
        with open(out_path, "a") as f:
            f.write(f"{sinr_db:.2f},{n_blocks},{n_errors},{bler_val:.6f}\n")
        if bler_val < 0.01 and n_errors > 0:
            break

    print(f"\nResults saved to {out_path}")

if __name__ == "__main__":
    main()
