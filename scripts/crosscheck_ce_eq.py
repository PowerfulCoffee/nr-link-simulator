#!/usr/bin/env python3
"""Cross-check C++ vs Sionna LS Channel Estimation and LMMSE Equalization.

Generates identical inputs, runs both implementations, and compares outputs.
C++ algorithms are faithfully reimplemented in Python based on:
- ChannelEstimator.cpp (LsChannelEstimator)
- Equalizer.cpp (MmseEqualizer, SISO fast path)
"""

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['CUDA_VISIBLE_DEVICES'] = ''

import sys
sys.path.insert(0, '/workspace/sionna/src')

import torch
torch.manual_seed(42)

import torch.compiler
_orig_compile = torch.compile
def _no_compile(*args, **kwargs):
    if len(args) > 0 and callable(args[0]):
        return args[0]
    return lambda f: f
torch.compile = _no_compile
for _attr in ['is_inductor_compiled', 'is_compiling', 'compile_pystate', 'codecache',
              'get_num_threads', 'set_num_threads', 'reset', 'is_torchdynamo_compiling']:
    if not hasattr(torch.compiler, _attr):
        setattr(torch.compiler, _attr, lambda *args, **kwargs: None)
if not hasattr(torch, 'is_dynamo_compiling'):
    torch.is_dynamo_compiling = lambda: False

import numpy as np
from sionna.phy.nr import PUSCHConfig, PUSCHTransmitter, PUSCHLSChannelEstimator
from sionna.phy.mimo import lmmse_equalizer
from sionna.phy.utils import flatten_last_dims

np.set_printoptions(precision=6, suppress=True)


def cpp_ls_channel_estimate(y, pilot_positions, pilot_values, n_sc, n_sym, n_rx_ant, n_layers):
    """Replicate C++ LsChannelEstimator::estimate() in Python.
    
    Args:
        y: numpy array shape [n_rx_ant, n_sym, n_sc] complex - received resource grid
        pilot_positions: list of (sym, sc) tuples - DMRS pilot locations
        pilot_values: numpy array shape [n_pilots] complex - known pilot values
        n_sc: number of subcarriers
        n_sym: number of OFDM symbols
        n_rx_ant: number of receive antennas
        n_layers: number of layers
    
    Returns:
        h_est: numpy array shape [n_sc, n_sym, n_rx_ant*n_layers] complex
    """
    n_ch = n_rx_ant * n_layers
    h_est = np.zeros((n_sc, n_sym, n_ch), dtype=complex)
    
    # Group pilot positions by symbol
    pilot_syms = sorted(set(sym for sym, sc in pilot_positions))
    if not pilot_syms:
        h_est.fill(1.0 + 0j)
        return h_est
    
    # Build a dict mapping (sym, sc) -> pilot index
    pilot_lookup = {}
    pilots_by_sym = {}
    for idx, (sym, sc) in enumerate(pilot_positions):
        pilot_lookup[(sym, sc)] = idx
        if sym not in pilots_by_sym:
            pilots_by_sym[sym] = []
        pilots_by_sym[sym].append(sc)
    
    for sym in pilots_by_sym:
        pilots_by_sym[sym].sort()
    
    ref_tx_ant = min(n_rx_ant, n_layers)
    
    for dmrs_sym in pilot_syms:
        for rx in range(n_rx_ant):
            for l in range(n_layers):
                ch_idx = rx * n_layers + l
                tx = l % ref_tx_ant
                
                dmrs_sc_list = pilots_by_sym[dmrs_sym]
                
                # Step 1: LS estimates at pilot positions
                for sc in dmrs_sc_list:
                    pidx = pilot_lookup[(dmrs_sym, sc)]
                    rx_val = y[rx, dmrs_sym, sc]
                    ref_val = pilot_values[pidx]
                    if abs(ref_val) > 1e-10:
                        h_est[sc, dmrs_sym, ch_idx] = rx_val / ref_val
                    else:
                        h_est[sc, dmrs_sym, ch_idx] = 1.0 + 0j
                
                # Step 2: CDM group frequency averaging (pairs of adjacent pilots)
                # C++ code: cdm_group_size=2, average (i,i+1) pairs if sc_b - sc_a <= 2
                for i in range(0, len(dmrs_sc_list) - 1, 2):
                    j = i + 1
                    sc_a = dmrs_sc_list[i]
                    sc_b = dmrs_sc_list[j]
                    if sc_b - sc_a > 2:
                        continue
                    avg = (h_est[sc_a, dmrs_sym, ch_idx] + h_est[sc_b, dmrs_sym, ch_idx]) / 2.0
                    h_est[sc_a, dmrs_sym, ch_idx] = avg
                    h_est[sc_b, dmrs_sym, ch_idx] = avg
                
                # Step 3: Frequency interpolation (linear)
                h_at_dmrs = [h_est[sc, dmrs_sym, ch_idx] for sc in dmrs_sc_list]
                first_dmrs_sc = dmrs_sc_list[0]
                last_dmrs_sc = dmrs_sc_list[-1]
                
                # Before first pilot: constant extrapolation
                for sc in range(0, first_dmrs_sc):
                    h_est[sc, dmrs_sym, ch_idx] = h_at_dmrs[0]
                
                # Between pilots: linear interpolation
                for i in range(len(dmrs_sc_list) - 1):
                    sc_a = dmrs_sc_list[i]
                    sc_b = dmrs_sc_list[i + 1]
                    h_a = h_at_dmrs[i]
                    h_b = h_at_dmrs[i + 1]
                    dist = sc_b - sc_a
                    if dist <= 0:
                        continue
                    for sc in range(sc_a, sc_b + 1):
                        alpha = (sc - sc_a) / dist
                        h_est[sc, dmrs_sym, ch_idx] = h_a * (1.0 - alpha) + h_b * alpha
                
                # After last pilot: constant extrapolation
                for sc in range(last_dmrs_sc, n_sc):
                    h_est[sc, dmrs_sym, ch_idx] = h_at_dmrs[-1]
    
    # Step 4: Time interpolation
    for rx in range(n_rx_ant):
        for l in range(n_layers):
            ch_idx = rx * n_layers + l
            
            first_dmrs_sym = pilot_syms[0]
            last_dmrs_sym = pilot_syms[-1]
            
            # Before first DMRS symbol: constant
            for sym in range(0, first_dmrs_sym):
                for sc in range(n_sc):
                    h_est[sc, sym, ch_idx] = h_est[sc, first_dmrs_sym, ch_idx]
            
            # Between DMRS symbols: linear
            for i in range(len(pilot_syms) - 1):
                sym_a = pilot_syms[i]
                sym_b = pilot_syms[i + 1]
                dist = sym_b - sym_a
                if dist <= 1:
                    continue
                for sym in range(sym_a + 1, sym_b):
                    alpha = (sym - sym_a) / dist
                    for sc in range(n_sc):
                        h_a = h_est[sc, sym_a, ch_idx]
                        h_b = h_est[sc, sym_b, ch_idx]
                        h_est[sc, sym, ch_idx] = h_a * (1.0 - alpha) + h_b * alpha
            
            # After last DMRS symbol: constant
            for sym in range(last_dmrs_sym + 1, n_sym):
                for sc in range(n_sc):
                    h_est[sc, sym, ch_idx] = h_est[sc, last_dmrs_sym, ch_idx]
    
    return h_est


def cpp_lmmse_equalize(y, h_est, noise_var, n_sc, n_sym, n_rx_ant, n_layers):
    """Replicate C++ MmseEqualizer::equalize() SISO fast path in Python.
    
    For 1x1 SISO:
        G = conj(h) / (|h|^2 + no)
        gy = G * y
        d = |h|^2 / (|h|^2 + no)
        x_hat = gy / d = y / h  (ZF/unbiased LMMSE)
        no_eff = 1/d - 1 = no / |h|^2
    
    Args:
        y: numpy array shape [n_sc*n_sym, n_rx_ant] complex - received data REs (flattened)
        h_est: numpy array shape [n_sc, n_sym, n_rx_ant*n_layers] complex - channel estimates
        noise_var: float - AWGN noise variance
        n_sc, n_sym, n_rx_ant, n_layers: dimensions
    
    Returns:
        x_hat: numpy array shape [n_sc*n_sym, n_layers] complex
        no_eff: numpy array shape [n_sc*n_sym*n_layers] float
    """
    n_re = n_sc * n_sym
    equalized = np.zeros((n_re, n_layers), dtype=complex)
    no_eff = np.zeros(n_re * n_layers, dtype=float)
    
    for i in range(n_re):
        sc = i % n_sc
        sym_idx = i // n_sc
        sym = sym_idx if sym_idx < n_sym else (n_sym - 1)
        
        H = np.zeros((n_rx_ant, n_layers), dtype=complex)
        for rx in range(n_rx_ant):
            for l in range(n_layers):
                ch_idx = rx * n_layers + l
                H[rx, l] = h_est[sc, sym, ch_idx]
        
        y_vec = np.zeros(n_rx_ant, dtype=complex)
        for rx in range(n_rx_ant):
            y_vec[rx] = y[i, rx]
        
        if n_rx_ant == 1 and n_layers == 1:
            h = H[0, 0]
            h_mag_sq = abs(h) ** 2
            gy = np.conj(h) * y_vec[0] / (h_mag_sq + noise_var)
            d = h_mag_sq / (h_mag_sq + noise_var)
            equalized[i, 0] = gy / d
            no_eff[i] = (1.0 / d - 1.0)
        else:
            # General MIMO case using C++ algorithm (matrix inverse)
            Hh = H.conj().T
            HHh = H @ Hh
            for rx in range(n_rx_ant):
                HHh[rx, rx] += noise_var
            G = Hh @ np.linalg.inv(HHh)
            gy = G @ y_vec
            GH = G @ H
            for l in range(n_layers):
                d = GH[l, l].real
                if abs(d) < 1e-12:
                    d = 1e-12
                equalized[i, l] = gy[l] / d
                no_eff[i * n_layers + l] = (1.0 / d - 1.0)
    
    return equalized, no_eff


def main():
    print("=" * 70)
    print("  C++ vs Sionna Cross-Check: LS Channel Estimation & LMMSE Equalization")
    print("=" * 70)
    
    # Configuration matching C++ pdsch_bler_simulation.cpp
    mcs = 7
    n_rb = 25
    num_antenna_ports = 1
    num_layers = 1
    
    cfg = PUSCHConfig()
    cfg.n_size_bwp = n_rb
    cfg.num_antenna_ports = num_antenna_ports
    cfg.num_layers = num_layers
    cfg.transform_precoding = False
    cfg.dmrs.config_type = 1
    cfg.dmrs.length = 1
    cfg.dmrs.additional_position = 0
    cfg.dmrs.num_cdm_groups_without_data = 1
    cfg.dmrs.dmrs_port_set = [0]
    cfg.tb.mcs_index = mcs
    cfg.tb.mcs_table = 1
    
    # Create transmitter in frequency domain
    pusch_tx = PUSCHTransmitter(cfg, output_domain='freq')
    rg = pusch_tx.resource_grid
    
    n_sc = rg.num_effective_subcarriers
    n_sym = rg.num_ofdm_symbols
    print(f"\nResource grid: {n_sc} subcarriers x {n_sym} symbols")
    print(f"Config: MCS={mcs}, n_rb={n_rb}, 1x1 SISO, DMRS Type1, single-symbol, pos=2")
    
    # Extract pilot pattern
    mask = rg.pilot_pattern.mask.numpy()  # [num_tx, num_ports, n_sym, n_sc]
    pilots = rg.pilot_pattern.pilots.numpy()  # [num_tx, num_ports, n_pilots]
    
    # For 1x1: mask[0,0,sym,sc] = 1 where pilots exist
    # Build list of (sym, sc) pilot positions
    pilot_positions = []
    pilot_values = []
    pilot_count = 0
    for sym in range(n_sym):
        for sc in range(n_sc):
            if mask[0, 0, sym, sc] > 0:
                pilot_positions.append((sym, sc))
                pilot_values.append(pilots[0, 0, pilot_count])
                pilot_count += 1
    pilot_values = np.array(pilot_values, dtype=complex)
    print(f"Number of DMRS pilots: {len(pilot_positions)}")
    
    # Find DMRS symbols
    dmrs_syms = sorted(set(sym for sym, sc in pilot_positions))
    print(f"DMRS symbol(s): {dmrs_syms}")
    
    # Test at multiple SNR points
    snr_db_list = [0, 5, 10, 20]
    batch_size = 10
    
    print("\n" + "=" * 70)
    print("  Part 1: Channel Estimation Cross-Check")
    print("=" * 70)
    
    for snr_db in snr_db_list:
        no_val = 1.0 / (10 ** (snr_db / 10.0))
        sigma = np.sqrt(no_val / 2.0)
        
        h_err_all = []
        h_corr_all = []
        
        for b in range(batch_size):
            torch.manual_seed(1000 + b)
            np.random.seed(1000 + b)
            
            # Generate transmit signal
            x, bits = pusch_tx(1)
            # x shape: [batch=1, num_tx=1, num_streams=1, n_sym, n_sc]
            x_np = x[0, 0, 0].detach().numpy()  # [n_sym, n_sc]
            
            # Add AWGN noise (matching C++: complex Gaussian with variance no)
            noise = (np.random.randn(n_sym, n_sc) + 1j * np.random.randn(n_sym, n_sc)) * sigma
            y_np = x_np + noise
            
            # === C++-style LS CE (Python reimplementation) ===
            # y shape: [n_rx_ant=1, n_sym, n_sc]
            y_for_cpp = y_np[np.newaxis, :, :]  # add rx_ant dimension
            h_cpp = cpp_ls_channel_estimate(
                y_for_cpp, pilot_positions, pilot_values,
                n_sc, n_sym, n_rx_ant=1, n_layers=1
            )
            # h_cpp shape: [n_sc, n_sym, 1]
            h_cpp_2d = h_cpp[:, :, 0]  # [n_sc, n_sym]
            
            # === Sionna LS CE ===
            ce = PUSCHLSChannelEstimator(
                rg,
                dmrs_length=1,
                dmrs_additional_position=0,
                num_cdm_groups_without_data=1,
                interpolation_type="lin"
            )
            # y input: [batch, num_rx, num_rx_ant, n_sym, n_sc]
            y_t = torch.tensor(y_np[np.newaxis, np.newaxis, np.newaxis, :, :], dtype=torch.complex64)
            no_t = torch.tensor(no_val, dtype=torch.float32)
            h_sionna, err_var_sionna = ce(y_t, no_t)
            # h_sionna: [batch, num_rx=1, num_rx_ant=1, num_tx=1, num_streams=1, n_sym, n_sc]
            h_sionna_np = h_sionna[0, 0, 0, 0, 0].detach().numpy()  # [n_sym, n_sc]
            h_sionna_2d = h_sionna_np.T  # transpose to [n_sc, n_sym]
            
            # Compare channel estimates
            # Compare at all REs
            diff = h_cpp_2d - h_sionna_2d
            err = np.abs(diff)
            h_mag_ref = np.abs(h_cpp_2d)
            rel_err = err / (h_mag_ref + 1e-12)
            
            h_err_all.append(err)
            h_corr_all.append(rel_err)
        
        h_err_all = np.concatenate([e.flatten() for e in h_err_all])
        h_corr_all = np.concatenate([e.flatten() for e in h_corr_all])
        
        print(f"\n  SNR = {snr_db} dB (no = {no_val:.6f}):")
        print(f"    |h_cpp - h_sionna| Mean:   {np.mean(h_err_all):.6e}")
        print(f"    |h_cpp - h_sionna| Max:    {np.max(h_err_all):.6e}")
        print(f"    |h_cpp - h_sionna| Median: {np.median(h_err_all):.6e}")
        print(f"    Relative error Mean:       {np.mean(h_corr_all)*100:.4f}%")
        print(f"    Relative error 99th pct:   {np.percentile(h_corr_all, 99)*100:.4f}%")
    
    print("\n" + "=" * 70)
    print("  Part 2: LMMSE Equalization Cross-Check")
    print("=" * 70)
    print("  (Using identical h_est from Sionna CE, then comparing EQ outputs)")
    
    for snr_db in snr_db_list:
        no_val = 1.0 / (10 ** (snr_db / 10.0))
        sigma = np.sqrt(no_val / 2.0)
        
        x_err_abs = []
        x_err_rel = []
        no_eff_err_abs = []
        
        for b in range(batch_size):
            torch.manual_seed(2000 + b)
            np.random.seed(2000 + b)
            
            x, bits = pusch_tx(1)
            x_np = x[0, 0, 0].detach().numpy()  # [n_sym, n_sc]
            noise = (np.random.randn(n_sym, n_sc) + 1j * np.random.randn(n_sym, n_sc)) * sigma
            y_np = x_np + noise
            
            # Get Sionna channel estimate first (we'll use this for both equalizers
            # to isolate equalizer differences)
            ce = PUSCHLSChannelEstimator(
                rg, dmrs_length=1, dmrs_additional_position=0,
                num_cdm_groups_without_data=1, interpolation_type="lin"
            )
            y_t = torch.tensor(y_np[np.newaxis, np.newaxis, np.newaxis, :, :], dtype=torch.complex64)
            no_t = torch.tensor(no_val, dtype=torch.float32)
            h_sionna, _ = ce(y_t, no_t)
            # h_sionna: [batch=1, num_rx=1, num_rx_ant=1, num_tx=1, num_streams=1, n_sym, n_sc]
            h_for_eq = h_sionna[0, 0, 0, 0, 0].T.detach().numpy()  # [n_sc, n_sym]
            
            # Flatten to RE order (C++ uses sc-first, then sym)
            # C++ Equalizer expects rx_symbols as n_re x n_rx_ant where n_re = n_sc*n_sym
            # and channel_est as n_sc x n_sym x n_ch
            y_flat = np.zeros((n_sc * n_sym, 1), dtype=complex)
            h_cube = np.zeros((n_sc, n_sym, 1), dtype=complex)
            x_flat = np.zeros(n_sc * n_sym, dtype=complex)
            
            for sym in range(n_sym):
                for sc in range(n_sc):
                    idx = sym * n_sc + sc  # Wait, C++ uses sc-first: i = sym*n_sc + sc?
                    # Let me check C++ indexing in Equalizer:
                    # sc = i % n_ch_sc (n_sc), sym = i / n_ch_sc
                    # So for i=0: sc=0, sym=0; i=1: sc=1, sym=0; ...
                    # This means REs are ordered sc first within each symbol!
                    pass
            
            # Actually looking at C++ Equalizer:
            # for (int i = 0; i < n_re; i++) {
            #     sc = i % n_ch_sc;  // n_sc
            #     sym = i / n_ch_sc;
            # }
            # And PdschProcessor passes rx_pdsch from extract_pdsch_from_grid
            # Let me check the RE ordering
            # In C++ PdschProcessor, rx_pdsch is built with n_re rows (PDSCH REs only)
            # But for the whole grid cross-check, let's use sc-major ordering:
            # i = sym * n_sc + sc -> for sym=0, sc=0..n_sc-1, then sym=1, etc.
            # Wait no: i = sym*n_sc + sc gives sc = i % n_sc, sym = i/n_sc.
            # Let's verify: i=0 -> sc=0, sym=0. i=n_sc -> sc=0, sym=1. Yes.
            
            for sym in range(n_sym):
                for sc in range(n_sc):
                    i = sym * n_sc + sc
                    y_flat[i, 0] = y_np[sym, sc]
                    h_cube[sc, sym, 0] = h_for_eq[sc, sym]
                    x_flat[i] = x_np[sym, sc]
            
            # === C++-style LMMSE equalization ===
            x_cpp, no_eff_cpp = cpp_lmmse_equalize(
                y_flat, h_cube, no_val, n_sc, n_sym, n_rx_ant=1, n_layers=1
            )
            x_cpp_flat = x_cpp[:, 0]  # [n_re]
            
            # === Sionna LMMSE equalization (per-RE) ===
            # Call lmmse_equalizer for each RE with h as scalar, s = no
            x_sionna_flat = np.zeros(n_sc * n_sym, dtype=complex)
            no_eff_sionna_flat = np.zeros(n_sc * n_sym, dtype=float)
            
            for i in range(n_sc * n_sym):
                sc = i % n_sc
                sym = i // n_sc
                h_scalar = h_for_eq[sc, sym]
                y_scalar = y_np[sym, sc]
                
                # lmmse_equalizer expects: y [...,M], h [...,M,K], s [...,M,M]
                y_t = torch.tensor([y_scalar], dtype=torch.complex64)
                h_t = torch.tensor([[h_scalar]], dtype=torch.complex64)
                s_t = torch.tensor([[no_val]], dtype=torch.complex64)
                
                x_hat_t, no_eff_t = lmmse_equalizer(y_t, h_t, s_t, whiten_interference=False, precision='single')
                x_sionna_flat[i] = x_hat_t[0].detach().numpy()
                no_eff_sionna_flat[i] = no_eff_t[0].detach().numpy()
            
            # Compare
            x_diff = x_cpp_flat - x_sionna_flat
            x_err = np.abs(x_diff)
            x_ref = np.abs(x_cpp_flat)
            x_rel = x_err / (x_ref + 1e-12)
            
            ne_diff = no_eff_cpp[:n_sc*n_sym] - no_eff_sionna_flat
            
            x_err_abs.append(x_err)
            x_err_rel.append(x_rel)
            no_eff_err_abs.append(np.abs(ne_diff))
        
        x_err_abs = np.concatenate(x_err_abs)
        x_err_rel = np.concatenate(x_err_rel)
        no_eff_err_abs = np.concatenate(no_eff_err_abs)
        
        print(f"\n  SNR = {snr_db} dB (no = {no_val:.6f}):")
        print(f"    Equalized symbols:")
        print(f"      |x_cpp - x_sionna| Mean:   {np.mean(x_err_abs):.6e}")
        print(f"      |x_cpp - x_sionna| Max:    {np.max(x_err_abs):.6e}")
        print(f"      Relative error Mean:       {np.mean(x_err_rel)*100:.6f}%")
        print(f"    Effective noise variance:")
        print(f"      |no_eff_cpp - no_eff_sionna| Mean: {np.mean(no_eff_err_abs):.6e}")
        print(f"      |no_eff_cpp - no_eff_sionna| Max:  {np.max(no_eff_err_abs):.6e}")
    
    print("\n" + "=" * 70)
    print("  Part 3: End-to-end CE+EQ Cross-Check (using both CE implementations)")
    print("=" * 70)
    
    for snr_db in snr_db_list:
        no_val = 1.0 / (10 ** (snr_db / 10.0))
        sigma = np.sqrt(no_val / 2.0)
        
        x_err_e2e = []
        
        for b in range(batch_size):
            torch.manual_seed(3000 + b)
            np.random.seed(3000 + b)
            
            x, bits = pusch_tx(1)
            x_np = x[0, 0, 0].detach().numpy()
            noise = (np.random.randn(n_sym, n_sc) + 1j * np.random.randn(n_sym, n_sc)) * sigma
            y_np = x_np + noise
            
            # C++ CE + C++ EQ
            y_for_cpp = y_np[np.newaxis, :, :]
            h_cpp = cpp_ls_channel_estimate(
                y_for_cpp, pilot_positions, pilot_values,
                n_sc, n_sym, 1, 1
            )
            y_flat = np.zeros((n_sc * n_sym, 1), dtype=complex)
            h_cube = h_cpp  # already [n_sc, n_sym, 1]
            for sym in range(n_sym):
                for sc in range(n_sc):
                    i = sym * n_sc + sc
                    y_flat[i, 0] = y_np[sym, sc]
            x_cpp, _ = cpp_lmmse_equalize(y_flat, h_cube, no_val, n_sc, n_sym, 1, 1)
            x_cpp_e2e = x_cpp[:, 0]
            
            # Sionna CE + Sionna EQ
            ce = PUSCHLSChannelEstimator(
                rg, dmrs_length=1, dmrs_additional_position=0,
                num_cdm_groups_without_data=1, interpolation_type="lin"
            )
            y_t = torch.tensor(y_np[np.newaxis, np.newaxis, np.newaxis, :, :], dtype=torch.complex64)
            no_t = torch.tensor(no_val, dtype=torch.float32)
            h_sionna, _ = ce(y_t, no_t)
            h_eq = h_sionna[0, 0, 0, 0, 0].T.detach().numpy()  # [n_sc, n_sym]
            
            x_sionna_e2e = np.zeros(n_sc * n_sym, dtype=complex)
            for i in range(n_sc * n_sym):
                sc = i % n_sc
                sym = i // n_sc
                h_s = h_eq[sc, sym]
                y_s = y_np[sym, sc]
                y_tt = torch.tensor([y_s], dtype=torch.complex64)
                h_tt = torch.tensor([[h_s]], dtype=torch.complex64)
                s_tt = torch.tensor([[no_val]], dtype=torch.complex64)
                x_hat_t, _ = lmmse_equalizer(y_tt, h_tt, s_tt, whiten_interference=False, precision='single')
                x_sionna_e2e[i] = x_hat_t[0].detach().numpy()
            
            diff_e2e = np.abs(x_cpp_e2e - x_sionna_e2e)
            x_err_e2e.append(diff_e2e)
        
        x_err_e2e = np.concatenate(x_err_e2e)
        print(f"\n  SNR = {snr_db} dB:")
        print(f"    E2E |x_cpp - x_sionna| Mean:   {np.mean(x_err_e2e):.6e}")
        print(f"    E2E |x_cpp - x_sionna| Max:    {np.max(x_err_e2e):.6e}")
        print(f"    E2E |x_cpp - x_sionna| Median: {np.median(x_err_e2e):.6e}")
    
    print("\n" + "=" * 70)
    print("  Summary")
    print("=" * 70)
    print("""
  Key implementation differences between C++ and Sionna:
  
  1. CHANNEL ESTIMATION (LS):
     - Both use identical algorithm: LS at pilots -> CDM pair averaging -> 
       linear frequency interpolation -> constant/linear time interpolation
     - C++ uses double precision (float64); Sionna uses single (complex64)
     - Pilot PN sequences may differ (C++ vs Sionna gold sequence), but LS
       division cancels this when reference matches transmitted pilots.
     - Expected error: ~1e-7 from float32 rounding.
  
  2. LMMSE EQUALIZATION:
     - For 1x1 SISO: identical algorithm: G=conj(h)/(|h|^2+no), x_hat=Gy/(G*h),
       no_eff=1/(G*h)-1, which simplifies to x_hat=y/h (ZF) with no_eff=no/|h|^2
     - Sionna by default whitens interference before LMMSE (whiten_interference=True),
       which is mathematically equivalent for SISO white noise.
     - Key difference: Sionna's receiver adds channel estimation error variance 
       (err_var) to the noise covariance S, while C++ passes only AWGN noise_var
       to the equalizer. For this cross-check, we use S = no*I (matching C++).
     - Expected error: ~1e-7 from float32 rounding.
  
  3. PRACTICAL PERFORMANCE IMPACT:
     - The Sionna err_var addition slightly increases effective noise for EQ,
       which marginally degrades LLR quality vs C++ which ignores CE error.
     - Precision difference (double vs float) gives C++ ~1e-15 vs Sionna ~1e-7
       numerical accuracy.
     - LLR computation: C++ uses max-log; Sionna may use different demapping.
     - LDPC decoding: C++ uses 25 iterations; Sionna defaults to 20.
""")


if __name__ == "__main__":
    main()
