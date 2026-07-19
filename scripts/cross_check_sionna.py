#!/usr/bin/env python3
"""Cross-check nr-link-simulator (C++) vs Sionna at each PDSCH chain module.

Loads intermediate results from C++ cross_check_tx program and compares
with Sionna's processing of the same input transport block bits.
"""

import os
import sys
import numpy as np
import torch

sys.path.insert(0, '/workspace/sionna/src')

from sionna.phy.fec.crc import CRCEncoder, CRCDecoder
from sionna.phy.fec.ldpc import LDPC5GEncoder, LDPC5GDecoder
from sionna.phy.fec.scrambling import TB5GScrambler
from sionna.phy.mapping import Mapper, Demapper, Constellation
from sionna.phy.nr import TBEncoder, TBDecoder
from sionna.phy.nr.utils import calculate_tb_size, generate_prng_seq

CROSS_DIR = '/workspace/nr-link-simulator/results/cross_check/'

def load_npy(name):
    return np.load(os.path.join(CROSS_DIR, name))

def compare_bits(name, cpp_bits, sionna_bits, max_diff=20):
    """Compare two bit arrays and report mismatches."""
    cpp = np.asarray(cpp_bits, dtype=np.int32).flatten()
    sio = np.asarray(sionna_bits, dtype=np.int32).flatten()
    n = min(len(cpp), len(sio))
    if len(cpp) != len(sio):
        print(f"  [LENGTH MISMATCH] {name}: C++={len(cpp)}, Sionna={len(sio)}")
        return False
    diffs = np.where(cpp[:n] != sio[:n])[0]
    n_diff = len(diffs)
    if n_diff == 0:
        print(f"  [OK] {name}: {len(cpp)} bits, exact match")
        return True
    else:
        print(f"  [MISMATCH] {name}: {n_diff}/{len(cpp)} bits differ")
        for idx in diffs[:max_diff]:
            print(f"    bit[{idx}]: C++={cpp[idx]}, Sionna={sio[idx]}")
        return False

def compare_complex(name, cpp_sym, sionna_sym, tol=1e-6):
    cpp = np.asarray(cpp_sym, dtype=np.complex128).flatten()
    sio = np.asarray(sionna_sym, dtype=np.complex128).flatten()
    n = min(len(cpp), len(sio))
    if len(cpp) != len(sio):
        print(f"  [LENGTH MISMATCH] {name}: C++={len(cpp)}, Sionna={len(sio)}")
        return False
    diff = np.abs(cpp[:n] - sio[:n])
    max_diff = np.max(diff)
    mean_diff = np.mean(diff)
    if max_diff < tol:
        print(f"  [OK] {name}: {len(cpp)} symbols, max|diff|={max_diff:.2e}")
        return True
    else:
        print(f"  [MISMATCH] {name}: {len(cpp)} symbols, max|diff|={max_diff:.6f}, mean|diff|={mean_diff:.6f}")
        bad_idx = np.where(diff > tol)[0]
        for idx in bad_idx[:10]:
            print(f"    sym[{idx}]: C++={cpp[idx]:.8f}, Sionna={sio[idx]:.8f}, diff={diff[idx]:.2e}")
        return False

def compare_float(name, cpp_val, sionna_val, tol=1e-6):
    cpp = np.asarray(cpp_val, dtype=np.float64).flatten()
    sio = np.asarray(sionna_val, dtype=np.float64).flatten()
    n = min(len(cpp), len(sio))
    if len(cpp) != len(sio):
        print(f"  [LENGTH MISMATCH] {name}: C++={len(cpp)}, Sionna={len(sio)}")
        return False
    diff = np.abs(cpp[:n] - sio[:n])
    max_diff = np.max(diff)
    if max_diff < tol:
        print(f"  [OK] {name}: {len(cpp)} values, max|diff|={max_diff:.2e}")
        return True
    else:
        print(f"  [MISMATCH] {name}: {len(cpp)} values, max|diff|={max_diff:.6f}")
        bad_idx = np.where(diff > tol)[0]
        for idx in bad_idx[:5]:
            print(f"    [{idx}]: C++={cpp[idx]:.8f}, Sionna={sio[idx]:.8f}, diff={diff[idx]:.2e}")
        return False

def main():
    mcs = 11
    n_prb = 3
    n_re_per_prb = 156
    n_layers = 1
    qm = 6
    r1024 = 466
    R = r1024 / 1024.0
    tbs = 1288
    G = n_prb * n_re_per_prb * qm * n_layers
    num_bits_per_symbol = qm
    n_sym = G // qm
    es_n0_db = 10.0
    rv = 0
    n_iter = 20

    print("=" * 70)
    print("PDSCH TX Chain Cross-Check: C++ (nr-link-simulator) vs Sionna")
    print(f"MCS={mcs}, Qm={qm}, R={R:.4f}, TBS={tbs}, G={G}, n_sym={n_sym}")
    print("=" * 70)

    print("\n--- Loading C++ intermediate results ---")
    with open(os.path.join(CROSS_DIR, 'cb_info.txt'), 'r') as f:
        print(f.read())

    cpp_tb_bits = load_npy('tb_bits.npy')
    cpp_tb_crc = load_npy('tb_crc_bits.npy')
    cpp_cb0_input = load_npy('cb0_ldpc_input.npy')
    cpp_cb0_output = load_npy('cb0_ldpc_output.npy')
    cpp_cb0_punctured = load_npy('cb0_punctured.npy')
    cpp_cb0_nofiller = load_npy('cb0_no_filler.npy')
    cpp_cb0_rm = load_npy('cb0_rate_matched.npy')
    cpp_rm_all = load_npy('rm_bits_all.npy')
    cpp_scrambled = load_npy('scrambled_bits.npy')
    cpp_modulated = load_npy('modulated_symbols.npy')

    print(f"TB bits shape: {cpp_tb_bits.shape}, CRC bits shape: {cpp_tb_crc.shape}")
    print(f"LDPC input shape: {cpp_cb0_input.shape}, LDPC output shape: {cpp_cb0_output.shape}")

    tb_bits_t = torch.tensor(cpp_tb_bits, dtype=torch.float32).unsqueeze(0)

    print("\n--- Step 1: TB CRC Attachment ---")
    crc_enc = CRCEncoder("CRC24A" if tbs > 3824 else "CRC16")
    crc_bits = crc_enc(tb_bits_t)
    crc_bits_np = crc_bits.numpy().flatten().astype(np.uint8)
    print(f"  CRC polynomial: {'CRC24A' if tbs > 3824 else 'CRC16'}")
    compare_bits("TB+CRC bits", cpp_tb_crc, crc_bits_np)

    print("\n--- Step 2: LDPC Encoding (Sionna LDPC5GEncoder) ---")
    k = tbs + (24 if tbs > 3824 else 16)
    enc = LDPC5GEncoder(k=k, n=G, num_bits_per_symbol=num_bits_per_symbol)
    print(f"  Sionna encoder: k={enc.k}, n={enc.n}, bg={enc._bg}, z={enc.z}, k_b={enc._k_b}")
    print(f"  k_filler={enc.k_filler}, n_ldpc={enc.n_ldpc}, n_cb={enc.n_cb}, n_cb_comp={enc.n_cb_comp}")
    print(f"  rv_starts={enc.rv_starts}")

    cw = enc(crc_bits, rv=[rv])
    cw_np = cw.numpy().flatten().astype(np.uint8)
    compare_bits("Rate-matched codeword (after BICM)", cpp_rm_all, cw_np)

    print("\n--- Step 2b: Compare LDPC encoder full output (before rate matching) ---")
    u_fill = torch.cat([crc_bits, torch.zeros(1, enc.k_filler, dtype=torch.float32)], dim=1)
    c_full = enc._encode_fast(u_fill.to(torch.float32))
    c_full = c_full.reshape(1, enc.n_ldpc)
    c_full_np = c_full.numpy().flatten().astype(np.uint8)
    compare_bits("Full LDPC codeword (N bits)", cpp_cb0_output, c_full_np)

    print("\n--- Step 2c: After puncturing first 2Z ---")
    c_punct = c_full[:, 2*enc.z:]
    c_punct_np = c_punct.numpy().flatten().astype(np.uint8)
    compare_bits("After puncturing 2Z", cpp_cb0_punctured, c_punct_np)

    print("\n--- Step 2d: After removing filler bits ---")
    c_nofill = torch.cat([c_full[:, :enc.k], c_full[:, enc.k_ldpc:]], dim=1)
    c_nofill_np = c_nofill.numpy().flatten().astype(np.uint8)
    compare_bits("After removing fillers (before puncturing)", cpp_cb0_output[:enc.k].astype(np.uint8), c_full[:, :enc.k].numpy().flatten().astype(np.uint8))

    c_nofill_punct = c_nofill[:, 2*enc.z:]
    c_nofill_punct_np = c_nofill_punct.numpy().flatten().astype(np.uint8)
    compare_bits("Compressed RM buffer (after 2Z puncture + filler removal)",
                 cpp_cb0_nofiller, c_nofill_punct_np)

    print("\n--- Step 3: Scrambling ---")
    n_rnti = 1
    n_id = 0
    scrambler = TB5GScrambler(n_rnti=n_rnti, n_id=n_id, channel_type='PUSCH', binary=True)
    scrambled_s = scrambler(torch.tensor(cpp_rm_all, dtype=torch.float32).unsqueeze(0))
    scrambled_s_np = scrambled_s.numpy().flatten().astype(np.uint8)

    cinit_cpp = 1
    cinit_sionna = n_rnti * (2**15) + n_id
    print(f"  C++ c_init = {cinit_cpp}")
    print(f"  Sionna c_init (PUSCH, n_rnti={n_rnti}, n_id={n_id}) = {cinit_sionna}")
    seq_cpp = generate_prng_seq(G, cinit_cpp)
    seq_sio = generate_prng_seq(G, cinit_sionna)
    seq_match = np.all(seq_cpp == seq_sio)
    print(f"  Scrambling sequences match: {seq_match}")
    compare_bits("Scrambled bits (Sionna scrambler)", cpp_scrambled, scrambled_s_np)

    scrambled_same_cinit = (cpp_rm_all.astype(np.int32) ^ seq_cpp.astype(np.int32)).astype(np.uint8)
    compare_bits("Scrambled bits (same c_init=1)", cpp_scrambled, scrambled_same_cinit)

    print("\n--- Step 4: Constellation Mapping (Modulation) ---")
    constellation = Constellation("qam", num_bits_per_symbol)
    mapper = Mapper(constellation=constellation)

    bits_for_map = torch.tensor(cpp_scrambled, dtype=torch.float32).reshape(1, -1)
    symbols_s = mapper(bits_for_map)
    symbols_s_np = symbols_s.numpy().flatten()

    cpp_pwr = np.mean(np.abs(cpp_modulated)**2)
    sio_pwr = np.mean(np.abs(symbols_s_np)**2)
    print(f"  C++ constellation power: {cpp_pwr:.6f}")
    print(f"  Sionna constellation power: {sio_pwr:.6f}")

    cpp_mod = cpp_modulated / np.sqrt(cpp_pwr) if abs(cpp_pwr - 1.0) > 0.01 else cpp_modulated
    compare_complex("Modulated symbols", cpp_mod, symbols_s_np, tol=1e-4)

    print("\n--- Step 4b: Check constellation point mapping ---")
    cpp_unique = np.unique(np.round(np.real(cpp_modulated), 6))
    sio_unique = np.unique(np.round(np.real(symbols_s_np), 6))
    print(f"  C++ real axis unique values: {cpp_unique}")
    print(f"  Sionna real axis unique values: {sio_unique}")

    print("\n--- Step 5: LLR Demodulation (AWGN, Es/N0 = 10 dB) ---")
    no_val = 1.0 / (10 ** (es_n0_db / 10.0))
    noise_std = np.sqrt(no_val / 2.0)

    np.random.seed(12345)
    noise = (np.random.randn(n_sym) + 1j*np.random.randn(n_sym)) * noise_std
    rx_sym = cpp_modulated + noise

    demapper = Demapper("app", constellation=constellation)
    no_t = torch.tensor(no_val, dtype=torch.float32)
    rx_t = torch.tensor(rx_sym, dtype=torch.complex64).unsqueeze(0)
    llr_s = demapper(rx_t, no_t)
    llr_s_np = llr_s.numpy().flatten()

    cpp_rx = load_npy('rx_symbols.npy')
    cpp_llr = load_npy('llr_values.npy')
    print(f"  Note: C++ uses its own noise (different seed), so LLRs won't match directly")
    print(f"  C++ LLR range: [{cpp_llr.min():.2f}, {cpp_llr.max():.2f}], mean|LLR|={np.mean(np.abs(cpp_llr)):.2f}")
    print(f"  Sionna LLR range: [{llr_s_np.min():.2f}, {llr_s_np.max():.2f}], mean|LLR|={np.mean(np.abs(llr_s_np)):.2f}")

    print("\n" + "=" * 70)
    print("RX Chain Analysis: LLR computation and LDPC decoding")
    print("=" * 70)

    print("\n--- Step 6: Verify Sionna LLR formula matches C++ ---")
    test_sym = torch.tensor([0.5+0.3j, -0.2+0.8j, 1.0-0.5j], dtype=torch.complex64).unsqueeze(0)
    test_no = torch.tensor(0.1, dtype=torch.float32)
    test_llr_s = demapper(test_sym, test_no).numpy().flatten()

    def cpp_qam64_llr(y, noise_var):
        m = 3
        M = 1 << m
        ms_pam = sum((2*k-M+1)**2 for k in range(M)) / M
        a_pam = 1.0 / np.sqrt(2.0 * ms_pam)
        var_per_dim = max(noise_var / 2.0, 1e-12)
        out = np.zeros(6)
        for dim, ydim in enumerate([y.real, y.imag]):
            for b in range(m):
                d0 = 1e30; d1 = 1e30
                for k in range(M):
                    s = (2.0*k - M + 1.0) * a_pam
                    dist = (ydim - s)**2
                    g = k ^ (k >> 1)
                    bit_val = (g >> (m - 1 - b)) & 1
                    if bit_val == 0:
                        d0 = min(d0, dist)
                    else:
                        d1 = min(d1, dist)
                out[dim*m + b] = (d1 - d0) / (2.0 * var_per_dim)
        return out

    test_sym_np = test_sym.numpy().flatten()
    for i, sym in enumerate(test_sym_np):
        cpp_llr_test = cpp_qam64_llr(sym, 0.1)
        sio_llr_test = test_llr_s[i*6:(i+1)*6]
        llr_diff = np.max(np.abs(cpp_llr_test - sio_llr_test))
        print(f"  Symbol {i} ({sym:.4f}): max LLR diff = {llr_diff:.6e}")

    print("\n--- Step 7: End-to-End with TBEncoder/TBDecoder (Sionna) ---")
    tb_enc = TBEncoder(
        target_tb_size=tbs, num_coded_bits=G, target_coderate=R,
        num_bits_per_symbol=qm, num_layers=n_layers, channel_type='PUSCH',
        use_scrambler=True, verbose=False
    )
    tb_dec = TBDecoder(tb_enc, num_iter=n_iter, hard_out=True, verbose=False,
                       constellation=constellation, num_bits_per_symbol=qm)

    n_test = 10
    bit_errors = 0
    block_errors = 0
    for snr_test in [8, 10, 12, 14]:
        no_test = 1.0 / (10**(snr_test/10.0))
        n_err = 0; n_blk = 0
        for _ in range(n_test):
            b = torch.randint(0, 2, (1, tbs), dtype=torch.float32)
            c = tb_enc(b)
            x = mapper(c)
            ns = np.sqrt(no_test/2.0)
            noise_t = torch.complex(
                torch.randn_like(x, dtype=torch.float32)*ns,
                torch.randn_like(x, dtype=torch.float32)*ns
            ).to(x.dtype)
            y = x + noise_t
            llr = demapper(y, torch.tensor(no_test, dtype=torch.float32))
            _, crc_ok = tb_dec(llr)
            if not crc_ok[0]:
                n_err += 1
            n_blk += 1
        bler = n_err / n_blk
        print(f"  SNR={snr_test}dB: BLER={bler:.2f} ({n_err}/{n_blk} errors)")

    print("\n" + "=" * 70)
    print("SUMMARY OF KEY DIFFERENCES")
    print("=" * 70)
    print("""
Key areas where C++ and Sionna may differ:
1. Scrambling c_init: C++ uses c_init = slot_idx+1 = 1; Sionna TB5GScrambler uses
   c_init = n_rnti*2^15 + n_id. This doesn't affect BLER (scrambling is self-inverse)
   BUT if descrambling uses same seed as scrambling, it's transparent.

2. LDPC Decoder algorithm:
   - C++: Offset Min-Sum (offset=0.5), layered? No - it's flooding schedule
   - Sionna: Belief Propagation (sum-product) with flooding schedule
   Min-Sum is an approximation of BP that causes ~0.1-0.5dB loss depending on offset.

3. LLR computation: Both use max-log-MAP approximation, check normalization.

4. Rate matching / BICM interleaving: Check bit ordering matches.

5. Noise variance scaling: Check that Es/N0 definition matches
   (C++ defines noise_var = 1/snr_lin, which is N0 for unit-power constellation).
""")

if __name__ == '__main__':
    main()
