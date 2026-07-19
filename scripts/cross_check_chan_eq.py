import numpy as np
import os

data_dir = "/workspace/nr-link-simulator/results/cross_check_chan/"

tx_grid = np.load(data_dir + "tx_grid.npy")
rx_grid = np.load(data_dir + "rx_grid.npy")
h_est_cpp = np.load(data_dir + "h_est.npy")
equalized_cpp = np.load(data_dir + "equalized.npy")
eff_noise_cpp = np.load(data_dir + "eff_noise_var.npy")
rx_pdsch = np.load(data_dir + "rx_pdsch.npy")
h_at_pdsch = np.load(data_dir + "h_at_pdsch.npy")
llr_cpp = np.load(data_dir + "llr_values.npy")

with open(data_dir + "config.txt") as f:
    config = {}
    for line in f:
        k, v = line.strip().split("=")
        config[k] = v

print("=== Config ===")
for k, v in config.items():
    print(f"  {k} = {v}")

n_ant, n_sym, n_sc = tx_grid.shape
n_rb = n_sc // 12
es_n0_db = float(config["es_n0_db"])
noise_var = float(config["noise_var"])
n_pdsch_re = int(config["n_pdsch_re"])
qm = int(config["qm"])

print(f"\n=== Grid shapes ===")
print(f"  tx_grid: {tx_grid.shape}")
print(f"  rx_grid: {rx_grid.shape}")
print(f"  h_est_cpp: {h_est_cpp.shape}")
print(f"  equalized_cpp: {equalized_cpp.shape}")
print(f"  rx_pdsch: {rx_pdsch.shape}")
print(f"  h_at_pdsch: {h_at_pdsch.shape}")

TYPE1_RE_POS = [0, 2, 4, 6, 8, 10]
dmrs_sym = 2

dmrs_sc = []
for prb in range(n_rb):
    for re_pos in TYPE1_RE_POS:
        dmrs_sc.append(prb * 12 + re_pos)
dmrs_sc = np.array(sorted(dmrs_sc))
print(f"\n=== DMRS Symbol (symbol {dmrs_sym}) ===")
print(f"  DMRS subcarriers ({len(dmrs_sc)}): {dmrs_sc}")

rx_dmrs = rx_grid[0, dmrs_sym, dmrs_sc]
tx_dmrs = tx_grid[0, dmrs_sym, dmrs_sc]

h_ls = rx_dmrs / tx_dmrs
print(f"\n  LS estimates at DMRS positions (first 6): {h_ls[:6]}")
print(f"  LS mean: {np.mean(h_ls):.6f} (should be ~1+0j for AWGN)")

num_cdm_groups_without_data = 1
n = 2 * num_cdm_groups_without_data
k_pairs = len(dmrs_sc) // n
print(f"\n  CDM freq averaging: n={n}, k_pairs={k_pairs}")

h_ls_reshaped = h_ls.reshape(k_pairs, n)
h_avg = h_ls_reshaped.sum(axis=-1, keepdims=True) / 2.0
h_ls_avg = np.repeat(h_avg, n, axis=-1).reshape(-1)
print(f"  After CDM averaging (first 6): {h_ls_avg[:6]}")
print(f"  Averaged mean: {np.mean(h_ls_avg):.6f}")

print(f"\n=== C++ h_est at DMRS positions ===")
h_cpp_at_dmrs = h_est_cpp[0, dmrs_sym, dmrs_sc]
print(f"  C++ h_est at DMRS (first 6): {h_cpp_at_dmrs[:6]}")
print(f"  C++ mean at DMRS: {np.mean(h_cpp_at_dmrs):.6f}")

print(f"\n=== Compare LS+CDM with C++ h_est at DMRS positions ===")
diff_dmrs = h_ls_avg - h_cpp_at_dmrs
print(f"  Max abs diff: {np.max(np.abs(diff_dmrs)):.2e}")
print(f"  Mean abs diff: {np.mean(np.abs(diff_dmrs)):.2e}")

print(f"\n=== C++ h_est at ALL subcarriers (DMRS sym) ===")
h_cpp_all = h_est_cpp[0, dmrs_sym, :]
print(f"  Mean: {np.mean(h_cpp_all):.6f}")
print(f"  Std: {np.std(h_cpp_all):.6f}")
print(f"  Min/Max |h|: {np.min(np.abs(h_cpp_all)):.4f} / {np.max(np.abs(h_cpp_all)):.4f}")

print(f"\n=== Check data symbols (non-DMRS) h_est ===")
data_sc = [sc for sc in range(n_sc) if sc not in dmrs_sc]
h_cpp_data = h_est_cpp[0, dmrs_sym, data_sc]
print(f"  At data SCs: mean={np.mean(h_cpp_data):.6f}, std={np.std(h_cpp_data):.6f}")

print(f"\n=== MMSE Equalization Check (SISO AWGN) ===")
print(f"  noise_var = {noise_var:.6f}")

h_pdsch = h_at_pdsch[0, 0, :]
y_pdsch = rx_pdsch[:, 0]

print(f"\n  h_at_pdsch: shape={h_pdsch.shape}, mean={np.mean(h_pdsch):.6f}")
print(f"  rx_pdsch: shape={y_pdsch.shape}")

gy = np.conj(h_pdsch) * y_pdsch / (np.abs(h_pdsch)**2 + noise_var)
d = np.abs(h_pdsch)**2 / (np.abs(h_pdsch)**2 + noise_var)
eq_py = gy / d
eff_noise_py = 1.0/d - 1.0

print(f"\n  Python MMSE eq: mean(|eq|)={np.mean(np.abs(eq_py)):.4f}")
print(f"  C++ equalized: mean(|eq|)={np.mean(np.abs(equalized_cpp[:,0])):.4f}")

eq_diff = eq_py - equalized_cpp[:, 0]
print(f"\n  Equalized max abs diff: {np.max(np.abs(eq_diff)):.2e}")
print(f"  Equalized mean abs diff: {np.mean(np.abs(eq_diff)):.2e}")

print(f"\n  Eff noise: py mean={np.mean(eff_noise_py):.6f}, cpp mean={np.mean(eff_noise_cpp):.6f}")
noise_diff = eff_noise_py - eff_noise_cpp
print(f"  Eff noise max diff: {np.max(np.abs(noise_diff)):.2e}")

print(f"\n=== Compare with original tx symbols (AWGN, h=1) ===")
from collections import defaultdict

pdsch_sc_list = []
pdsch_sym_list = []
for sym in range(n_sym):
    if sym == dmrs_sym:
        continue
    for sc in range(n_sc):
        pdsch_sc_list.append(sc)
        pdsch_sym_list.append(sym)

tx_pdsch = []
for sym, sc in zip(pdsch_sym_list, pdsch_sc_list):
    tx_pdsch.append(tx_grid[0, sym, sc])
tx_pdsch = np.array(tx_pdsch)

print(f"  tx_pdsch count: {len(tx_pdsch)}, n_pdsch_re: {n_pdsch_re}")

eq_cmp = equalized_cpp[:, 0]
nmse = np.mean(np.abs(eq_cmp - tx_pdsch)**2) / np.mean(np.abs(tx_pdsch)**2)
print(f"  NMSE of C++ equalized vs tx: {nmse:.6f} ({10*np.log10(nmse):.2f} dB)")

nmse_py = np.mean(np.abs(eq_py - tx_pdsch)**2) / np.mean(np.abs(tx_pdsch)**2)
print(f"  NMSE of Python equalized vs tx: {nmse_py:.6f} ({10*np.log10(nmse_py):.2f} dB)")

print(f"\n=== LLR comparison (check first 20) ===")
for i in range(min(20, len(llr_cpp))):
    print(f"  [{i}] llr={llr_cpp[i]:.4f}")
