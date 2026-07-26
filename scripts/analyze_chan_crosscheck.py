import numpy as np
import os

output_dir = "/workspace/nr-link-simulator/results/cross_check_chan/"

def load_npy_complex(filename):
    return np.load(os.path.join(output_dir, filename))

tx_grid = load_npy_complex("tx_grid.npy")
rx_grid = load_npy_complex("rx_grid.npy")
h_est = load_npy_complex("h_est.npy")
rx_pdsch = load_npy_complex("rx_pdsch.npy")
h_at_pdsch = load_npy_complex("h_at_pdsch.npy")
equalized = load_npy_complex("equalized.npy")
eff_noise = np.load(os.path.join(output_dir, "eff_noise_var.npy"))
llr = np.load(os.path.join(output_dir, "llr_values.npy"))

with open(os.path.join(output_dir, "config.txt")) as f:
    config = {}
    for line in f:
        k, v = line.strip().split("=")
        try:
            config[k] = float(v) if "." in v else int(v)
        except:
            config[k] = v

n_ant, n_sym, n_sc = tx_grid.shape
noise_var = config['noise_var']
n_pdsch_re = int(config['n_pdsch_re'])

print("=== DMRS positions ===")
for sym in range(n_sym):
    dmrs_sc = []
    for sc in range(n_sc):
        if np.abs(tx_grid[0, sym, sc]) > 1e-10:
            # Check if it's data or DMRS: on sym=2 check SC positions
            dmrs_sc.append(sc)
    print(f"  sym={sym}: {len(dmrs_sc)} non-zero REs, first few: {dmrs_sc[:6]}")

print("\n=== h_est stats at DMRS sym=2 ===")
h_dmrs = h_est[0, 2, :]
print(f"  Mean h_est at DMRS sym: {np.mean(h_dmrs):.4f}")
print(f"  Mean |h_est|: {np.mean(np.abs(h_dmrs)):.4f}")
print(f"  MSE vs ideal h=1: {np.mean(np.abs(h_dmrs - 1)**2):.6f}")
print(f"  First 6 h_est values: {h_dmrs[:6]}")

print(f"\n=== h_est at data symbols (all should be same as sym=2 in AWGN) ===")
for sym in [0, 1, 2, 3, 12, 13]:
    h = h_est[0, sym, :]
    max_diff = np.max(np.abs(h - h_dmrs))
    print(f"  sym={sym}: mean|h|={np.mean(np.abs(h)):.4f}, MSE(1)={np.mean(np.abs(h-1)**2):.6f}, max|diff from sym2|={max_diff:.6f}")

print(f"\n=== h_at_pdsch shape ===")
print(f"  h_at_pdsch shape: {h_at_pdsch.shape}")
print(f"  rx_pdsch shape: {rx_pdsch.shape}")
print(f"  equalized shape: {equalized.shape}")
print(f"  n_pdsch_re = {n_pdsch_re}")

# Extract tx data symbols properly
# Per is_pdsch_symbol: DMRS sym=2 is skipped entirely, all other syms carry data
tx_data = []
h_est_data = []
for sym in range(n_sym):
    if sym == 2:
        continue  # DMRS symbol - no PDSCH
    for sc in range(n_sc):
        tx_data.append(tx_grid[0, sym, sc])
        h_est_data.append(h_est[0, sym, sc])
tx_data = np.array(tx_data)
h_est_data = np.array(h_est_data)
print(f"\n  Extracted {len(tx_data)} tx symbols (expected {n_pdsch_re})")

if len(tx_data) == n_pdsch_re:
    eq = equalized[:, 0]
    print(f"  Mean |h_est| at data REs: {np.mean(np.abs(h_est_data)):.4f}")
    print(f"  MSE of h_est vs 1: {np.mean(np.abs(h_est_data - 1)**2):.6f}")
    
    # Check equalization
    error = eq - tx_data
    mse_eq = np.mean(np.abs(error)**2)
    print(f"  MSE(eq - tx): {mse_eq:.6f} (noise_var={noise_var:.6f})")
    print(f"  Mean eff_noise_var: {np.mean(eff_noise):.6f}")
    
    # Check if eq has correct scaling
    # E[eq * conj(tx)] should be ~1 for unit gain
    corr = np.mean(eq * np.conj(tx_data))
    print(f"  E[eq * conj(tx)]: {corr:.4f} (should be ~1)")
    
    print(f"\n  First 12 constellation points:")
    for i in range(12):
        print(f"    {i}: tx={tx_data[i]:.4f}, eq={eq[i]:.4f}, err={error[i]:.4f}, eff_nv={eff_noise[i]:.4f}")

print(f"\n=== LLR stats ===")
print(f"  LLR mean: {np.mean(llr):.4f}")
print(f"  LLR std: {np.std(llr):.4f}")
print(f"  LLR min/max: {np.min(llr):.2f} / {np.max(llr):.2f}")
