import torch
import numpy as np
torch.manual_seed(42); np.random.seed(42)

qm=2; n_prb=30; n_sc=360; n_ofdm_sym=14
n_layers=4; n_rx=4; n_tx=4
pdsch_symbols=[s for s in range(n_ofdm_sym) if s not in [2,7]]
n_pdsch = len(pdsch_symbols)
M = n_layers * n_pdsch * n_sc
print(f"M={M}")

syms = torch.arange(M, dtype=torch.float32).to(torch.complex64).unsqueeze(0)
print(f"syms[0,:10] = {syms[0,:10]}")

syms_grid = syms.reshape(1, n_layers, n_sc, n_pdsch).permute(0,1,3,2)
print(f"After reshape sc-first then permute: syms_grid[0,0,0,:10] = {syms_grid[0,0,0,:10]}")
x_hat = syms_grid.permute(0,2,3,1).reshape(1, n_pdsch*n_sc, n_layers)
rx = x_hat.reshape(1, n_pdsch, n_sc, n_layers).permute(0,3,1,2).reshape(1,-1)
print(f"rx[0,:10] = {rx[0,:10]}")
print(f"match: {torch.equal(syms, rx)}")

print("\nDirect reshape to (n_layers, n_pdsch, n_sc) [layer, sym, sc]:")
syms_grid2 = syms.reshape(1, n_layers, n_pdsch, n_sc)
print(f"syms_grid2[0,0,0,:10] = {syms_grid2[0,0,0,:10]}")
rx2 = syms_grid2.permute(0,2,3,1).reshape(1, n_pdsch*n_sc, n_layers)
rx2 = rx2.reshape(1, n_pdsch, n_sc, n_layers).permute(0,3,1,2).reshape(1,-1)
print(f"rx2[0,:10] = {rx2[0,:10]}")
print(f"match: {torch.equal(syms, rx2)}")
