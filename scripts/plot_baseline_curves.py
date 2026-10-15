#!/usr/bin/env python3
"""Plot baseline link curves (SNR vs BLER, SNR vs Throughput) for the
pre-fix PDSCH simulation. Reads results/baseline/baseline_bler_awgn.csv
and writes results/baseline/baseline_snr_vs_bler.png and
results/baseline/baseline_snr_vs_tput.png."""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
CSV = os.path.join(ROOT, "results", "baseline", "baseline_bler_awgn.csv")
OUTDIR = os.path.join(ROOT, "results", "baseline")

df = pd.read_csv(CSV)
# slot duration at SCS=15kHz with 14 symbols = 1 ms
df["tput_mbps"] = df["tbs"] * (1.0 - df["bler"]) / 1.0e-3 / 1.0e6

mcs_list = sorted(df["mcs"].unique())
cmap = plt.get_cmap("viridis")
colors = {m: cmap(i / max(1, len(mcs_list) - 1)) for i, m in enumerate(mcs_list)}

mod_name = {2: "QPSK", 4: "16QAM", 6: "64QAM", 8: "256QAM"}

# ---- SNR vs BLER ----
fig, ax = plt.subplots(figsize=(10, 6.5))
for m in mcs_list:
    d = df[df["mcs"] == m].sort_values("snr_db")
    qm = d["qm"].iloc[0]
    r = d["rate"].iloc[0]
    ax.semilogy(d["snr_db"], d["bler"].clip(lower=1.0 / (d["n_blocks"].max() * 2)),
                "o-", ms=4, lw=1.4, color=colors[m],
                label=f"MCS {m} ({mod_name[qm]},R={r:.2f})")
ax.set_xlabel("SNR (dB)")
ax.set_ylabel("BLER")
ax.set_title("PDSCH Baseline (pre-fix): SNR vs BLER\n"
             "3 PRB, AWGN, perfect CSI, SISO, slot=1ms@15kHz")
ax.grid(True, which="both", alpha=0.4)
ax.set_ylim(1e-3, 1.5)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "baseline_snr_vs_bler.png"), dpi=150)

# ---- SNR vs Throughput ----
fig, ax = plt.subplots(figsize=(10, 6.5))
for m in mcs_list:
    d = df[df["mcs"] == m].sort_values("snr_db")
    qm = d["qm"].iloc[0]
    r = d["rate"].iloc[0]
    ax.plot(d["snr_db"], d["tput_mbps"], "o-", ms=4, lw=1.4, color=colors[m],
            label=f"MCS {m} ({mod_name[qm]},R={r:.2f})")
ax.set_xlabel("SNR (dB)")
ax.set_ylabel("Throughput (Mbps)")
ax.set_title("PDSCH Baseline (pre-fix): SNR vs Throughput\n"
             "3 PRB x 156 RE, TBS per 38.214, BLER-weighted, slot=1ms")
ax.grid(True, alpha=0.4)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "baseline_snr_vs_tput.png"), dpi=150)

# peak throughput summary
peak = df.groupby("mcs").apply(
    lambda d: d.loc[d["tput_mbps"].idxmax(), ["snr_db", "tput_mbps", "tbs"]])
print("=== Peak throughput per MCS (10% BLER reference) ===")
b10 = df[df["bler"] <= 0.10].groupby("mcs").first()
print(b10[["snr_db", "bler", "tbs", "tput_mbps"]].to_string())
print("Saved:", os.path.join(OUTDIR, "baseline_snr_vs_bler.png"))
print("Saved:", os.path.join(OUTDIR, "baseline_snr_vs_tput.png"))
