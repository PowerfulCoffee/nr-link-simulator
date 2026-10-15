#!/usr/bin/env python3
"""Compare pre-fix vs post-fix PDSCH link curves (AWGN, 3 PRB, SISO).

Loads the 4 sweep CSVs (perfect/LS CSI x prefix/postfix), computes the SNR at
10% BLER per MCS via linear interpolation in the waterfall region, plots
overlaid SNR-vs-BLER and SNR-vs-Throughput curves, and writes a delta table.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

BASE = os.path.join(os.path.dirname(__file__), "..", "results", "baseline")

def load(name):
    df = pd.read_csv(os.path.join(BASE, name))
    # throughput in Mbps: TBS bits per 1 ms slot
    df["tput"] = df["tbs"] * (1.0 - df["bler"]) / 1e-3 / 1e6
    return df

CURVES = {
    ("perfect", "pre"):  load("bler_prefix_perfect.csv"),
    ("perfect", "post"): load("bler_postfix_perfect.csv"),
    ("ls", "pre"):       load("bler_prefix_ls.csv"),
    ("ls", "post"):      load("bler_postfix_ls.csv"),
}

def snr_at_bler(df, mcs, target=0.1):
    d = df[df["mcs"] == mcs].sort_values("snr_db")
    snr, bler = d["snr_db"].values, d["bler"].values
    # find first crossing below target
    for i in range(len(bler) - 1):
        if bler[i] >= target > bler[i + 1]:
            # linear interpolation in log(bler)
            lb0, lb1 = np.log10(max(bler[i], 1e-4)), np.log10(max(bler[i + 1], 1e-4))
            frac = (np.log10(target) - lb0) / (lb1 - lb0) if lb1 != lb0 else 0.0
            return snr[i] + frac * (snr[i + 1] - snr[i])
    return np.nan

mcs_list = sorted(CURVES[("perfect", "pre")]["mcs"].unique())

rows = []
for mcs in mcs_list:
    qm = CURVES[("perfect", "pre")][CURVES[("perfect", "pre")]["mcs"] == mcs]["qm"].iloc[0]
    row = {"mcs": mcs, "qm": qm}
    for mode in ("perfect", "ls"):
        s_pre = snr_at_bler(CURVES[(mode, "pre")], mcs)
        s_post = snr_at_bler(CURVES[(mode, "post")], mcs)
        row[f"{mode}_pre"] = s_pre
        row[f"{mode}_post"] = s_post
        row[f"{mode}_delta"] = s_post - s_pre if (np.isfinite(s_pre) and np.isfinite(s_post)) else np.nan
    rows.append(row)

tab = pd.DataFrame(rows)
print("=== SNR @ 10% BLER (dB): post - pre (negative = improvement) ===")
print(tab.to_string(index=False, float_format=lambda v: f"{v:6.2f}"))
tab.to_csv(os.path.join(BASE, "fix_impact_snr10_table.csv"), index=False)

# --- Plots ---
for mode, fname in (("perfect", "compare_perfect_csi.png"), ("ls", "compare_ls_csi.png")):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for mcs in mcs_list:
        for ver, ls, alpha in (("pre", "--", 0.55), ("post", "-", 1.0)):
            d = CURVES[(mode, ver)]
            dd = d[d["mcs"] == mcs]
            lbl = f"MCS{mcs} ({ver})"
            axes[0].semilogy(dd["snr_db"], np.maximum(dd["bler"], 1e-4), ls,
                             alpha=alpha, label=lbl if mode == "perfect" else lbl, lw=1.4)
            axes[1].plot(dd["snr_db"], dd["tput"], ls, alpha=alpha, lw=1.4)
    axes[0].set_xlabel("SNR (dB)"); axes[0].set_ylabel("BLER")
    axes[0].set_title(f"SNR vs BLER ({'perfect CSI' if mode == 'perfect' else 'LS CSI'})")
    axes[0].set_ylim(1e-4, 1.2); axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=6, ncol=2)
    axes[1].set_xlabel("SNR (dB)"); axes[1].set_ylabel("Throughput (Mbps)")
    axes[1].set_title(f"SNR vs Throughput ({'perfect CSI' if mode == 'perfect' else 'LS CSI'})")
    axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(BASE, fname)
    plt.savefig(out, dpi=130)
    plt.close()
    print("Saved:", out)

# --- Summary stats ---
print("\n=== Summary (10% BLER SNR delta, post - pre) ===")
for mode in ("perfect", "ls"):
    d = tab[f"{mode}_delta"].dropna()
    print(f"{mode:8s}: mean {d.mean():+.3f} dB, min {d.min():+.3f}, max {d.max():+.3f}, "
          f"|delta|>0.5dB count {(d.abs() > 0.5).sum()}/{len(d)}")
