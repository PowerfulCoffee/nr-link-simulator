#!/usr/bin/env python3
"""Per-MCS grid plots: each MCS gets its own subplot with SNR-vs-BLER (log,
left axis) and SNR-vs-Throughput (linear, right axis), pre-fix vs post-fix
overlaid. One figure per CSI mode (perfect / LS)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

BASE = os.path.join(os.path.dirname(__file__), "..", "results", "baseline")

def load(name):
    df = pd.read_csv(os.path.join(BASE, name))
    df["tput"] = df["tbs"] * (1.0 - df["bler"]) / 1e-3 / 1e6  # Mbps per 1 ms slot
    return df

CURVES = {
    ("perfect", "pre"):  load("bler_prefix_perfect.csv"),
    ("perfect", "post"): load("bler_postfix_perfect.csv"),
    ("ls", "pre"):       load("bler_prefix_ls.csv"),
    ("ls", "post"):      load("bler_postfix_ls.csv"),
}

MODE_TITLE = {"perfect": "Perfect CSI", "ls": "LS CSI (non-ideal)"}
NCOLS, NROWS = 4, 4

for mode in ("perfect", "ls"):
    mcs_list = sorted(CURVES[(mode, "pre")]["mcs"].unique())
    fig, axes = plt.subplots(NROWS, NCOLS, figsize=(18, 12))
    axes = axes.ravel()

    for idx, mcs in enumerate(mcs_list):
        ax_bler = axes[idx]
        ax_tput = ax_bler.twinx()
        color = f"C{idx % 10}"

        for ver, ls, alpha in (("pre", "--", 0.6), ("post", "-", 1.0)):
            d = CURVES[(mode, ver)]
            dd = d[d["mcs"] == mcs].sort_values("snr_db")
            ax_bler.semilogy(dd["snr_db"], np.maximum(dd["bler"], 1e-4), ls,
                             color=color, alpha=alpha, lw=1.6,
                             label="pre-fix" if ver == "pre" else "post-fix")
            ax_tput.plot(dd["snr_db"], dd["tput"], ls, color="gray", alpha=alpha, lw=1.2)

        d0 = CURVES[(mode, "pre")]
        dd = d0[d0["mcs"] == mcs]
        qm = dd["qm"].iloc[0]
        tbs = dd["tbs"].iloc[0]
        ax_bler.set_title(f"MCS {mcs}  (Qm={qm}, TBS={tbs})", fontsize=10)
        ax_bler.set_ylim(1e-4, 1.2)
        ax_bler.set_xlabel("SNR (dB)", fontsize=8)
        ax_bler.set_ylabel("BLER", fontsize=8, color=color)
        ax_tput.set_ylabel("Tput (Mbps)", fontsize=8, color="gray")
        ax_tput.tick_params(axis="y", labelsize=7, colors="gray")
        ax_bler.tick_params(labelsize=7)
        ax_bler.grid(True, which="both", alpha=0.25)
        ax_bler.axhline(0.1, color="red", lw=0.6, alpha=0.5)

    # one legend on the first subplot
    axes[0].legend(fontsize=8, loc="lower left")
    for idx in range(len(mcs_list), NROWS * NCOLS):
        axes[idx].axis("off")

    fig.suptitle(f"PDSCH Link Curves per MCS - {MODE_TITLE[mode]} "
                 "(AWGN, 3 PRB, SISO; solid=post-fix, dashed=pre-fix; gray=Tput, right axis)",
                 fontsize=13)
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(BASE, f"per_mcs_grid_{mode}.png")
    plt.savefig(out, dpi=130)
    plt.close()
    print("Saved:", out)
