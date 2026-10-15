# PDSCH 基线数据（修改前参考版本）

> 用途：代码修复（速率匹配 k0、LDPC parity interleaver、DMRS comb-2 等）后的对照基线。
> 采集时间：2026-09-19，代码版本：main 分支（修复前）。

## 仿真配置

| 参数 | 值 |
|---|---|
| 信道 | AWGN（频域快路径），理想 CSI（perfect CSI） |
| 天线 | 1x1 SISO，单位阵预编码 |
| 带宽 | 3 PRB，156 RE/PRB，SCS=15 kHz（slot=1 ms） |
| LDPC | 25 次迭代，Sum-Product，early termination |
| 统计 | 每 SNR 点 target=60 错误块 / 上限 2000 块；0 错误 50 块即停 |
| MCS | 0,2,4,6,9,11,13,16,19,21,24,26,27（Table 1，覆盖 QPSK/16QAM/64QAM） |
| 种子 | 123（每个 SNR 点固定起始种子） |

## 产物

- `baseline_bler_awgn.csv` — 原始数据：`mcs,qm,rate,tbs,snr_db,bler,n_blocks,n_errors`
- `baseline_snr_vs_bler.png` — SNR vs BLER（对数轴）
- `baseline_snr_vs_tput.png` — SNR vs 吞吐量（`Tput = TBS·(1−BLER)/1ms`，Mbps）
- 复现命令：在 WSL 中 `build_wsl/examples/pdsch_bler_baseline`（源码 `examples/pdsch_bler_baseline.cpp`）
- 绘图：`python3 scripts/plot_baseline_curves.py`

## 关键参考点（10% BLER 门限附近）

| MCS | 调制 | 码率 | TBS(bit) | 10% BLER SNR (dB) | 峰值 Tput (Mbps) |
|---|---|---|---|---|---|
| 0 | QPSK | 0.117 | 104 | ≈ −3.0 | 0.094 |
| 2 | QPSK | 0.188 | 176 | ≈ −1.0 | 0.176 |
| 4 | QPSK | 0.301 | 288 | ≈ 0.2 | 0.288 |
| 6 | QPSK | 0.438 | 408 | ≈ 2.2 | 0.408 |
| 9 | QPSK | 0.663 | 640 | ≈ 3.3 | 0.640 |
| 11 | 16QAM | 0.369 | 704 | ≈ 4.9 | 0.699 |
| 13 | 16QAM | 0.479 | 888 | ≈ 6.2 | 0.888 |
| 16 | 16QAM | 0.643 | 1224 | ≈ 9.1 | 1.140 |
| 19 | 64QAM | 0.505 | 1416 | ≈ 10.9 | 1.394 |
| 21 | 64QAM | 0.602 | 1736 | ≈ 12.9 | 1.566 |
| 24 | 64QAM | 0.754 | 2152 | ≈ 15.5 | 2.152 |
| 26 | 64QAM | 0.853 | 2408 | ≈ 17.5 | 2.408 |
| 27 | 64QAM | 0.889 | 2472 | ≈ 18.3 | 2.472 |

## 备注

- BLER≈0 的点最多统计 50 块，曲线下端仅表示"无错误"；对比修改效果时以瀑布区（10%~90% BLER）为准。
- 修改速率匹配 k0 / parity interleaver / DMRS 后，预期瀑布区 SNR 会移动（k0 修正主要影响 RV>0 的场景；对 RV=0 的本基线影响可能较小）。若曲线整体移位明显，应检查是否引入回归。
